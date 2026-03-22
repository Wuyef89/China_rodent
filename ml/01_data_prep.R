#!/usr/bin/env Rscript
## =============================================================================
## 01_data_prep.R  —  数据预处理（优化版：混合背景采样 + 空间块强化）
## =============================================================================
## 更新：引入 Hybrid Background (SRE + Random) 以降低模型过拟合风险
## 更新：强制 Block Size 最小值为 400km 以应对干旱区空间粘滞性
## =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(blockCV)
  library(dismo)
  library(gbm)
  library(gstat)
  library(sp)
  library(future.apply)
})

# ---- 加载配置 ---------------------------------------------------------------
if (!exists("BASE_DIR")) {
  script_path <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) {
      args <- commandArgs(trailingOnly = FALSE)
      f <- grep("--file=", args, value = TRUE)
      if (length(f) > 0) dirname(sub("--file=", "", f[1])) else getwd()
    }
  )
  config_file <- file.path(script_path, "00_config.R")
  if (file.exists(config_file)) source(config_file)
  else stop("未找到 00_config.R")
}

buffer_km <- BUFFER_KM
set.seed(20260211 + buffer_km)

# ---- 输出目录 ---------------------------------------------------------------
OUTDIR <- buffer_out_dir(buffer_km)
dir.create(file.path(OUTDIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "figures"), recursive = TRUE, showWarnings = FALSE)

tmpdir <- file.path(DIR_OUT, "tmp", paste0("prep_b", buffer_km))
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmpdir)

cat("══════════════════════════════════════════════════════════\n")
cat("[01] 数据预处理 — 缓冲区 M =", buffer_km, "km\n")
cat("══════════════════════════════════════════════════════════\n")

# =============================================================================
# 辅助函数
# =============================================================================

make_M_buffer <- function(pres_sf, buffer_m) {
  aggregate(buffer(vect(pres_sf), width = buffer_m))
}

# 优化后的混合采样逻辑
sample_hybrid_bg <- function(env_sre, env_model, pres_sf, M_v, n_target, 
                             oversample = 8, seed = 1L) {
  set.seed(seed)
  # 50/50 分配：一半环境极端点，一半纯空间随机点
  n_sre <- floor(n_target * 0.5)
  n_rnd <- n_target - n_sre
  
  cat(sprintf("[01]    混合采样策略: SRE(%d) + Random(%d)\n", n_sre, n_rnd))
  
  # 1. 采样 SRE 背景
  bg_sre_v <- sample_bg_SRE_core(env_sre, env_model, pres_sf, M_v, n_sre, oversample, seed)
  
  # 2. 采样 纯随机 背景
  pres_cell <- cellFromXY(env_model[[1]], crds(vect(pres_sf)))
  bg_rnd_v <- spatSample(env_model[[1]], size = n_rnd * 5, method = "random",
                         na.rm = TRUE, as.points = TRUE, ext = M_v)
  bg_rnd_cell <- cellFromXY(env_model[[1]], crds(bg_rnd_v))
  bg_rnd_v <- bg_rnd_v[!(bg_rnd_cell %in% pres_cell)] # 排除存在点所在像元
  
  if (nrow(bg_rnd_v) > n_rnd) bg_rnd_v <- bg_rnd_v[1:n_rnd]
  
  return(rbind(bg_sre_v, bg_rnd_v))
}

# 内部 SRE 计算逻辑
sample_bg_SRE_core <- function(env_sre, env_model, pres_sf, M_v, n_target, oversample, seed) {
  pres_v <- vect(pres_sf)
  pres_x <- terra::extract(env_sre, pres_v)
  if ("ID" %in% names(pres_x)) pres_x <- pres_x[, -1, drop = FALSE]
  lo <- apply(pres_x, 2, min, na.rm = TRUE)
  hi <- apply(pres_x, 2, max, na.rm = TRUE)
  pres_cell <- cellFromXY(env_model[[1]], crds(pres_v))
  
  n_cand <- n_target * oversample
  cand_v <- spatSample(env_model[[1]], size = n_cand, method = "random",
                       na.rm = TRUE, as.points = TRUE, ext = M_v)
  cand_cell <- cellFromXY(env_model[[1]], crds(cand_v))
  cand_v <- cand_v[!(cand_cell %in% pres_cell)]
  
  cand_x <- terra::extract(env_sre, cand_v)
  if ("ID" %in% names(cand_x)) cand_x <- cand_x[, -1, drop = FALSE]
  
  inside <- rep(TRUE, nrow(cand_x))
  for (vn in names(cand_x)) {
    inside <- inside & (cand_x[[vn]] >= lo[[vn]]) & (cand_x[[vn]] <= hi[[vn]])
  }
  bg_v <- cand_v[!inside]
  
  if (nrow(bg_v) < n_target) bg_v <- rbind(bg_v, cand_v[inside][1:min(n_target-nrow(bg_v), sum(inside))])
  bg_v[1:min(n_target, nrow(bg_v))]
}

# [其他辅助函数 get_range_from_res, range_by_variogram 保持不变...]
get_range_from_res <- function(res) {
  if (is.null(res)) return(NA_real_)
  for (f in c("range", "range_size", "median_range"))
    if (!is.null(res[[f]])) return(as.numeric(res[[f]]))
  if (!is.null(res$results) && "range" %in% names(res$results))
    return(as.numeric(res$results$range[1]))
  NA_real_
}

range_by_variogram <- function(sf_pts, value_col = "z") {
  sp_pts <- as(sf_pts, "Spatial")
  vgm_emp <- gstat::variogram(stats::as.formula(paste0(value_col, "~1")), sp_pts)
  psill0 <- stats::var(sf_pts[[value_col]], na.rm = TRUE)
  if (!is.finite(psill0) || psill0 <= 0) return(NA_real_)
  maxd <- max(vgm_emp$dist, na.rm = TRUE)
  model0 <- gstat::vgm(psill = psill0, model = "Sph", range = maxd / 3, nugget = psill0 * 0.1)
  fit <- try(gstat::fit.variogram(vgm_emp, model0), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  rg <- fit$range[which(fit$model != "Nug")]
  as.numeric(rg[1])
}

# =============================================================================
# 1) 加载协变量 + 存在点
# =============================================================================
cat("[01] 1. 加载协变量栅格 ...\n")
env_files <- list.files(DIR_ENV_TRAIN, pattern = "\\.tif$", full.names = TRUE)
env <- rast(env_files)
names(env) <- tools::file_path_sans_ext(basename(env_files))
all_var_names <- names(env)

cat("[01] 2. 提取存在点 ...\n")
pr <- rast(PRES_RASTER)
cells <- which(!is.na(values(pr, mat=F)) & values(pr, mat=F) == 1)
xy <- xyFromCell(pr, cells)
pres_sf <- st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = crs(pr))
pres_sf <- st_transform(pres_sf, crs = crs(env))

# =============================================================================
# 2) 构建缓冲区 M + 混合背景采样
# =============================================================================
cat("[01] 3. 缓冲区 M =", buffer_km, "km + Hybrid 伪缺失 ...\n")

buffer_m <- buffer_km * 1000
n_pres <- nrow(pres_sf)
n_bg   <- n_pres * PB_RATIO
M_v <- make_M_buffer(pres_sf, buffer_m)

sre_names <- intersect(SRE_VARS, all_var_names)
env_sre <- env[[sre_names]]

# 使用优化后的采样函数
bg_v <- sample_hybrid_bg(
  env_sre = env_sre, env_model = env, pres_sf = pres_sf, M_v = M_v,
  n_target = n_bg, oversample = SRE_OVERSAMPLE, seed = 9000 + buffer_km
)

bg_sf <- st_as_sf(bg_v) %>% st_transform(crs = crs(env))

# 合并与 1:1 平衡
pres_sf$y <- 1L; bg_sf$y <- 0L
pts_all <- bind_rows(pres_sf[, "y"], bg_sf[, "y"])
pts_all$row_id <- seq_len(nrow(pts_all))

xdat <- terra::extract(env, vect(pts_all))
if ("ID" %in% names(xdat)) xdat <- xdat[, -1, drop = FALSE]
df_all <- cbind(row_id = pts_all$row_id, y = pts_all$y, xdat) %>% filter(complete.cases(.))
pts_all <- pts_all[df_all$row_id, ]

set.seed(20260211)
pres_idx <- which(df_all$y == 1L); bg_idx <- which(df_all$y == 0L)
n_min <- min(length(pres_idx), length(bg_idx))
keep_idx <- sort(c(sample(pres_idx, n_min), sample(bg_idx, n_min)))

df0 <- df_all[keep_idx, ]
pts_sf <- pts_all[keep_idx, ]
pts_sf$row_id <- seq_len(nrow(pts_sf))
df0$row_id <- pts_sf$row_id

# =============================================================================
# 3) 空间自相关估计 → Block Size 强化 (Min 300km)
# =============================================================================
cat("[01] 4. 空间自相关估计 (目标: 强化空间隔离) ...\n")

# [变量分组与 range 估计逻辑维持原样...]
ECO_VARS <- c("Elevation",  "TEM", "GST", "PRE", "EVP", "RHU", "PRS", "SSD", "WIN", "Mammal_richness")
SOCIO_VARS <- c("GDP", "PD", "ReportEffort")

pts_env <- terra::extract(env, vect(pts_sf))
if ("ID" %in% names(pts_env)) pts_env <- pts_env[, -1, drop = FALSE]
for (vn in all_var_names) pts_sf[[vn]] <- pts_env[[vn]]

auto_ranges <- data.frame(var = character(), group = character(), range_m = double(), stringsAsFactors = FALSE)
for (vn in all_var_names) {
  df1 <- pts_sf[, vn, drop = FALSE]; names(df1)[1] <- "z"; df1 <- df1[!is.na(df1$z), ]
  if (nrow(df1) < 200) next
  rg <- range_by_variogram(df1, "z")
  auto_ranges <- rbind(auto_ranges, data.frame(var = vn, group = ifelse(vn %in% ECO_VARS, "Ecological", "Other"), range_m = rg))
}

med_range <- median(auto_ranges$range_m[auto_ranges$group == "Ecological"], na.rm = TRUE)
if (!is.finite(med_range)) med_range <- 200000

#修改后：
if (exists("BLOCK_FORCE_FIXED") && BLOCK_FORCE_FIXED) {
  # 既然是强制实验，就直接给值，不要再套 max(300000) 了
  BLOCK_SIZE_M <- BLOCK_CLAMP_MAX_M
  cat("[01]    ！！！实验模式：直接使用空间块大小 =", BLOCK_SIZE_M / 1000, "km\n")
} else {
  # 生产模式：自适应计算，且不小于 300km (为了稳健性)
  BLOCK_SIZE_M <- max(300000, min(BLOCK_CLAMP_MAX_M, med_range))
  cat("[01]    自适应模式：Block size =", round(BLOCK_SIZE_M / 1000), "km\n")
}

# =============================================================================
# 4) 变量筛选（3×5 空间 CV，在全部变量上做）
#    目的: 不是删除变量来避免共线性，
#    而是识别"完全无信号"的变量（贡献 < 2% 且不稳定）
# =============================================================================
cat("[01] 5. 变量重要性筛选（", VARSEL_REPEATS, "×", VARSEL_FOLDS, "）...\n")
cat("[01]   注意: 输入全部", length(all_var_names), "个变量（不做共线性筛选）\n")

# 并行
if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = min(3, future::availableCores()))
} else {
  future::plan(future::multicore, workers = min(3, future::availableCores()))
}

run_one_repeat <- function(r, buffer_km, pts_sf, df0,
                           BLOCK_SIZE_M, VARSEL_FOLDS,
                           BRT_TC, BRT_LR, BRT_BF, BRT_MAXT) {
  set.seed(1000 + r * 17 + buffer_km)
  
  sb <- blockCV::cv_spatial(
    x = pts_sf, column = "y",
    size = BLOCK_SIZE_M, k = VARSEL_FOLDS,
    selection = "random", iteration = 100, progress = FALSE
  )
  folds <- sb$folds_list
  
  contrib_list <- vector("list", VARSEL_FOLDS)
  meta_list    <- vector("list", VARSEL_FOLDS)
  
  for (k in seq_len(VARSEL_FOLDS)) {
    train_idx <- folds[[k]][[1]]
    test_idx  <- folds[[k]][[2]]
    train_df  <- df0[train_idx, , drop = FALSE]
    test_df   <- df0[test_idx, , drop = FALSE]
    
    pred_cols <- setdiff(names(train_df), c("row_id", "y"))
    brt_data  <- train_df[, c("y", pred_cols), drop = FALSE]
    
    mod <- try(gbm.step(
      data = brt_data, gbm.x = 2:ncol(brt_data), gbm.y = 1,
      family = "bernoulli", tree.complexity = BRT_TC,
      learning.rate = BRT_LR, bag.fraction = BRT_BF,
      max.trees = BRT_MAXT, silent = TRUE, plot.main = FALSE
    ), silent = TRUE)
    
    if (inherits(mod, "try-error") || is.null(mod)) next
    
    contrib_list[[k]] <- mod$contributions %>% mutate(Repeat = r, Fold = k)
    
    pred_test <- predict(mod, test_df, n.trees = mod$gbm.call$best.trees, type = "response")
    ev <- dismo::evaluate(p = pred_test[test_df$y == 1], a = pred_test[test_df$y == 0])
    
    meta_list[[k]] <- data.frame(
      Stage = "varsel", Repeat = r, Fold = k,
      n_train = nrow(train_df), n_test = nrow(test_df), AUC = ev@auc)
  }
  
  list(contrib = bind_rows(contrib_list), meta = bind_rows(meta_list))
}

res_list <- future.apply::future_lapply(
  X = seq_len(VARSEL_REPEATS),
  FUN = run_one_repeat,
  buffer_km = buffer_km, pts_sf = pts_sf, df0 = df0,
  BLOCK_SIZE_M = BLOCK_SIZE_M, VARSEL_FOLDS = VARSEL_FOLDS,
  BRT_TC = BRT_TC, BRT_LR = BRT_LR, BRT_BF = BRT_BF, BRT_MAXT = BRT_MAXT,
  future.seed = TRUE
)

contrib_df   <- bind_rows(lapply(res_list, `[[`, "contrib"))
fold_meta_df <- bind_rows(lapply(res_list, `[[`, "meta"))

# 稳定性汇总
stab <- contrib_df %>%
  group_by(var) %>%
  summarise(
    mean_imp     = mean(rel.inf, na.rm = TRUE),
    sd_imp       = sd(rel.inf, na.rm = TRUE),
    prop_nonzero = mean(rel.inf > 0, na.rm = TRUE),
    n_models     = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_imp))

# 筛选: 平均贡献 >= 2% 且 >= 60% 折叠中有非零贡献，两个条件都要满足
vars_keep <- stab %>%
  filter(mean_imp >= IMP_THRESH_MEAN, prop_nonzero >= KEEP_PROP) %>%
  pull(var)

# 如果筛选后不足 3 个，保留全部（树模型可以处理）
if (length(vars_keep) < 3) {
  cat("[01]   ⚠ 筛选后仅", length(vars_keep), "个变量，保留全部\n")
  vars_keep <- all_var_names
}

write.csv(contrib_df,   file.path(OUTDIR, "tables", "varsel_contrib.csv"), row.names = FALSE)
write.csv(fold_meta_df, file.path(OUTDIR, "tables", "varsel_fold_structure.csv"), row.names = FALSE)
write.csv(stab,         file.path(OUTDIR, "tables", "varsel_variable_stability.csv"), row.names = FALSE)

cat("[01]   输入:", length(all_var_names), "个 → 保留:", length(vars_keep), "个\n")
cat("[01]  ", paste(vars_keep, collapse = ", "), "\n")

# 被剔除的变量（如果有）
vars_dropped <- setdiff(all_var_names, vars_keep)
if (length(vars_dropped) > 0) {
  cat("[01]   剔除（真正无信号）:", paste(vars_dropped, collapse = ", "), "\n")
} else {
  cat("[01]   全部变量保留（没有变量被剔除）\n")
}

# 构建建模数据（含坐标）
coords_df <- as.data.frame(st_coordinates(pts_sf))
names(coords_df) <- c("x", "y_coord")
df_model <- cbind(
  df0[, c("row_id", "y"), drop = FALSE],
  x = coords_df$x, y_coord = coords_df$y_coord,
  df0[, vars_keep, drop = FALSE]
)

# =============================================================================
# 5) 生成 10×10 空间 CV 索引
# =============================================================================
cat("[01] 6. 生成", FINAL_REPEATS, "×", FINAL_FOLDS, "空间 CV 索引 ...\n")

fold_list <- list()
for (r in 1:FINAL_REPEATS) {
  set.seed(5000 + r * 31 + buffer_km)
  sb <- cv_spatial(
    x = pts_sf, column = "y",
    size = BLOCK_SIZE_M, k = FINAL_FOLDS,
    selection = "random", iteration = 100, progress = FALSE
  )
  fold_list[[r]] <- sb$folds_list
}
cat("[01]   完毕\n")

# =============================================================================
# 6) 保存
# =============================================================================
cat("[01] 7. 保存 ...\n")

saveRDS(df_model,     file.path(OUTDIR, "df_model.rds"))
saveRDS(pts_sf,       file.path(OUTDIR, "pts_sf.rds"))
saveRDS(vars_keep,    file.path(OUTDIR, "selected_vars.rds"))
saveRDS(fold_list,    file.path(OUTDIR, "fold_list.rds"))
saveRDS(BLOCK_SIZE_M, file.path(OUTDIR, "block_size_m.rds"))

env2_files <- file.path(DIR_ENV_TRAIN, paste0(vars_keep, ".tif"))
saveRDS(env2_files, file.path(OUTDIR, "env_files_selected.rds"))

# 测试集栅格路径（用于时间外推验证）
# 变量名→文件名映射（兼容 Dem_1km ↔ Elevation 等重命名情况）
var_to_filename <- function(vname, dir) {
  # 直接匹配
  f <- file.path(dir, paste0(vname, ".tif"))
  if (file.exists(f)) return(f)
  # 别名映射
  alias_map <- c("Dem_1km" = "Elevation", "Elevation" = "Dem_1km")
  if (vname %in% names(alias_map)) {
    f2 <- file.path(dir, paste0(alias_map[vname], ".tif"))
    if (file.exists(f2)) {
      cat("[01]   ℹ 变量", vname, "→ 文件", basename(f2), "\n")
      return(f2)
    }
  }
  return(f)  # 返回原路径（后面会检测不存在）
}

env_test_files <- sapply(vars_keep, var_to_filename, dir = DIR_ENV_TEST)
names(env_test_files) <- vars_keep
test_exist <- file.exists(env_test_files)
if (all(test_exist)) {
  saveRDS(env_test_files, file.path(OUTDIR, "env_files_test.rds"))
  cat("[01]   ✓ 测试集栅格:", length(env_test_files), "个变量均存在\n")
  cat("[01]   文件:", paste(basename(env_test_files), collapse = ", "), "\n")
} else {
  missing <- vars_keep[!test_exist]
  cat("[01]   ⚠ 测试集缺少变量:", paste(missing, collapse = ", "), "\n")
  cat("[01]   ⚠ 将跳过时间外推验证\n")
}

# 保存完整工作空间（以防后续复用）
save(df_model, df0, pts_sf, vars_keep, fold_list,
     BLOCK_SIZE_M, buffer_km, all_var_names, auto_ranges,
     stab, contrib_df, fold_meta_df,
     file = file.path(OUTDIR, "data_prep_results.RData"))
cat("[01]   ✓ data_prep_results.RData 已保存\n")

cat("══════════════════════════════════════════════════════════\n")
cat("[01] 完成！\n")
cat("[01]   变量: 全部", length(all_var_names), "→ 保留", length(vars_keep), "\n")
cat("[01]   Block:", round(BLOCK_SIZE_M / 1000), "km\n")
cat("[01]   输出:", OUTDIR, "\n")
cat("══════════════════════════════════════════════════════════\n")