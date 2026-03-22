#!/usr/bin/env Rscript
## =============================================================================
## 02_brt.R  —  BRT 建模
## =============================================================================
## 策略（同疾病模型）:
##   循环内（100次）: 只在样点数据上做 OOF 预测 → AUC/TSS/重要性（瞬间完成）
##   循环外（10次）:  每 repeat 训练全量模型 → 1 次全国栅格预测
##   最后: 10 张预测图 → Welford 求 mean / sd
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dismo)
  library(gbm)
  library(dplyr)
})

# ---- 加载配置 ---------------------------------------------------------------
if (!exists("BASE_DIR")) {
  script_path <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) {
      args0 <- commandArgs(trailingOnly = FALSE)
      f <- grep("--file=", args0, value = TRUE)
      if (length(f) > 0) dirname(sub("--file=", "", f[1])) else getwd()
    }
  )
  source(file.path(script_path, "00_config.R"))
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}
buffer_km <- as.numeric(get_arg("--buffer_km", NA))
if (!is.finite(buffer_km)) buffer_km <- BUFFER_KM

OUTDIR <- buffer_out_dir(buffer_km)
tmpdir <- file.path(DIR_OUT, "tmp", paste0("brt_b", buffer_km))
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmpdir, memmax = 0.5, progress = 0)

cat("══════════════════════════════════════════════════════════\n")
cat("[BRT] BRT 建模 — M =", buffer_km, "km\n")
cat("══════════════════════════════════════════════════════════\n")

# ---- 加载数据 ---------------------------------------------------------------
df_model  <- readRDS(file.path(OUTDIR, "df_model.rds"))
fold_list <- readRDS(file.path(OUTDIR, "fold_list.rds"))
vars_keep <- readRDS(file.path(OUTDIR, "selected_vars.rds"))
env_files <- readRDS(file.path(OUTDIR, "env_files_selected.rds"))

cat("[BRT] 数据:", nrow(df_model), "行 |", length(vars_keep), "个变量 |",
    FINAL_REPEATS, "rep ×", FINAL_FOLDS, "fold\n")

# ---- 辅助函数 ---------------------------------------------------------------
max_tss_threshold <- function(pred, obs) {
  thr <- sort(unique(pred))
  thr <- thr[thr > 0 & thr < 1]
  if (length(thr) < 10) thr <- seq(0.01, 0.99, by = 0.01)
  best <- list(threshold = 0.5, TSS = -Inf, sens = NA, spec = NA)
  for (t in thr) {
    tp <- sum(pred >= t & obs == 1); tn <- sum(pred < t & obs == 0)
    fp <- sum(pred >= t & obs == 0); fn <- sum(pred < t & obs == 1)
    sens <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
    spec <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
    tss <- sens + spec - 1
    if (is.finite(tss) && tss > best$TSS)
      best <- list(threshold = t, TSS = tss, sens = sens, spec = spec)
  }
  best
}

# =============================================================================
# 阶段 1: CV 循环（只在样点上运算，不做栅格预测）
# =============================================================================
cat("\n[BRT] ── 阶段 1: 10×10 空间 CV（样点级）──\n")
t0 <- Sys.time()

oof_list  <- list()
perf_list <- list()
imp_list  <- list()
run_id <- 0L

for (r in 1:FINAL_REPEATS) {
  folds <- fold_list[[r]]
  
  for (k in 1:FINAL_FOLDS) {
    run_id <- run_id + 1L
    
    train_idx <- folds[[k]][[1]]
    test_idx  <- folds[[k]][[2]]
    train_df  <- df_model[train_idx, , drop = FALSE]
    test_df   <- df_model[test_idx, , drop = FALSE]
    
    brt_data <- train_df[, c("y", vars_keep), drop = FALSE]
    brt <- tryCatch(
      gbm.step(
        data = brt_data, gbm.x = 2:ncol(brt_data), gbm.y = 1,
        family = "bernoulli", tree.complexity = BRT_TC,
        learning.rate = BRT_LR, bag.fraction = BRT_BF,
        max.trees = BRT_MAXT, silent = TRUE, plot.main = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(brt)) {
      cat(sprintf("[BRT] ⚠ rep=%d fold=%d: 失败，跳过\n", r, k))
      next
    }
    
    # OOF（只在样点上预测，几千行数据瞬间完成）
    pred_test <- predict(brt, test_df[, c("y", vars_keep), drop = FALSE],
                         n.trees = brt$gbm.call$best.trees, type = "response")
    
    oof_list[[run_id]] <- data.frame(
      row_id = test_df$row_id, x = test_df$x, y_coord = test_df$y_coord,
      y = test_df$y, brt = as.numeric(pred_test)
    )
    
    # 性能（保护：如果测试折中只有一个类，跳过评估）
    n_pos <- sum(test_df$y == 1); n_neg <- sum(test_df$y == 0)
    if (n_pos == 0 || n_neg == 0) {
      cat(sprintf("[BRT] rep=%d fold=%d ⚠ 测试折仅含 %d/%d (pos/neg)，跳过评估\n",
                  r, k, n_pos, n_neg))
      rm(brt); gc(); next
    }
    
    ev  <- dismo::evaluate(p = pred_test[test_df$y == 1], a = pred_test[test_df$y == 0])
    tss <- max_tss_threshold(pred_test, test_df$y)
    perf_list[[run_id]] <- data.frame(
      Repeat = r, Fold = k, Run = run_id,
      AUC = ev@auc, TSS = tss$TSS,
      Threshold = tss$threshold, Sens = tss$sens, Spec = tss$spec
    )
    
    imp_list[[run_id]] <- brt$contributions %>% mutate(Repeat = r, Fold = k, Run = run_id)
    
    cat(sprintf("[BRT] rep=%d fold=%d AUC=%.4f TSS=%.4f\n", r, k, ev@auc, tss$TSS))
    
    rm(brt); gc()
  }
}

t_cv <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat("[BRT] 阶段 1 完成！耗时:", t_cv, "min\n")

# 保存 CV 结果
dir.create(file.path(OUTDIR, "tables"), recursive = TRUE, showWarnings = FALSE)
oof_df  <- bind_rows(oof_list) %>% arrange(row_id)
perf_df <- bind_rows(perf_list)
imp_df  <- bind_rows(imp_list)

write.csv(oof_df,  file.path(OUTDIR, "tables", "brt_oof.csv"), row.names = FALSE)
write.csv(perf_df, file.path(OUTDIR, "tables", "brt_performance.csv"), row.names = FALSE)
write.csv(imp_df,  file.path(OUTDIR, "tables", "brt_importance_allruns.csv"), row.names = FALSE)

# =============================================================================
# 阶段 2: 每 repeat 训全量模型 → 测试集栅格预测（10 次）
# =============================================================================
cat("\n[BRT] ── 阶段 2: 10 个全量模型 + 测试集栅格预测 ──\n")
t1 <- Sys.time()

# 加载测试集栅格
env_test_file <- file.path(OUTDIR, "env_files_test.rds")
if (!file.exists(env_test_file)) stop("[BRT] 未找到 env_files_test.rds，请检查 01_data_prep")
env_test_files <- readRDS(env_test_file)
env_test <- rast(env_test_files)
names(env_test) <- vars_keep
cat("[BRT] ✓ 测试集栅格已加载（2020 年环境）\n")

brt_data_full <- df_model[, c("y", vars_keep), drop = FALSE]
brt_models_for_pdp <- list()

# Welford: 测试集栅格
n_rast <- 0L; mean_all <- NULL; M2_all <- NULL

for (r in 1:FINAL_REPEATS) {
  cat(sprintf("[BRT] 全量模型 %d/%d ...", r, FINAL_REPEATS))
  
  set.seed(3000 + r * 13 + buffer_km)
  
  brt_full <- tryCatch(
    gbm.step(
      data = brt_data_full, gbm.x = 2:ncol(brt_data_full), gbm.y = 1,
      family = "bernoulli", tree.complexity = BRT_TC,
      learning.rate = BRT_LR, bag.fraction = BRT_BF,
      max.trees = BRT_MAXT, silent = TRUE, plot.main = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(brt_full)) {
    cat(" 失败，跳过\n"); next
  }
  
  # 测试集栅格预测
  pred_r <- predict(env_test, brt_full, n.trees = brt_full$gbm.call$best.trees,
                    type = "response")
  n_rast <- n_rast + 1L
  if (n_rast == 1L) {
    mean_all <- pred_r; M2_all <- pred_r * 0
  } else {
    dx <- pred_r - mean_all
    mean_all <- mean_all + dx / n_rast
    M2_all   <- M2_all + dx * (pred_r - mean_all)
  }
  
  brt_models_for_pdp[[r]] <- brt_full
  
  cat(sprintf(" 完成 (%d 棵树)\n", brt_full$gbm.call$best.trees))
  rm(brt_full, pred_r); gc()
}

t_rast <- round(difftime(Sys.time(), t1, units = "mins"), 1)
cat("[BRT] 阶段 2 完成！耗时:", t_rast, "min\n")

if (n_rast < 2) stop("[BRT] 有效全量模型不足")

# 保存栅格（在测试集环境上的预测）
sd_brt <- sqrt(M2_all / (n_rast - 1))
writeRaster(mean_all, file.path(OUTDIR, "BRT_mean.tif"), overwrite = TRUE)
writeRaster(sd_brt,   file.path(OUTDIR, "BRT_sd.tif"),   overwrite = TRUE)

# 时间外推 AUC（训练集样点坐标 × 测试集栅格预测值）
pts_xy <- df_model[, c("x", "y_coord")]
test_vals <- terra::extract(mean_all, pts_xy)[, 2]
valid <- is.finite(test_vals) & is.finite(df_model$y)
if (sum(df_model$y[valid] == 1) > 0 && sum(df_model$y[valid] == 0) > 0) {
  ev_test <- dismo::evaluate(
    p = test_vals[valid & df_model$y == 1],
    a = test_vals[valid & df_model$y == 0]
  )
  cat(sprintf("[BRT] 测试集 AUC (时间外推): %.4f\n", ev_test@auc))
  
  # 训练集 AUC（来自阶段1 CV）vs 测试集 AUC → 差值
  train_auc <- mean(perf_df$AUC, na.rm = TRUE)
  test_auc  <- ev_test@auc
  delta_auc <- train_auc - test_auc
  
  auc_compare <- data.frame(
    model = "BRT",
    AUC_train_cv = round(train_auc, 4),
    AUC_test_temporal = round(test_auc, 4),
    delta_AUC = round(delta_auc, 4),
    robust = ifelse(abs(delta_auc) < 0.05, "YES", "NO"),
    n_models = n_rast
  )
  write.csv(auc_compare, file.path(OUTDIR, "tables", "brt_test_auc.csv"), row.names = FALSE)
  
  cat(sprintf("[BRT] 训练集 CV AUC: %.4f | 测试集 AUC: %.4f | ΔAUC: %.4f %s\n",
              train_auc, test_auc, delta_auc,
              ifelse(abs(delta_auc) < 0.05, "✓ 稳健", "⚠ 差距较大")))
}

rm(mean_all, M2_all, sd_brt, env_test); gc()

# PDP 子模型
brt_models_for_pdp <- Filter(Negate(is.null), brt_models_for_pdp)
saveRDS(brt_models_for_pdp, file.path(OUTDIR, "brt_cv_models.rds"))

# 性能汇总
perf_summary <- perf_df %>%
  summarise(AUC_mean = mean(AUC, na.rm = TRUE), AUC_sd = sd(AUC, na.rm = TRUE),
            TSS_mean = mean(TSS, na.rm = TRUE), TSS_sd = sd(TSS, na.rm = TRUE),
            n_runs = n())
write.csv(perf_summary, file.path(OUTDIR, "tables", "brt_performance_summary.csv"),
          row.names = FALSE)

# RData
elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
save(oof_df, perf_df, imp_df, perf_summary,
     brt_models_for_pdp, n_rast, elapsed,
     file = file.path(OUTDIR, "brt_results.RData"))

cat("══════════════════════════════════════════════════════════\n")
cat("[BRT] 完成！CV:", t_cv, "min + 栅格:", t_rast, "min\n")
cat("[BRT] AUC:", round(perf_summary$AUC_mean, 4), "±", round(perf_summary$AUC_sd, 4), "\n")
cat("[BRT] TSS:", round(perf_summary$TSS_mean, 4), "±", round(perf_summary$TSS_sd, 4), "\n")
cat("[BRT] 栅格预测:", n_rast, "次 (全量模型)\n")
terra::tmpFiles(remove = TRUE)
cat("══════════════════════════════════════════════════════════\n")