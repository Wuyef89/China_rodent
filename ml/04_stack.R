#!/usr/bin/env Rscript
## =============================================================================
## 04_stack.R  —  Stacking 集成（地理坐标硬对齐 + QP & Logistic 双方法）
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(dismo)
  library(quadprog)
  library(readr)
  library(stringr)
})

source("/data/home/test01/rodent/script/00_config.R")

# ---- 解析参数 ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}
buffer_km <- as.numeric(get_arg("--buffer_km", NA))
if (!is.finite(buffer_km)) stop("用法: Rscript 04_stack.R --buffer_km 100|200|300")

OUTDIR <- buffer_out_dir(buffer_km)
TABDIR <- file.path(OUTDIR, "tables")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

cat("══════════════════════════════════════════════════════════\n")
cat("[STACK] Stacking 集成 — 缓冲区 M =", buffer_km, "km\n")
cat("══════════════════════════════════════════════════════════\n")

# =============================================================================
# 辅助函数
# =============================================================================

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
    tss  <- sens + spec - 1
    if (is.finite(tss) && tss > best$TSS)
      best <- list(threshold = t, TSS = tss, sens = sens, spec = spec)
  }
  best
}

calc_auc <- function(pred, obs) {
  ev <- tryCatch(dismo::evaluate(p = pred[obs == 1], a = pred[obs == 0]),
                 error = function(e) NULL)
  if (!is.null(ev)) ev@auc else NA_real_
}

calc_suitable_area_km2 <- function(binary_r) {
  a <- cellSize(binary_r, unit = "km")
  v <- values(binary_r, mat = FALSE)
  sum(values(a, mat = FALSE)[which(v == 1)], na.rm = TRUE)
}

read_block_info_safe <- function(outdir) {
  f_csv <- file.path(outdir, "tables", "block_size_choice.csv")
  f_rds <- file.path(outdir, "block_size_m.rds")

  block_km <- NA_real_
  method   <- NA_character_
  reason   <- NA_character_

  if (file.exists(f_csv)) {
    x <- tryCatch(readr::read_csv(f_csv, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(x) && nrow(x) > 0) {
      nm <- names(x)
      col_method <- intersect(nm, c("method","Method"))
      col_bkm    <- intersect(nm, c("block_size_km","block_km"))
      col_reason <- intersect(nm, c("reason","Reason"))
      if (length(col_bkm) >= 1) {
        block_km <- as.numeric(x[[col_bkm[1]]][1])
        method   <- if (length(col_method) >= 1) as.character(x[[col_method[1]]][1]) else "csv"
        reason   <- if (length(col_reason) >= 1) as.character(x[[col_reason[1]]][1]) else "provided"
      }
    }
  }

  if (!is.finite(block_km)) {
    if (!file.exists(f_rds)) stop("[STACK] No block_size_choice.csv and no block_size_m.rds found.")
    block_m <- readRDS(f_rds)
    block_km <- as.numeric(block_m) / 1000
    method <- "rds_fallback"
    reason <- "read from block_size_m.rds"
  }

  data.frame(method = method, block_size_km = block_km, reason = reason, stringsAsFactors = FALSE)
}

read_data_summary_safe <- function(outdir) {
  f <- file.path(outdir, "tables", "data_summary.csv")
  if (!file.exists(f)) {
    return(list(n_presence = NA_integer_, n_pseudoabsence = NA_integer_, n_total = NA_integer_))
  }
  x <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) {
    return(list(n_presence = NA_integer_, n_pseudoabsence = NA_integer_, n_total = NA_integer_))
  }

  nm <- names(x)

  # Case A: already wide columns
  if (all(c("n_presence","n_pseudoabsence") %in% nm)) {
    return(list(
      n_presence = as.integer(x$n_presence[1]),
      n_pseudoabsence = as.integer(x$n_pseudoabsence[1]),
      n_total = as.integer(if ("n_total" %in% nm) x$n_total[1] else NA)
    ))
  }

  # Case B: Metric/Value long format
  if (all(c("Metric","Value") %in% nm)) {
    getv <- function(key) {
      v <- x$Value[match(key, x$Metric)]
      if (length(v) == 0 || is.na(v)) NA_integer_ else as.integer(v)
    }
    n_total <- getv("Total_Samples")
    n_presence <- getv("n_presence")
    n_pseudoabsence <- getv("n_pseudoabsence")
    return(list(n_presence = n_presence, n_pseudoabsence = n_pseudoabsence, n_total = n_total))
  }

  list(n_presence = NA_integer_, n_pseudoabsence = NA_integer_, n_total = NA_integer_)
}

# =============================================================================
# 1) 加载 OOF 并以坐标为键严格对齐
# =============================================================================
cat("[STACK] 1. 加载 OOF 并地理坐标对齐 ...\n")

brt_oof <- read.csv(file.path(TABDIR, "brt_oof.csv"))
rf_oof  <- read.csv(file.path(TABDIR, "rf_oof.csv"))
cat(sprintf("  BRT OOF: %d 行 | RF OOF: %d 行\n", nrow(brt_oof), nrow(rf_oof)))

brt_avg <- brt_oof %>%
  group_by(row_id, x, y_coord) %>%
  summarise(y = first(y), brt = mean(brt, na.rm = TRUE),
            n_brt = n(), .groups = "drop")

rf_avg <- rf_oof %>%
  group_by(row_id, x, y_coord) %>%
  summarise(y = first(y), rf = mean(rf, na.rm = TRUE),
            n_rf = n(), .groups = "drop")

full_data <- inner_join(
  brt_avg %>% dplyr::select(row_id, x, y_coord, y, brt, n_brt),
  rf_avg  %>% dplyr::select(row_id, x, y_coord, rf, n_rf),
  by = c("row_id", "x", "y_coord")
) %>%
  filter(complete.cases(.)) %>%
  arrange(row_id)

cat(sprintf("  对齐后样本数: %d | 阳性: %d | 阴性: %d\n",
            nrow(full_data), sum(full_data$y == 1), sum(full_data$y == 0)))

write.csv(full_data, file.path(TABDIR, "oof_aligned.csv"), row.names = FALSE)

# =============================================================================
# 2) 方法 A: QP 非负约束权重 (w>=0, sum(w)=1)
# =============================================================================
cat("[STACK] 2. 方法A: QP 非负权重 ...\n")

X <- as.matrix(full_data[, c("brt", "rf")])
y_vec <- full_data$y
k <- ncol(X)

# min ||Xw - y||^2  s.t. w>=0, 1'w = 1
Dmat <- t(X) %*% X + diag(1e-8, k)
dvec <- t(X) %*% y_vec

# quadprog form: min 1/2 w'Dw - d'w
# constraints: A^T w >= b
# w>=0:
Amat1 <- diag(k)
bvec1 <- rep(0, k)
# sum(w)=1 as equality: a^T w = 1 -> implement with meq=1, include a in Amat
Amat  <- cbind(rep(1, k), Amat1)
bvec  <- c(1, bvec1)

qp_sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
w_qp <- qp_sol$solution
w_qp[w_qp < 0] <- 0
if (sum(w_qp) == 0) w_qp <- rep(1/k, k)
w_qp <- w_qp / sum(w_qp)

cat(sprintf("  QP 权重: wBRT = %.4f, wRF = %.4f\n", w_qp[1], w_qp[2]))

full_data$stack_qp <- w_qp[1] * full_data$brt + w_qp[2] * full_data$rf
auc_qp <- calc_auc(full_data$stack_qp, full_data$y)
tss_qp <- max_tss_threshold(full_data$stack_qp, full_data$y)

write.csv(data.frame(wBRT = w_qp[1], wRF = w_qp[2], AUC = auc_qp, TSS = tss_qp$TSS),
          file.path(TABDIR, "stacking_qp_weights.csv"), row.names = FALSE)

cat(sprintf("  QP Stacking: AUC = %.4f | TSS = %.4f\n", auc_qp, tss_qp$TSS))

# =============================================================================
# 3) 方法 B: Logistic 元学习器（Bootstrap 系数分布）
# =============================================================================
cat("[STACK] 3. 方法B: Logistic 元学习器 (100 次 Bootstrap) ...\n")

N_BOOT <- 100
set.seed(42 + buffer_km)

stacking_coefs <- data.frame()
std_coefs_all  <- data.frame()
boot_perf      <- data.frame()

for (b in 1:N_BOOT) {
  idx_1 <- sample(which(full_data$y == 1), replace = TRUE)
  idx_0 <- sample(which(full_data$y == 0), replace = TRUE)
  train_idx <- c(idx_1, idx_0)
  oob_idx   <- setdiff(1:nrow(full_data), unique(train_idx))
  if (length(oob_idx) < 10) next

  train_b <- full_data[train_idx, ]
  oob_b   <- full_data[oob_idx, ]

  glm_fit <- tryCatch(glm(y ~ brt + rf, data = train_b, family = binomial()),
                      error = function(e) NULL)
  if (is.null(glm_fit)) next

  coefs <- coef(glm_fit)
  stacking_coefs <- rbind(stacking_coefs, data.frame(
    Boot = b, Intercept = coefs[1], Coef_BRT = coefs[2], Coef_RF = coefs[3]
  ))

  train_std <- train_b %>%
    mutate(brt_std = as.numeric(scale(brt)), rf_std = as.numeric(scale(rf)))
  glm_std <- tryCatch(glm(y ~ brt_std + rf_std, data = train_std, family = binomial()),
                      error = function(e) NULL)
  if (!is.null(glm_std)) {
    sc <- coef(glm_std)
    std_coefs_all <- rbind(std_coefs_all, data.frame(
      Boot = b, Coef_BRT_std = as.numeric(sc["brt_std"]), Coef_RF_std = as.numeric(sc["rf_std"])
    ))
  }

  pred_stack <- predict(glm_fit, oob_b, type = "response")
  pred_qp    <- w_qp[1] * oob_b$brt + w_qp[2] * oob_b$rf

  auc_s <- calc_auc(pred_stack, oob_b$y)
  auc_q <- calc_auc(pred_qp, oob_b$y)
  auc_b <- calc_auc(oob_b$brt, oob_b$y)
  auc_r <- calc_auc(oob_b$rf, oob_b$y)

  tss_s <- max_tss_threshold(pred_stack, oob_b$y)
  tss_b <- max_tss_threshold(oob_b$brt, oob_b$y)
  tss_r <- max_tss_threshold(oob_b$rf, oob_b$y)

  boot_perf <- rbind(boot_perf, data.frame(
    Boot = b,
    AUC_BRT = auc_b, AUC_RF = auc_r,
    AUC_Stacking_Logistic = auc_s, AUC_Stacking_QP = auc_q,
    TSS_BRT = tss_b$TSS, TSS_RF = tss_r$TSS, TSS_Stacking = tss_s$TSS
  ))
}

write.csv(stacking_coefs, file.path(TABDIR, "stacking_logistic_coefficients.csv"), row.names = FALSE)

mean_coefs <- colMeans(stacking_coefs[, c("Intercept", "Coef_BRT", "Coef_RF")], na.rm = TRUE)
cat(sprintf("  Logistic 平均系数: β0=%.4f, β1(BRT)=%.4f, β2(RF)=%.4f\n",
            mean_coefs[1], mean_coefs[2], mean_coefs[3]))

# ---- 标准化权重：优先用 std_coefs_all，失败则回退到 mean_coefs ----
weight_brt <- NA_real_
weight_rf  <- NA_real_
if (nrow(std_coefs_all) > 0) {
  std_mean <- colMeans(std_coefs_all[, c("Coef_BRT_std", "Coef_RF_std")], na.rm = TRUE)
  abs_b <- abs(std_mean[1]); abs_r <- abs(std_mean[2])
  if (is.finite(abs_b + abs_r) && (abs_b + abs_r) > 0) {
    weight_brt <- abs_b / (abs_b + abs_r)
    weight_rf  <- abs_r / (abs_b + abs_r)
  }
}
if (!is.finite(weight_brt) || !is.finite(weight_rf)) {
  abs_b <- abs(mean_coefs[2]); abs_r <- abs(mean_coefs[3])
  if (is.finite(abs_b + abs_r) && (abs_b + abs_r) > 0) {
    weight_brt <- abs_b / (abs_b + abs_r)
    weight_rf  <- abs_r / (abs_b + abs_r)
  } else {
    weight_brt <- 0.5; weight_rf <- 0.5
  }
}

cat(sprintf("  标准化/回退权重: BRT = %.1f%%, RF = %.1f%%\n", weight_brt * 100, weight_rf * 100))
write.csv(data.frame(weight_BRT = weight_brt, weight_RF = weight_rf),
          file.path(TABDIR, "stacking_weights_standardized.csv"), row.names = FALSE)

# =============================================================================
# 4) 方法比较
# =============================================================================
cat("[STACK] 4. 方法比较 ...\n")

mean_auc_logistic <- mean(boot_perf$AUC_Stacking_Logistic, na.rm = TRUE)
mean_auc_qp       <- mean(boot_perf$AUC_Stacking_QP, na.rm = TRUE)
cat(sprintf("  Logistic AUC: %.4f | QP AUC: %.4f\n", mean_auc_logistic, mean_auc_qp))
cat("  → 最终采用 Logistic 方法（有截距、可校准）\n")

# =============================================================================
# 5) 全样本 Logistic 预测 + MaxTSS
# =============================================================================
cat("[STACK] 5. 全样本 Stacking 预测 ...\n")

logit_full <- mean_coefs[1] + mean_coefs[2] * full_data$brt + mean_coefs[3] * full_data$rf
full_data$stack <- 1 / (1 + exp(-logit_full))

tss_stack <- max_tss_threshold(full_data$stack, full_data$y)
stack_auc <- calc_auc(full_data$stack, full_data$y)

cat(sprintf("  全样本: AUC=%.4f | MaxTSS=%.4f (阈值=%.4f)\n",
            stack_auc, tss_stack$TSS, tss_stack$threshold))

write.csv(data.frame(threshold = tss_stack$threshold, TSS = tss_stack$TSS,
                     Sens = tss_stack$sens, Spec = tss_stack$spec, AUC = stack_auc),
          file.path(TABDIR, "stack_threshold.csv"), row.names = FALSE)

# =============================================================================
# 6) 性能对比表
# =============================================================================
cat("[STACK] 6. 性能对比 ...\n")

boot_perf_out <- boot_perf %>% mutate(AUC_Stacking = AUC_Stacking_Logistic)
write.csv(boot_perf_out, file.path(TABDIR, "performance_comparison.csv"), row.names = FALSE)

overall_performance <- data.frame(
  Model    = c("BRT", "RF", "Stacking_Logistic", "Stacking_QP"),
  Mean_AUC = c(mean(boot_perf$AUC_BRT, na.rm = TRUE),
               mean(boot_perf$AUC_RF, na.rm = TRUE),
               mean_auc_logistic, mean_auc_qp),
  SD_AUC   = c(sd(boot_perf$AUC_BRT, na.rm = TRUE),
               sd(boot_perf$AUC_RF, na.rm = TRUE),
               sd(boot_perf$AUC_Stacking_Logistic, na.rm = TRUE),
               sd(boot_perf$AUC_Stacking_QP, na.rm = TRUE)),
  Mean_TSS = c(mean(boot_perf$TSS_BRT, na.rm = TRUE),
               mean(boot_perf$TSS_RF, na.rm = TRUE),
               mean(boot_perf$TSS_Stacking, na.rm = TRUE), NA),
  SD_TSS   = c(sd(boot_perf$TSS_BRT, na.rm = TRUE),
               sd(boot_perf$TSS_RF, na.rm = TRUE),
               sd(boot_perf$TSS_Stacking, na.rm = TRUE), NA)
)
write.csv(overall_performance, file.path(TABDIR, "overall_performance.csv"), row.names = FALSE)

cat("[STACK] 性能汇总:\n")
print(overall_performance)

# =============================================================================
# 7) Stacking 栅格（Logistic 变换）
# =============================================================================
cat("[STACK] 7. 生成 Stacking 栅格 ...\n")

mean_brt <- rast(file.path(OUTDIR, "BRT_mean.tif"))
mean_rf  <- rast(file.path(OUTDIR, "RF_mean.tif"))
sd_brt   <- rast(file.path(OUTDIR, "BRT_sd.tif"))
sd_rf    <- rast(file.path(OUTDIR, "RF_sd.tif"))

logit_rast <- mean_coefs[1] + mean_coefs[2] * mean_brt + mean_coefs[3] * mean_rf
mean_stack <- 1 / (1 + exp(-logit_rast))

p_1_p <- mean_stack * (1 - mean_stack)
var_logit <- mean_coefs[2]^2 * sd_brt^2 + mean_coefs[3]^2 * sd_rf^2
sd_stack <- p_1_p * sqrt(var_logit)

writeRaster(mean_stack, file.path(OUTDIR, "STACK_mean.tif"), overwrite = TRUE)
writeRaster(sd_stack,   file.path(OUTDIR, "STACK_sd.tif"),   overwrite = TRUE)

# =============================================================================
# 8) MaxTSS 二值化 + 面积
# =============================================================================
cat("[STACK] 8. 二值化 + 面积 ...\n")

bin_stack <- mean_stack >= tss_stack$threshold
writeRaster(bin_stack, file.path(OUTDIR, "STACK_binary_MaxTSS.tif"), overwrite = TRUE)

area_km2 <- calc_suitable_area_km2(bin_stack)
write.csv(data.frame(suitable_area_km2 = area_km2),
          file.path(TABDIR, "suitable_area_km2.csv"), row.names = FALSE)
cat(sprintf("  适宜区面积: %.0f km²\n", area_km2))

# =============================================================================
# 9) 集成变量重要性
# =============================================================================
cat("[STACK] 9. 集成变量重要性 ...\n")

tryCatch({
  brt_imp <- read.csv(file.path(TABDIR, "brt_importance_allruns.csv"))
  rf_imp  <- read.csv(file.path(TABDIR, "rf_importance_allruns.csv"))

  brt_sum <- brt_imp %>%
    group_by(var) %>%
    summarise(Importance_BRT = mean(rel.inf, na.rm = TRUE), .groups = "drop")

  rf_sum <- rf_imp %>%
    group_by(var) %>%
    summarise(Importance_RF = mean(perm_imp, na.rm = TRUE), .groups = "drop") %>%
    mutate(Importance_RF = pmax(Importance_RF, 0),
           Importance_RF = Importance_RF / sum(Importance_RF) * 100)

  importance_all <- full_join(brt_sum, rf_sum, by = "var") %>%
    mutate(Importance_BRT = replace(Importance_BRT, is.na(Importance_BRT), 0),
           Importance_RF  = replace(Importance_RF,  is.na(Importance_RF), 0),
           Importance_Stacking = weight_brt * Importance_BRT + weight_rf * Importance_RF) %>%
    arrange(desc(Importance_Stacking))

  write.csv(importance_all, file.path(TABDIR, "integrated_variable_importance.csv"), row.names = FALSE)
  cat("  ✓ 集成变量重要性已保存\n")
}, error = function(e) cat("  ⚠ 失败:", e$message, "\n"))

# =============================================================================
# 10) 保存关键数据（供画图）
# =============================================================================
save(full_data, mean_coefs, stacking_coefs, boot_perf,
     overall_performance, weight_brt, weight_rf, w_qp,
     tss_stack, area_km2,
     file = file.path(OUTDIR, "stacking_results.RData"))

# =============================================================================
# 11) 方法学汇总行（稳健写出）
# =============================================================================
data_info  <- read_data_summary_safe(OUTDIR)
block_info <- read_block_info_safe(OUTDIR)
vars_keep  <- readRDS(file.path(OUTDIR, "selected_vars.rds"))

cat(sprintf("[STACK] Block info: %s | %.1f km | %s\n",
            block_info$method[1], block_info$block_size_km[1], block_info$reason[1]))

methods_row <- data.frame(
  buffer_km = buffer_km,
  block_km = block_info$block_size_km[1],
  n_presence = data_info$n_presence,
  n_pseudoabsence = data_info$n_pseudoabsence,
  n_total = data_info$n_total,
  n_predictors = length(vars_keep),
  predictors = paste(vars_keep, collapse = ";"),
  AUC_BRT  = overall_performance$Mean_AUC[1],
  AUC_RF   = overall_performance$Mean_AUC[2],
  AUC_Stack_Logistic = mean_auc_logistic,
  AUC_Stack_QP       = mean_auc_qp,
  TSS_Stack = tss_stack$TSS,
  Thr_MaxTSS = tss_stack$threshold,
  suitable_area_km2 = area_km2,
  weight_BRT = weight_brt,
  weight_RF  = weight_rf,
  wQP_BRT = w_qp[1],
  wQP_RF  = w_qp[2],
  beta0 = mean_coefs[1],
  beta_BRT = mean_coefs[2],
  beta_RF = mean_coefs[3]
)

write.csv(methods_row, file.path(TABDIR, "methods_report_row.csv"), row.names = FALSE)

cat("══════════════════════════════════════════════════════════\n")
cat("[STACK] 完成！\n")
cat("  Logistic AUC:", round(mean_auc_logistic, 4),
    "| QP AUC:", round(mean_auc_qp, 4), "\n")
cat("  输出:", OUTDIR, "\n")
cat("══════════════════════════════════════════════════════════\n")
