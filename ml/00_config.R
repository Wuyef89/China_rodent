## =============================================================================
## 00_config.R  —  全局配置（所有脚本共用）
## =============================================================================

# ---- 线程安全（HPC 必需）----
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

# ---- 1. 路径 ----------------------------------------------------------------
BASE_DIR    <- "/data/home/test01/rodent"
DIR_DATA    <- file.path(BASE_DIR, "data")
DIR_SCRIPT  <- file.path(BASE_DIR, "script")
DIR_OUT     <- file.path(BASE_DIR, "out")
DIR_LOG     <- file.path(BASE_DIR, "logs")

# 栅格输入
DIR_ENV_TRAIN <- file.path(DIR_DATA, "rasters_clean", "train")
DIR_ENV_TEST  <- file.path(DIR_DATA, "rasters_clean", "test")
PRES_RASTER   <- file.path(DIR_DATA, "rasters_clean", "response_variable", "RodentsPoints.tif")

# 中国边界（栅格已预处理裁剪对齐，建模阶段不再使用）
# CHINA_SHP      <- file.path(DIR_DATA, "mask", "China", "China.shp")
# USE_CHINA_MASK <- file.exists(CHINA_SHP)
USE_CHINA_MASK <- FALSE

# 输出目录
buffer_out_dir <- function(bkm) file.path(DIR_OUT, paste0("buffer_", bkm, "km"))

dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)

# ---- 2. 缓冲区（固定 200 km）----
BUFFER_KM <- 200

# ---- 3. 伪缺失 / SRE 参数 ----
PB_RATIO       <- 1          # presence:background = 1:1
SRE_OVERSAMPLE <- 8          # 超采样倍数
SRE_VARS       <- c("TEM", "PRE", "RHU", "PRS", "WIN", "EVP", "SSD", "GST")
SRE_N_VARS     <- 8          # fallback

# ---- 4. 共线性筛选：已取消 ----
# BRT/RF 对共线性天然容忍，保留全部变量以获得更好的生态可解释性
# 如需恢复，取消下面的注释即可
# CORR_CUTOFF   <- 0.8
# VIF_CUTOFF    <- 10
# COLL_SAMPLE_N <- 50000

# ---- 5. 空间自相关 / block size ----
AUTORANGE_MAX_LAYERS <- 12
BLOCK_CLAMP_MIN_M    <- 50000    # 最小 50 km

# 获取环境变量（假设传入的是 200, 400, 600）
env_val <- as.numeric(Sys.getenv("FORCE_BLOCK_SIZE_KM"))

if (!is.na(env_val)) {  # 确保这里用的是 env_val
  BLOCK_CLAMP_MAX_M <- env_val * 1000
  BLOCK_FORCE_FIXED <- TRUE  
  cat("[CONFIG] 强制锁定 Block Size 为:", env_val, "km\n")
} else {
  BLOCK_CLAMP_MAX_M <- 400000  
  BLOCK_FORCE_FIXED <- FALSE
}

# 输出目录函数也要跟着变，否则文件夹名字会错乱
buffer_out_dir <- function(bkm) {
  path <- file.path(DIR_OUT, paste0("buffer_", bkm, "km"))
  if (!is.na(env_val)) {
    # 确保文件夹名显示的是 km
    path <- paste0(path, "_block", env_val, "km")
  }
  return(path)
}

# ---- 6. 变量筛选阶段（Stage A: 3×5 空间 CV）----
VARSEL_REPEATS  <- 3
VARSEL_FOLDS    <- 5
IMP_THRESH_MEAN <- 2.0       # 平均 rel.inf > 2%
KEEP_PROP       <- 0.60      # >= 60% 折叠中贡献 > 0
# 不设 fallback，筛出几个就用几个

# ---- 7. 最终建模阶段（Stage B: 10×10 空间 CV）----
FINAL_REPEATS <- 10
FINAL_FOLDS   <- 10

# ---- 8. BRT 超参数 ----
BRT_TC   <- 5
BRT_LR   <- 0.001
BRT_BF   <- 0.5
BRT_MAXT <- 20000

# ---- 9. RF 超参数 ----
RF_TREES <- 1500

# ---- 10. 并行设置 ----
# N_WORKERS 控制 future_lapply 的 worker 数
# 注意：每个 worker 会独立加载完整栅格做空间预测，内存占用 = N_WORKERS × 栅格大小
# 因此不宜全开，建议 3-4 个（HPC 64GB 内存下安全）
MAX_WORKERS <- 4   # ← 安全上限，根据 HPC 内存调整
N_WORKERS <- min(
  MAX_WORKERS,
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
)
if (!is.finite(N_WORKERS) || N_WORKERS < 1) N_WORKERS <- 1

cat("[CONFIG] 配置加载完毕。BUFFER_KM =", BUFFER_KM, "\n")