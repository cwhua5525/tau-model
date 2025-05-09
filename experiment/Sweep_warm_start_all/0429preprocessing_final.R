# ==============================================
# 載入函數
# ==============================================
source("utils.R") #生成相關
source("algorithms_v4.R") #AMA/ADMM使用函數
# ==============================================
# 主流程：生成多亞型資料與合併訓練集
# ==============================================

set.seed(8054)  # 方便重現結果

num_subtype <- 3
all_subtypes <- vector("list", num_subtype)

for (subtype_idx in 1:num_subtype) {
  cat("Generating subtype", subtype_idx, "...\n")
  
  p.dim <- 68
  
  tmp_res <- generate_subtype_data(
    d = 3,
    p = p.dim,
    target_density = 0.2,
    mu = 0.2*p.dim,
    sigma2 = 0.5,
    n_data = 20
  )
  all_subtypes[[subtype_idx]] <- tmp_res
}

# -------------------------------------------------------------
# (E) 將所有資料合成訓練集
# -------------------------------------------------------------
# 將所有 subtype 的資料合併成一個大的 list
combined_data <- list()

obs_idx <- 1  # 觀測值索引

for (subtype_id in seq_along(all_subtypes)) {
  subtype_data <- all_subtypes[[subtype_id]]$data
  for (i in seq_along(subtype_data)) {
    # 單筆資料取出
    entry <- subtype_data[[i]]
    
    # 建立新的觀測值 list
    combined_data[[obs_idx]] <- list(
      hGamma  = entry$hGamma,
      y       = entry$y,
      Tval       = entry$Tval,
      subtype = subtype_id
    )
    
    obs_idx <- obs_idx + 1
  }
}

# combined_data 現在是長度 1000 的 list，每筆是一個 list，包含 4 個欄位
data <- combined_data[[1]]
str(data)

data$y %>% sum()

# -------------------------------------------------------------
# (D) 結果結構說明
# -------------------------------------------------------------
# all_subtypes 是一個長度=10的 list，每一個元素對應一個亞型，包含:
#   - $Omega :  該亞型的 baseline 連接性矩陣  (p x p)
#   - $G     :  該亞型的 adjacency matrix (p x p)
#   - $s0    :  該亞型的傳播初始節點
#   - $data  :  list(100) ，其中每個元素是一筆 (hGamma, y, Tval)

# -----------------------------------------------------------------
# 假設您已經有 all_subtypes 結構:
#   - all_subtypes[[j]]$Omega  (p x p 矩陣)
#   - all_subtypes[[j]]$data[[i]]$hGamma (p x p 矩陣)
#
# 調用函數:
res <- compute_within_between(all_subtypes)
cat("組內距離 (平均) =", res$within_mean, "\n",
    "組間距離 (平均) =", res$between_mean, "\n",
    "比值 (within/between) =", res$ratio, "\n")


# ==============================================
# ✅ Network Lasso (ADMM) - Full Graph with Convergence and Clustering Count
# ==============================================


# ==============================================
# (6) 資料前處理：生成 X_raw , rho_matrix , rho_matrix_adjusted
# ==============================================
N <- length(combined_data)
p <- nrow(combined_data[[1]]$hGamma)
d <- p * p
d_adjusted <- p * (p - 1) / 2 # 調整過的維度

X_raw <- matrix(NA, nrow = N, ncol = d)
Y_mat <- matrix(NA, nrow = N, ncol = p)
true_labels <- integer(N)

for (i in 1:N) {
  hGamma <- combined_data[[i]]$hGamma
  y <- combined_data[[i]]$y
  y_diff <- abs(outer(y, y, "-"))
  y_sum  <- abs(outer(y, y, "+"))
  adjusted <- hGamma - y_diff - y_sum
  X_raw[i, ] <- as.vector(adjusted)
  
  Y_mat[i, ] <- combined_data[[i]]$y
  true_labels[i] <- combined_data[[i]]$subtype
}


##### 相似度矩陣計算
rho_matrix <- matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  for (j in 1:N) {
    yi <- Y_mat[i, ]
    yj <- Y_mat[j, ]
    denom <- max(sum(yi), sum(yj))
    if (denom > 0) {
      rho_matrix[i, j] <- sum(yi * yj) / denom
    }
  }
}

##### 相似度矩陣視覺化
image(rho_matrix)

# 建構 X_raw：將每筆資料的 adjusted Gamma 上三角展平
X_raw_upper <- t(sapply(1:N, function(i) {
  hGamma <- combined_data[[i]]$hGamma
  y      <- combined_data[[i]]$y
  y_diff <- abs(outer(y, y, "-"))
  y_sum  <- abs(outer(y, y, "+"))
  adjusted <- hGamma - y_diff - y_sum
  get_upper_triangle(adjusted)
}))

# 構圖（block 模式下做健全性測試）
#rho_matrix_adjusted <- build_graph_matrix(rho_matrix, type = "block", num_block = 3)
#rho_matrix_adjusted <- build_graph_matrix(rho_matrix, type = "full")
rho_matrix_adjusted <- build_graph_matrix(rho_matrix, type = "thres",thres = quantile(rho_matrix,0.5))  
#rho_matrix_adjusted <- build_graph_matrix(rho_matrix, type = "topk", k = 5) 

image(rho_matrix_adjusted)
