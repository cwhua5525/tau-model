# ============================================================================
# 📄 AMA Sweep Script
# Description: Sweep over a sequence of lambdas using AMA solver for Network Lasso
# ============================================================================

# 先設定 lambda 範圍
lambda_min <- 1
lambda_max <- 1000
lambda_list <- exp(seq(log(lambda_min), log(lambda_max), length.out = 100))

# 用來儲存結果的 list
sweep_results_ama <- list()

# Summary dataframe 的初始化
summary_df_ama <- data.frame(
  lambda = numeric(0),
  n_iter = integer(0),
  exec_time = numeric(0),
  n_cluster = integer(0),
  mean_intercluster_dist = numeric(0),
  mean_distance_to_true = numeric(0)
)

# 提取 baseline true Omega
Omega_true_list <- lapply(all_subtypes, function(sub) sub$Omega)

# 開始 sweep
for (lambda_val in lambda_list) {
  cat("======= Running AMA lambda =", lambda_val, "=======\n")
  
  start_time <- Sys.time()
  
  ama_result <- network_lasso_ama_single(
    X_raw_upper = X_raw_upper,
    rho_matrix_adjusted = rho_matrix_adjusted,
    lambda = lambda_val,
    max_iter = 10000,
    tol = 1e-5,
    quiet = TRUE
  )
  
  end_time <- Sys.time()
  exec_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  est_array <- ama_result$U
  N <- nrow(est_array)
  pairwise_dist <- matrix(0, N, N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      pairwise_dist[i, j] <- sqrt(sum((est_array[i, ] - est_array[j, ])^2))
      pairwise_dist[j, i] <- pairwise_dist[i, j]
    }
  }
  
  cluster_labels <- rep(NA, N)
  current_label <- 1
  for (i in 1:N) {
    if (is.na(cluster_labels[i])) {
      cluster_labels[i] <- current_label
      cluster_labels[pairwise_dist[i, ] < 1e-4] <- current_label
      current_label <- current_label + 1
    }
  }
  n_cluster <- length(unique(cluster_labels))
  
  cluster_centroids <- t(sapply(unique(cluster_labels), function(cl) {
    colMeans(est_array[cluster_labels == cl, , drop=FALSE])
  }))
  
  n_centroids <- nrow(cluster_centroids)
  
  if (n_centroids > 1) {
    inter_dist_list <- c()
    for (i in 1:(n_centroids-1)) {
      for (j in (i+1):n_centroids) {
        inter_dist_list <- c(inter_dist_list, sqrt(sum((cluster_centroids[i, ] - cluster_centroids[j, ])^2)))
      }
    }
    mean_intercluster_dist <- mean(inter_dist_list)
  } else {
    mean_intercluster_dist <- NA
  }
  
  restore_full_matrix <- function(upper_vec, p) {
    mat <- matrix(0, p, p)
    mat[upper.tri(mat)] <- upper_vec
    mat <- mat + t(mat)
    diag(mat) <- 0
    return(mat)
  }
  
  p <- (1 + sqrt(1 + 8 * ncol(est_array))) / 2
  Omega_est_list <- lapply(1:N, function(i) restore_full_matrix(est_array[i, ], p))
  
  dist_to_true <- matrix(0, nrow = length(Omega_true_list), ncol = N)
  for (i in 1:length(Omega_true_list)) {
    for (j in 1:N) {
      dist_to_true[i, j] <- sqrt(sum((get_upper_triangle(Omega_true_list[[i]]) - get_upper_triangle(Omega_est_list[[j]]))^2))
    }
  }
  mean_distance_to_true <- mean(apply(dist_to_true, 2, min))
  
  summary_df_ama <- rbind(summary_df_ama, data.frame(
    lambda = lambda_val,
    n_iter = ama_result$n_iter,
    exec_time = exec_time,
    n_cluster = n_cluster,
    mean_intercluster_dist = mean_intercluster_dist,
    mean_distance_to_true = mean_distance_to_true
  ))
  
  sweep_results_ama[[as.character(lambda_val)]] <- list(
    ama_result = ama_result,
    cluster_labels = cluster_labels
  )
  
  # (3) ✅ 新增：儲存目前 sweep 進度
  saveRDS(summary_df_ama, file = "sweep_summary_df_ama.rds")
  saveRDS(sweep_results_ama, file = "sweep_results_ama.rds")
}

print(summary_df_ama)

# ============================================================================
# 📈 做圖
# ============================================================================
library(ggplot2)
library(dplyr)

# -------------------------------
# 📈 群數 vs Mean Distance to True
# -------------------------------

# 計算縮放比例
scale_factor_dist <- max(summary_df_ama$n_cluster, na.rm = TRUE) * 1.2 / max(summary_df_ama$mean_distance_to_true, na.rm = TRUE)

# 建資料
df_cluster <- summary_df_ama %>%
  mutate(curve_type = "Number of Clusters",
         value = n_cluster)

df_dist <- summary_df_ama %>%
  mutate(curve_type = "Mean Distance to True",
         value = mean_distance_to_true * scale_factor_dist)

# 合併
plot_df_dist <- rbind(df_cluster, df_dist)
plot_df_dist$curve_type <- factor(plot_df_dist$curve_type, levels = c("Number of Clusters", "Mean Distance to True"))

# 畫圖
P_ama_dist <- ggplot(plot_df_dist, aes(x = lambda, y = value, color = curve_type, linetype = curve_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Number of Clusters",
    sec.axis = sec_axis(~ . / scale_factor_dist, name = "Mean Distance to True")
  ) +
  scale_color_manual(
    values = c("Number of Clusters" = "#E41A1C", "Mean Distance to True" = "#377EB8")
  ) +
  labs(
    title = "[AMA] Clusters and Estimation Error vs Lambda",
    x = expression(lambda),
    color = "Curve",
    linetype = "Curve"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(P_ama_dist)

# -------------------------------
# 📈 群數 vs Iterations
# -------------------------------

# 計算縮放比例
scale_factor_iter <- max(summary_df_ama$n_cluster, na.rm = TRUE) * 1.2 / max(summary_df_ama$n_iter, na.rm = TRUE)

# 建資料
df_cluster <- summary_df_ama %>%
  mutate(curve_type = "Number of Clusters",
         value = n_cluster)

df_iter <- summary_df_ama %>%
  mutate(curve_type = "Number of Iterations",
         value = n_iter * scale_factor_iter)

# 合併
plot_df_iter <- rbind(df_cluster, df_iter)
plot_df_iter$curve_type <- factor(plot_df_iter$curve_type, levels = c("Number of Clusters", "Number of Iterations"))

# 畫圖
P_ama_iter <- ggplot(plot_df_iter, aes(x = lambda, y = value, color = curve_type, linetype = curve_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Number of Clusters",
    sec.axis = sec_axis(~ . / scale_factor_iter, name = "Number of Iterations")
  ) +
  scale_color_manual(
    values = c("Number of Clusters" = "#E41A1C", "Number of Iterations" = "#377EB8")
  ) +
  labs(
    title = "[AMA] Clusters and Iterations vs Lambda",
    x = expression(lambda),
    color = "Curve",
    linetype = "Curve"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(P_ama_iter)

# -------------------------------
# 📈 群數 vs Execution Time
# -------------------------------

# 計算縮放比例
scale_factor_time <- max(summary_df_ama$n_cluster, na.rm = TRUE) * 1.2 / max(summary_df_ama$exec_time, na.rm = TRUE)

# 建資料
df_cluster <- summary_df_ama %>%
  mutate(curve_type = "Number of Clusters",
         value = n_cluster)

df_time <- summary_df_ama %>%
  mutate(curve_type = "Execution Time",
         value = exec_time * scale_factor_time)

# 合併
plot_df_time <- rbind(df_cluster, df_time)
plot_df_time$curve_type <- factor(plot_df_time$curve_type, levels = c("Number of Clusters", "Execution Time"))

# 畫圖
P_ama_time <- ggplot(plot_df_time, aes(x = lambda, y = value, color = curve_type, linetype = curve_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Number of Clusters",
    sec.axis = sec_axis(~ . / scale_factor_time, name = "Execution Time (seconds)")
  ) +
  scale_color_manual(
    values = c("Number of Clusters" = "#E41A1C", "Execution Time" = "#377EB8")
  ) +
  labs(
    title = "[AMA] Clusters and Execution Time vs Lambda",
    x = expression(lambda),
    color = "Curve",
    linetype = "Curve"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(P_ama_time)