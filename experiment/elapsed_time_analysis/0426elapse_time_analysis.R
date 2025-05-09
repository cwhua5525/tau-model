lambda <- 400

ama_result <- network_lasso_ama_single(
  X_raw_upper = X_raw_upper,
  rho_matrix_adjusted = rho_matrix_adjusted,
  lambda = lambda,
  tol = 1e-6,
  max_iter = 100000,
  quiet = TRUE
)

admm_result <- network_lasso_admm_single(
  X_raw_upper = X_raw_upper,
  rho_matrix_adjusted = rho_matrix_adjusted,
  lambda = lambda,
  max_iter = 1000,
  tol_primal = 1e-4,
  tol_dual = 1e-4,
  rho_admm = 1.0,
  quiet = TRUE  # 如果想要安靜模式，就設 TRUE
)

# 整理 AMA 資料
df_ama <- data.frame(
  time = cumsum(ama_result$loss$time),  # ✅ 累積時間
  obj = ama_result$loss$all,
  method = "AMA"
)

# 整理 ADMM 資料
df_admm <- data.frame(
  time = cumsum(admm_result$loss$time),  # ✅ 累積時間
  obj = admm_result$loss$all,
  method = "ADMM"
)

# 合併
df_plot <- rbind(df_ama, df_admm)

# 畫圖
library(ggplot2)

P.elapsed<- ggplot(df_plot, aes(x = time, y = obj, color = method)) +
  geom_line(size = 1) +
  scale_y_log10() +  # 目標值通常跨數量級，用 logy 比較好看
  labs(
    title = "Objective vs Time",
    x = "Time (seconds)",
    y = "Objective (log scale)",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(P.elapsed)






P2 <- P.elapsed +
      geom_line(size = 1.4) +
      xlim(0,15) +
      labs(
        title = "Objective vs Time, "~lambda~"= 400",
        x = "Time (seconds)",
        y = "Objective (log scale)",
        color = "Method"
      )
print(P2)

