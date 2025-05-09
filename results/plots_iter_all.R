library(ggplot2)
library(dplyr)

########################################
##### 讀取資料
########################################

df_admm_all <- readRDS("D:/MEGA/PhD/2025_Q1_spring/8054/project/results/sweep_summary_df_admm_all.rds")
df_ama_all <- readRDS("D:/MEGA/PhD/2025_Q1_spring/8054/project/results/sweep_summary_df_ama_all.rds")

########################################
##### 計算縮放比例 (iter標準化：用log10)
########################################

max_iter_log10 <- max(log10(c(df_admm_all$n_iter, df_ama_all$n_iter)), na.rm = TRUE)
scale_factor_iter <- max(df_admm_all$n_cluster, na.rm = TRUE) * 1.2 / max_iter_log10

########################################
##### 建立資料框
########################################

# Cluster 數（ADMM）
df_cluster <- df_admm_all %>%
  mutate(curve_type = "Number of Clusters",
         value = n_cluster)

# ADMM iter (log10)
df_admm_iter <- df_admm_all %>%
  mutate(curve_type = "ADMM Iterations",
         value = log10(n_iter) * scale_factor_iter)

# AMA iter (log10)
df_ama_iter <- df_ama_all %>%
  mutate(curve_type = "AMA Iterations",
         value = log10(n_iter) * scale_factor_iter)

# 合併
plot_df <- bind_rows(df_cluster, df_admm_iter, df_ama_iter)
plot_df$curve_type <- factor(plot_df$curve_type, levels = c("Number of Clusters", "ADMM Iterations", "AMA Iterations"))

########################################
##### 雙軸圖: Cluster & log10(Iterations) vs Lambda
########################################

# 設定右側 ticks (log scale)
iter_ticks <- c(1, 10, 100, 1000, 10000)
iter_ticks_log <- log10(iter_ticks)

P_compare_log <- ggplot(plot_df, aes(x = lambda, y = value, color = curve_type)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Number of Clusters",
    sec.axis = sec_axis(
      ~ . / scale_factor_iter,
      name = "Iteration (log scale)",
      breaks = iter_ticks_log,
      labels = iter_ticks
    )
  ) +
  scale_color_manual(
    values = c(
      "Number of Clusters" = "#666666",
      "ADMM Iterations" = "#E41A1C",
      "AMA Iterations" = "#00CED1"
    )
  ) +
  guides(
    color = guide_legend(title = "Curve"),
    linetype = "none",
    linewidth = "none",
    size = "none"
  ) +
  labs(
    title = "[ADMM vs AMA] Cluster Count and Iterations vs"~ lambda,
    x = expression(lambda)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "lines")
  )

print(P_compare_log)
