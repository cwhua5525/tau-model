library(ggplot2)
library(dplyr)

########################################
##### 讀取資料
########################################

df_admm_all <- readRDS("D:/MEGA/PhD/2025_Q1_spring/8054/project/results/sweep_summary_df_admm_all.rds")

########################################
##### 計算縮放比例 (標準化用)
########################################

scale_factor_dist <- max(df_admm_all$n_cluster, na.rm = TRUE) * 1.2 / max(df_admm_all$mean_distance_to_true, na.rm = TRUE)

########################################
##### 建立雙軸資料框
########################################

df_cluster <- df_admm_all %>%
  mutate(curve_type = "Number of Clusters",
         value = n_cluster)

df_dist <- df_admm_all %>%
  mutate(curve_type = "Mean Distance to True",
         value = mean_distance_to_true * scale_factor_dist)

plot_df <- rbind(df_cluster, df_dist)
plot_df$curve_type <- factor(plot_df$curve_type, levels = c("Number of Clusters", "Mean Distance to True"))

########################################
##### 雙軸圖: Cluster & Dist vs Lambda (線條版)
########################################

P_admm_dist <- ggplot(plot_df, aes(x = lambda, y = value, color = curve_type, linetype = curve_type)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Number of Clusters",
    sec.axis = sec_axis(~ . / scale_factor_dist, name = "Mean Distance to True")
  ) +
  scale_color_manual(
    values = c("Number of Clusters" = "#E41A1C", "Mean Distance to True" = "#377EB8")
  ) +
  labs(
    title = "[ADMM] Clusters and Estimation Error vs Lambda",
    x = expression(lambda),
    color = "Curve",
    linetype = "Curve"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(P_admm_dist)
