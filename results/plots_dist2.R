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
##### 建立雙軸資料框 (三條線)
########################################

df_cluster <- df_admm_all %>%
  mutate(curve_type = "Number of Clusters",
         value = n_cluster)

df_dist <- df_admm_all %>%
  mutate(curve_type = "Mean Distance to True",
         value = mean_distance_to_true * scale_factor_dist)

df_inter <- df_admm_all %>%
  mutate(curve_type = "Mean Inter-cluster Distance",
         value = mean_intercluster_dist * scale_factor_dist)

plot_df <- bind_rows(df_cluster, df_dist, df_inter)
plot_df$curve_type <- factor(plot_df$curve_type, levels = c("Number of Clusters", "Mean Distance to True", "Mean Inter-cluster Distance"))

########################################
##### 找第24行作為標記點
########################################

point_24 <- df_inter[24, ]

########################################
##### 雙軸圖: Cluster, Dist, Inter-dist vs Lambda (清理圖例版)
########################################

P_admm_dist2 <- ggplot(plot_df, aes(x = lambda, y = value, color = curve_type)) +
  geom_line(linewidth = 1.2) +  # 不再用 aes() 控制線型或粗細
  # 加標記點
  geom_point(data = point_24, aes(x = lambda, y = mean_intercluster_dist * scale_factor_dist),
             color = "#238B45", size = 3) +
  scale_x_log10() +
  scale_y_continuous(
    name = "Number of Clusters",
    sec.axis = sec_axis(~ . / scale_factor_dist, name = "Distance")
  ) +
  scale_color_manual(
    values = c(
      "Number of Clusters" = "#E41A1C",
      "Mean Distance to True" = "#377EB8",
      "Mean Inter-cluster Distance" = "#238B45"
    )
  ) +
  guides(
    color = guide_legend(title = ""),
    linetype = "none",
    linewidth = "none",
    size = "none"
  ) +
  labs(
    title = "Clustering Behaviors",
    x = expression(lambda)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "lines")
  )

print(P_admm_dist2)
