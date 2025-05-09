library(ggplot2)
library(dplyr)

########################################
##### 讀取資料
########################################

df_raw <- readRDS("D:/MEGA/PhD/2025_Q1_spring/8054/project/results/sweep_summary_df_admm_raw.rds")
df_primal <- readRDS("D:/MEGA/PhD/2025_Q1_spring/8054/project/results/sweep_summary_df_admm_primal.rds")
df_all <- readRDS("D:/MEGA/PhD/2025_Q1_spring/8054/project/results/sweep_summary_df_admm_all.rds")

########################################
##### 加上 Method 標籤
########################################

df_raw <- df_raw %>%
  mutate(Method = "Raw")

df_primal <- df_primal %>%
  mutate(Method = "Primal Warm Start")

df_all <- df_all %>%
  mutate(Method = "Full Warm Start")

########################################
##### 合併資料
########################################

plot_df <- bind_rows(df_raw, df_primal, df_all)

########################################
##### 繪圖: Iterations vs Lambda
########################################

P_compare_warmstart <- ggplot(plot_df, aes(x = lambda, y = n_iter, color = Method)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  labs(
    title = "[ADMM] Iterations vs Lambda: Warm Start Comparison",
    x = expression(lambda),
    y = "Number of Iterations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "lines")
  )

print(P_compare_warmstart)
