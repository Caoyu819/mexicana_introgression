#!/usr/bin/Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
f_corr <- args[1]
f_fig <- args[2]

# 读入数据
df <- read.table(f_corr, header = TRUE, sep = "\t")

# 绘制Spearman_rho的直方图
p_dist <- ggplot(df, aes(x = Spearman_rho)) +
  geom_histogram(binwidth = 0.02, fill = "#006400", color = "black", alpha = 0.8) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  labs(x = "Spearman's rho", y = "Number of samples")

ggsave(f_fig, p_dist, width = 6, height = 4)
