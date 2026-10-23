#!/usr/bin/env Rscript

# step2_plot_rec_admix.R
# Usage: Rscript step2_plot_rec_admix.R <f_merged> <out_dir> <out_prefix>

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript step2_plot_rec_admix.R <f_merged> <out_dir> <out_prefix>")
}

f_merged   <- args[1]
out_dir    <- args[2]
out_prefix <- args[3]

if (!file.exists(f_merged)) stop("Merged file not found: ", f_merged)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Read merged dataset
combined_df <- fread(f_merged)

# 2. Filter non-zero, finite rates and extreme quantiles (0.5% - 99.5%)
analysis_df <- combined_df[!is.na(rec_cmMb) & is.finite(rec_cmMb) &
                           !is.na(mexicana_dosage) & is.finite(mexicana_dosage) &
                           rec_cmMb != 0]

if (nrow(analysis_df) < 10) stop("Insufficient data points for analysis.")

q_bounds <- quantile(analysis_df$rec_cmMb, probs = c(0.005, 0.995), na.rm = TRUE)
data_filtered <- analysis_df[rec_cmMb >= q_bounds[1] & rec_cmMb <= q_bounds[2]]

# 3. Correlation tests
pear_test  <- cor.test(data_filtered$rec_cmMb, data_filtered$mexicana_dosage, method = "pearson")
spear_test <- cor.test(data_filtered$rec_cmMb, data_filtered$mexicana_dosage, method = "spearman")

stat_out <- file.path(out_dir, paste0(out_prefix, ".spearman_result_allsamples.txt"))
fwrite(data.table(
  pearson_r    = pear_test$estimate,
  pearson_p    = pear_test$p.value,
  spearman_rho = spear_test$estimate,
  spearman_p   = spear_test$p.value,
  N_windows    = nrow(data_filtered)
), stat_out, sep = "\t")

# 4. Correlation hexbin scatter plot
p_corr <- ggplot(data_filtered, aes(x = rec_cmMb, y = mexicana_dosage)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(trans = "log10") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +
  labs(
    x = "Recombination rate (cM/Mb)",
    y = "Mexicana dosage",
    title = "Recombination Rate vs. Mexicana Dosage",
    subtitle = sprintf("Pearson r = %.3f (p = %.2e) | Spearman rho = %.3f (p = %.2e)",
                       pear_test$estimate, pear_test$p.value, spear_test$estimate, spear_test$p.value)
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(out_dir, paste0(out_prefix, ".CorrRecAdmix_99percentiles.pdf")), p_corr, width = 8, height = 7)
ggsave(file.path(out_dir, paste0(out_prefix, ".CorrRecAdmix_99percentiles.png")), p_corr, width = 8, height = 7, dpi = 300)

# 5. Decile boxplot
valid_dec <- copy(analysis_df)
decile_breaks <- unique(quantile(valid_dec$rec_cmMb, probs = seq(0, 1, by = 0.1), na.rm = TRUE))
decile_labels <- sapply(seq_len(length(decile_breaks) - 1), function(i) {
  paste0(round(decile_breaks[i], 3), " - ", round(decile_breaks[i + 1], 3))
})
valid_dec[, rec_decile := cut(rec_cmMb, breaks = decile_breaks, include.lowest = TRUE, labels = decile_labels)]

y_max <- valid_dec %>%
  group_by(rec_decile) %>%
  summarise(upper = quantile(mexicana_dosage, 0.75, na.rm = TRUE) + 1.5 * IQR(mexicana_dosage, na.rm = TRUE), .groups = "drop") %>%
  pull(upper) %>% max(na.rm = TRUE)
if (!is.finite(y_max) || is.na(y_max)) y_max <- max(valid_dec$mexicana_dosage, na.rm = TRUE)

p_decile <- ggplot(valid_dec[!is.na(rec_decile)], aes(x = rec_decile, y = mexicana_dosage)) +
  geom_boxplot(outlier.shape = NA, fill = "coral") +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(x = "Recombination Rate Decile (cM/Mb)", y = "Mexicana Dosage") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

ggsave(file.path(out_dir, paste0(out_prefix, ".DecileBoxplot_RecAdmix.pdf")), p_decile, width = 9, height = 6)
ggsave(file.path(out_dir, paste0(out_prefix, ".DecileBoxplot_RecAdmix.png")), p_decile, width = 9, height = 6, dpi = 300)

message("Completed plotting.")
