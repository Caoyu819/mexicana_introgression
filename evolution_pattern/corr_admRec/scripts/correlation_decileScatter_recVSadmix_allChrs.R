#!/usr/bin/env Rscript

# correlation_rec_admix.R
# Usage: Rscript correlation_rec_admix.R <dir_elai> <dir_out> <sample_id> [mg] [sw]
#
# Example:
# Rscript correlation_rec_admix.R /public/.../out_elai /public/.../out SAMPLE001 ave10runs ps21

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpointdensity)
  library(viridis)
})

# ---------------------------
# Parse arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: correlation_rec_admix.R <dir_elai> <rec_map_file> <chrlen_file> <dir_out> <sample_id> [mg] [sw]")
}

dir_elai  <- args[1]
rec_map_file <- args[2]
chrlen_file <- args[3]
dir_out   <- args[4]
sample_id <- args[5]
mg        <- args[6]
sw        <- args[7]

if (FALSE){
dir_elai <- "${MEXICANA_INTROGRESSION_ROOT}/0_elai/out_elai_845maize"
rec_map_file <- "${MEXICANA_INTROGRESSION_ROOT}/1_ahmm_mls/input_ahmms/rec_map_B73V4.txt"
chrlen_file  <- "${MEXICANA_INTROGRESSION_ROOT}/0_list_map/B73_v4.genome.index"
dir_out <- "${MEXICANA_INTROGRESSION_ROOT}/0_elai/out_elai_845maize/visualization/s10w10_mg1000_unmask"
sample_id <- "B73.amp"
mg <- "mg1000"
sw <- "s10w10"
}

message("Running with:")
message(sprintf(" dir_elai = %s\n rec_map_file = %s\n chrlen_file = %s\n dir_out = %s\n sample_id = %s\n mg = %s\n sw = %s",
                dir_elai, rec_map_file, chrlen_file, dir_out, sample_id, mg, sw))

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Helper: safe column selector
# ---------------------------
get_col_name <- function(dt, possible_names) {
  for (nm in possible_names) {
    if (nm %in% names(dt)) return(nm)
  }
  return(NULL)
}

# ---------------------------
# Read maps used for interpolation
# ---------------------------
#rec_map_file <- "${MEXICANA_INTROGRESSION_ROOT}/1_ahmm_mls/input_ahmms/rec_map_B73V4.txt"
#chrlen_file  <- "${MEXICANA_INTROGRESSION_ROOT}/0_list_map/B73_v4.genome.index"

if (!file.exists(rec_map_file)) stop("Recombination map file not found: ", rec_map_file)
if (!file.exists(chrlen_file)) stop("Chromosome length file not found: ", chrlen_file)

newmap <- fread(rec_map_file, header = FALSE)
# Expecting columns: S_id, M_id, cM, chr, POS  (as in original)
# If your file has different format, adjust the setnames accordingly
setnames(newmap, c("S_id", "M_id", "cM", "chr", "POS"))

chrlen <- fread(chrlen_file, header = FALSE)
setnames(chrlen, c("chr", "length"))

# Ensure types
newmap[, POS := as.numeric(POS)]
newmap[, cM := as.numeric(cM)]

# Build interpolation functions per chromosome
interp_cm <- list()
for (i in seq_len(nrow(chrlen))) {
  ch <- chrlen$chr[i]
  map_chr <- newmap[chr == ch]
  if (nrow(map_chr) >= 2) {
    # ensure ordering by POS
    map_chr <- map_chr[order(POS)]
    interp_cm[[as.character(ch)]] <- approxfun(x = map_chr$POS, y = map_chr$cM, rule = 2)
  } else {
    warning("Chromosome ", ch, " has <2 markers; skipping interpolation for this chromosome.")
  }
}

# ---------------------------
# Read admixture per chromosome and compute rec rate
# ---------------------------
chromosomes <- 1:10
dat_list <- vector("list", length(chromosomes))

for (chr_id in chromosomes) {
  fpath <- file.path(dir_elai, paste0("chr", chr_id),
                     mg, sw,
                     paste0(sample_id, ".chr", chr_id, ".", mg, ".", sw, ".ave10runs.ps21.mex.addRec.txt"))
  if (!file.exists(fpath)) {
    warning("File not found for Chr", chr_id, ": ", fpath, " -- skipping this chromosome.")
    next
  }
  message("Reading: ", fpath)
  admixmap <- fread(fpath)

  # Try to detect important columns: chromosome, POS, mexicana dosage
  chr_col <- get_col_name(admixmap, c("CHROM", "Chr", "chr", "CH"))
  pos_col <- get_col_name(admixmap, c("POS", "pos", "POSITION", "position"))
  mex_col <- get_col_name(admixmap, c("mexicana_dosage", "mexicana", "Mean", "mean", "dosage"))

  if (is.null(pos_col) || is.null(mex_col)) {
    stop("Input file ", fpath, " does not contain expected POS or mexicana dosage columns. Found columns: ", paste(names(admixmap), collapse = ", "))
  }

  # Standardize column names
  setnames(admixmap, old = pos_col, new = "POS", skip_absent = TRUE)
  setnames(admixmap, old = mex_col, new = "mexicana_dosage", skip_absent = TRUE)
  if (!("CHROM" %in% names(admixmap))) {
    # create CHROM column from chr_id if absent
    admixmap[, CHROM := as.character(chr_id)]
  } else {
    setnames(admixmap, old = chr_col, new = "CHROM", skip_absent = TRUE)
  }

  # Ensure numeric types
  admixmap[, POS := as.numeric(POS)]
  admixmap[, mexicana_dosage := as.numeric(mexicana_dosage)]

  # Compute total_cm using interpolation function if available
  if (!is.null(interp_cm[[as.character(chr_id)]])) {
    admixmap[, total_cm := interp_cm[[as.character(chr_id)]](POS)]
  } else {
    admixmap[, total_cm := NA_real_]
  }

  # Compute diffs by row order (assumes data sorted by POS)
  admixmap <- admixmap[order(POS)]
  admixmap[, diff_cm := c(NA, diff(total_cm))]
  admixmap[, diff_bp := c(NA, diff(POS))]
  admixmap[, rec_cmMb := as.numeric( diff_cm / (diff_bp / 1e6) )]  # may produce NA/Inf
  # Replace infinite/NaN with NA
  admixmap[!is.finite(rec_cmMb), rec_cmMb := NA_real_]

  # Add Chr factor column used for faceting/combining
  admixmap[, Chr := paste0("Chr", chr_id)]

  # Save per-chromosome file (optional, into dir_out)
  out_chr_file <- file.path(dir_elai, paste0("chr", chr_id),
                     mg, sw,
                     paste0(sample_id, ".chr", chr_id, ".", mg, ".", sw, ".ave10runs.ps21.mex.addRec.txt"))

  fwrite(admixmap, file = out_chr_file, sep = "\t", quote = FALSE, col.names = TRUE)
  message("Wrote: ", out_chr_file)

  dat_list[[as.character(chr_id)]] <- admixmap
}

# Combine
combined_admixmap <- rbindlist(dat_list, use.names = TRUE, fill = TRUE)
if (nrow(combined_admixmap) == 0) stop("No chromosome files were read successfully. Exiting.")

combined_admixmap$Chr <- factor(combined_admixmap$Chr, levels = paste0("Chr", chromosomes))

# ---------------------------
# Correlation analysis
# ---------------------------
# Keep only rows with finite rec_cmMb and mexicana_dosage
analysis_df <- combined_admixmap[!is.na(rec_cmMb) & is.finite(rec_cmMb) & !is.na(mexicana_dosage) & is.finite(mexicana_dosage)]

if (nrow(analysis_df) < 10) stop("Too few datapoints for correlation analysis after filtering rec_cmMb and mexicana_dosage.")

# Remove exact zeros in rec_cmMb if desired (original code excluded diff_cmMb != 0)
analysis_df <- analysis_df[rec_cmMb != 0]

# Filter to central quantiles (0.5% - 99.5%) to remove extreme tails
lower_bound <- quantile(analysis_df$rec_cmMb, probs = 0.005, na.rm = TRUE)
upper_bound <- quantile(analysis_df$rec_cmMb, probs = 0.995, na.rm = TRUE)
data_filtered <- analysis_df[rec_cmMb >= lower_bound & rec_cmMb <= upper_bound]

# Spearman correlation
spearman_test <- cor.test(data_filtered$rec_cmMb, data_filtered$mexicana_dosage,
                          method = "spearman", exact = FALSE)
spearman_rho <- unname(spearman_test$estimate)
p_val <- spearman_test$p.value
p_label <- ifelse(p_val < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(p_val, 3)))
label_text <- sprintf("Spearman's rho = %.2f\n%s", spearman_rho, p_label)

out_df <- data.frame(
  Sample = sample_id,
  Spearman_rho = spearman_rho,
  P_value = p_label
)

out_file <- file.path(dir_out, paste0("spearman_result_allsamples.txt"))
write.table(out_df, file = out_file, sep = "\t", row.names = FALSE, col.names = !file.exists(out_file), append = TRUE, quote = FALSE)

if (FALSE) {
# Plot: density-colored scatter
p_corr <- ggplot(data_filtered, aes(x = rec_cmMb, y = mexicana_dosage)) +
  geom_pointdensity(adjust = 4) +
  scale_color_viridis(option = "D") +
  annotate("text",
           x = max(data_filtered$rec_cmMb, na.rm = TRUE) * 0.8,
           y = max(data_filtered$mexicana_dosage, na.rm = TRUE) * 0.8,
           label = label_text, hjust = 1, vjust = 0, size = 4, color = "black") +
  theme_minimal() +
  labs(x = "Recombination rate (cM/Mb)",
       y = "Mexicana dosage",
       title = "Correlation between Recombination Rate and Mexicana Dosage",
       subtitle = "Using 0.5% - 99.5% quantile range for rec_cmMb") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        panel.grid = element_blank())

# Save correlation plot
out_pdf_corr <- file.path(dir_out, paste0(sample_id, ".CorrRecAdmix_99percentiles.pdf"))
out_png_corr <- file.path(dir_out, paste0(sample_id, ".CorrRecAdmix_99percentiles.png"))
ggsave(out_pdf_corr, p_corr, width = 9, height = 6)
ggsave(out_png_corr, p_corr, width = 9, height = 6, dpi = 300)
message("Saved correlation plots to: ", out_pdf_corr, " and ", out_png_corr)
}

# ---------------------------
# Decile boxplots
# ---------------------------
# Create decile breaks (0.1..0.9 quantiles -> 10 bins)
valid_for_deciles <- analysis_df  # using rec_cmMb not NA
deciles <- quantile(valid_for_deciles$rec_cmMb, probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE)
decile_breaks <- c(min(valid_for_deciles$rec_cmMb, na.rm = TRUE), deciles, max(valid_for_deciles$rec_cmMb, na.rm = TRUE))
decile_labels <- sapply(seq_len(length(decile_breaks)-1), function(i) {
  paste0(round(decile_breaks[i], 3), " - ", round(decile_breaks[i+1], 3))
})

valid_for_deciles$rec_decile <- cut(valid_for_deciles$rec_cmMb,
                                    breaks = decile_breaks,
                                    include.lowest = TRUE,
                                    labels = decile_labels)

# compute upper whisker to set y axis
upper_whiskers <- valid_for_deciles %>%
  group_by(rec_decile) %>%
  summarise(Q1 = quantile(mexicana_dosage, 0.25, na.rm = TRUE),
            Q3 = quantile(mexicana_dosage, 0.75, na.rm = TRUE),
            IQR = Q3 - Q1,
            upper_whisker = Q3 + 1.5 * IQR,
            .groups = "drop") %>%
  pull(upper_whisker)

y_max <- max(upper_whiskers, na.rm = TRUE)
if (!is.finite(y_max) || is.na(y_max)) y_max <- max(valid_for_deciles$mexicana_dosage, na.rm = TRUE)

p_decile_scatter <- ggplot(valid_for_deciles, aes(x = rec_decile, y = mexicana_dosage)) +
  geom_boxplot(outlier.shape = NA, fill = "coral") +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(x = "Recombination Rate Decile (cM/Mb)", y = "Mexicana Dosage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        panel.grid = element_blank())

# Save decile plot
out_pdf_decile <- file.path(dir_out, paste0(sample_id, ".DecileBoxplot_RecAdmix.pdf"))
out_png_decile <- file.path(dir_out, paste0(sample_id, ".DecileBoxplot_RecAdmix.png"))
ggsave(out_pdf_decile, p_decile_scatter, width = 9, height = 6)
ggsave(out_png_decile, p_decile_scatter, width = 9, height = 6, dpi = 300)
message("Saved decile plots to: ", out_pdf_decile, " and ", out_png_decile)

# End
message("All done.")
