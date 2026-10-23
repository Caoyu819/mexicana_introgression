#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript step1_add_rec_per_chr.R <f_admix> <f_rec> <f_out> <chr_name>")
}

f_admix  <- args[1]
f_rec    <- args[2]
f_out    <- args[3]
chr_name <- args[4]

if (!file.exists(f_admix)) stop("Admixture file not found: ", f_admix)
if (!file.exists(f_rec))   stop("Recombination map file not found: ", f_rec)

dir.create(dirname(f_out), recursive = TRUE, showWarnings = FALSE)

# Helper: Match first existing column from candidates
get_col_name <- function(dt, candidates) {
  for (nm in candidates) {
    if (nm %in% names(dt)) return(nm)
  }
  return(NULL)
}

# Helper: Normalize chromosome identifiers (e.g., "chr1", "Chr1", "1" -> "1")
norm_chr <- function(x) {
  sub("^(chr|chrom|chromosome)", "", trimws(as.character(x)), ignore.case = TRUE)
}

target_chr <- norm_chr(chr_name)

# 1. Process recombination map
rec_dt  <- fread(f_rec)
pos_col <- get_col_name(rec_dt, c("Position(bp)", "POS", "Position", "pos", "bp"))
cm_col  <- get_col_name(rec_dt, c("Map(cM)", "cM", "Map", "cm"))
chr_col <- get_col_name(rec_dt, c("Chr", "chr", "#CHROM", "CHROM", "Chromosome"))

if (is.null(pos_col) || is.null(cm_col)) {
  stop("Columns 'Position(bp)' or 'Map(cM)' not found in: ", f_rec)
}

# Filter by chromosome
if (!is.null(chr_col)) {
  rec_dt <- rec_dt[norm_chr(get(chr_col)) == target_chr]
}

# Remove NA and cast numeric
rec_dt <- rec_dt[!is.na(get(pos_col)) & !is.na(get(cm_col)),
                 .(POS = as.numeric(get(pos_col)), cM = as.numeric(get(cm_col)))]

# Aggregate duplicated POS caused by liftover by taking the mean cM
rec_dt <- rec_dt[, .(cM = mean(cM)), by = POS][order(POS)]

if (nrow(rec_dt) < 2) {
  stop("Fewer than 2 valid markers found for chromosome ", chr_name, " in: ", f_rec)
}

# Apply Isotonic Regression to strictly enforce monotonic non-decreasing cM
# This resolves local coordinate inversions caused by liftover (PAVA algorithm)
iso_fit   <- isoreg(x = rec_dt$POS, y = rec_dt$cM)
rec_dt[, cM := iso_fit$yf]

min_rec_pos <- min(rec_dt$POS)
max_rec_pos <- max(rec_dt$POS)
max_rec_cm  <- max(rec_dt$cM)
interp_cm   <- approxfun(x = rec_dt$POS, y = rec_dt$cM, rule = 2)

# 2. Process admixture data and calculate recombination rate
admixmap <- fread(f_admix)
mex_col  <- get_col_name(admixmap, c("mexicana_dosage", "mexicana", "Mean", "mean", "dosage", "All_mean"))
adm_chr  <- get_col_name(admixmap, c("#CHROM", "CHROM", "Chr", "chr", "Chromosome"))

if (is.null(mex_col)) stop("Dosage column not found in: ", f_admix)

if (!is.null(adm_chr)) {
  admixmap <- admixmap[norm_chr(get(adm_chr)) == target_chr]
}

if (nrow(admixmap) == 0) {
  stop("No records found for chromosome ", chr_name, " in: ", f_admix)
}

admixmap[, mexicana_dosage := as.numeric(get(mex_col))]
admixmap[, POS := as.numeric(POS)]
admixmap <- admixmap[order(POS)]

# Interpolate cM
admixmap[, total_cm := interp_cm(POS)]
admixmap[POS < min_rec_pos, total_cm := 0]
admixmap[POS > max_rec_pos, total_cm := max_rec_cm]

# Calculate stepwise differences
admixmap[, diff_cm  := c(NA_real_, diff(total_cm))]
admixmap[, diff_bp  := c(NA_real_, diff(POS))]
admixmap[, rec_cmMb := as.numeric(diff_cm / (diff_bp / 1e6))]
admixmap[!is.finite(rec_cmMb), rec_cmMb := NA_real_]

# Enforce boundary rules outside map coverage
admixmap[POS < min_rec_pos, `:=`(total_cm = 0, rec_cmMb = 0)]
admixmap[POS > max_rec_pos, `:=`(total_cm = max_rec_cm, rec_cmMb = 0)]
admixmap[, Chr := chr_name]

# 3. Export chromosome output
fwrite(admixmap[, !c("diff_cm", "diff_bp")], file = f_out, sep = "\t", quote = FALSE)
message("Successfully wrote: ", f_out)
