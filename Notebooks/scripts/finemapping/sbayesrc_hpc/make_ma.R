#!/usr/bin/env Rscript
# Purpose: Convert the PGC ADHD daner summary statistics into the GCTB .ma format
#          required by SBayesRC. Runs on the HPC only (it reads the raw, DUA-bound
#          GWAS file). It NEVER prints data rows, it writes straight to disk.
# Author:  course materials
# Date:    2026-06-17
# Usage:   Rscript make_ma.R <daner.meta> <out.ma> [N_eff]
#
# .ma columns (GCTB): SNP A1 A2 freq b se p N
#   daner schema: CHR SNP BP A1 A2 FRQ_A_* FRQ_U_* INFO OR SE P Direction Nca Nco
#   A1 is the effect allele; b = log(OR); freq = case/control-weighted A1 freq.
#   N is the study's reported EFFECTIVE sample size (default 51568, Demontis 2023),
#   used as a constant. The daner Nca/Nco are pooled totals, so per-SNP
#   4/(1/Nca+1/Nco) (~128k) overestimates the true meta-analytic effective N;
#   the paper's figure keeps this consistent with the locus exercise.
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("Usage: Rscript make_ma.R <daner.meta> <out.ma> [N_eff]")
daner_file <- args[[1]]
out_file   <- args[[2]]
n_eff      <- if (length(args) >= 3L) as.numeric(args[[3]]) else 51568

g <- fread(daner_file)

frq_a <- grep("^FRQ_A", names(g), value = TRUE)[1]
frq_u <- grep("^FRQ_U", names(g), value = TRUE)[1]
stopifnot(!is.na(frq_a), !is.na(frq_u), all(c("OR","SE","P","Nca","Nco") %in% names(g)))

# Case/control-weighted A1 frequency (population proxy).
freq <- (g$Nca * g[[frq_a]] + g$Nco * g[[frq_u]]) / (g$Nca + g$Nco)

ma <- data.table(
  SNP  = g$SNP,
  A1   = g$A1,
  A2   = g$A2,
  freq = round(freq, 6),
  b    = round(log(g$OR), 6),
  se   = g$SE,
  p    = g$P,
  N    = round(n_eff)
)
ma <- ma[is.finite(b) & is.finite(se) & se > 0]

fwrite(ma, out_file, sep = " ")
cat("Wrote", nrow(ma), "SNPs to", out_file, "(N =", round(n_eff), ")\n")
