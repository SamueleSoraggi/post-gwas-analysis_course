#!/usr/bin/env Rscript
# Check whether a built locus demonstrates the "effect of L" lesson:
#   L=1 should collapse to 1 credible set; L=10 should recover 2 (high-purity).
suppressMessages({library(data.table); library(susieR)})
a <- commandArgs(trailingOnly = TRUE); outdir <- a[1]; lab <- a[2]
ss <- fread(file.path(outdir, "sumstats.tsv"))
R  <- readRDS(file.path(outdir, "ld_matched.RDS"))
n  <- ss$N[1]
f1  <- suppressWarnings(susie_rss(z = ss$Z, R = R, n = n, L = 1))
f10 <- suppressWarnings(susie_rss(z = ss$Z, R = R, n = n, L = 10))
cs <- f10$sets$cs
lead_pos <- function(idx) ss$POS[idx[which.max(susie_get_pip(f10)[idx])]]
info <- if (length(cs)) paste(sprintf("CS%d(size=%d,purity=%.2f,pos=%d)",
          seq_along(cs), lengths(cs),
          f10$sets$purity[, "min.abs.corr"],
          vapply(cs, lead_pos, numeric(1))), collapse = "  ") else "none"
cat(sprintf("%-8s  L1_nCS=%d  L10_nCS=%d  | %s\n",
            lab, length(f1$sets$cs), length(cs), info))
