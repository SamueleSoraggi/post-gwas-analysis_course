#!/usr/bin/env Rscript
# Orchestrator: build all teaching loci for the fine-mapping exercise from real
# GIANT BMI summary statistics + real 1000G LD (EUR matched / EAS mismatched),
# then write a mock SBayesRC comparison for locus 1.
#
# Each locus is built by build_locus() in build_locus_1000g.R, which streams the
# 1000G region, computes signed LD with PLINK, and harmonises alleles so that z
# and R are mutually consistent. Region downloads are cached under data/genotypes/.
#
# Author:  course materials
# Date:    2026-06-10
# Usage:   pixi run prepare-data        (from session-finemapping/)
set.seed(42)
suppressMessages(library(data.table))

BASE <- Sys.getenv("PGC_COURSE_BASE",
                   unset = "/home/dipe/projects/projects/dimitris/post-gwas-course/session-finemapping")
source(file.path(BASE, "build_locus_1000g.R"))

# Teaching locus (ADHD; Demontis et al. 2023, Nat Genet)
# SORCS3 / rs11596214 (chr10:106.45 Mb): the paper's flagship gene, implicated by
# both common AND rare variants, reported as a SINGLE genome-wide signal (COJO
# found no robust second independent signal). One locus carries the whole exercise:
#   Block 1 fine-map it, Block 2 EUR/EAS LD mismatch, Block 3 the effect of L (and
#   a reality check that L / reference-LD can manufacture an extra credible set the
#   paper does not support), Block 4 SBayesRC + FLAMES export.
# NB: no separate "two-signal" locus. The sweep over the paper's multi-IndSigSNP
# loci (PTPRF, FOXP2, KAT2B, CDH8) showed none separate cleanly in SuSiE on 1000G
# EUR, consistent with the paper's COJO conclusion. See DECISIONS.md (2026-06-17).
# n_eff = 51,568 is the study's reported effective sample size (Demontis 2023).
# The daner Nca/Nco are pooled totals, so per-SNP 4/(1/Nca+1/Nco) would overestimate
# (~128k); we use the paper's figure to stay consistent with it.
build_locus(10, 105953832, 106953832, "sorcs3", file.path(BASE, "data/locus1_clean"),
            n_eff = 51568)

# Mock SBayesRC output for locus 1
# SBayesRC concentrates PIP more sharply than SuSiE; approximated here as a
# softmax over Z^2 / 15. This is a teaching stand-in, NOT a real SBayesRC run.
message("Building mock SBayesRC output for locus 1 ...")
ss1      <- fread(file.path(BASE, "data/locus1_clean", "sumstats.tsv"))
pip_sbrc <- exp(ss1$Z^2 / 15) / sum(exp(ss1$Z^2 / 15))

sbrc_dir <- file.path(BASE, "data/sbayesrc_precomputed")
dir.create(sbrc_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(
  data.table(
    Index = seq_len(nrow(ss1)),
    Name = ss1$SNP, Chrom = ss1$CHR, Position = ss1$POS, A1 = ss1$A1, A2 = ss1$A2,
    A1Frq = round(runif(nrow(ss1), 0.1, 0.9), 3),
    A1Effect = ss1$BETA, SE = ss1$SE,
    VarExplained = round(pip_sbrc * ss1$SE^2, 8),
    PEP = round(1 - pip_sbrc, 6),
    Pi1 = 0, Pi2 = 0, Pi3 = 0, Pi4 = round(pip_sbrc, 6),
    PIP = round(pip_sbrc, 6),
    GelmanRubin_R = 1.0
  ),
  file.path(sbrc_dir, "locus1.snpRes"), sep = "\t"
)
fwrite(
  data.table(Par = c("Pi","Sigma2","h2","NumSnps"),
             Estimate = c(0.00102, 3.18e-5, 0.274, nrow(ss1)),
             SE = c(0.00018, 8.2e-7, 0.012, NA_real_)),
  file.path(sbrc_dir, "locus1.parRes"), sep = "\t"
)
message("  SBayesRC mock files saved")
message("\nData preparation complete")
