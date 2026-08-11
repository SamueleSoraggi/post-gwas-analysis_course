#!/usr/bin/env Rscript
# Purpose: Build a realistic example FLAMES output for the SORCS3 / ADHD locus,
#          in the REAL FLAMES output schema, for the teaching exercise. This is a
#          hand-crafted teaching stand-in (FLAMES needs a ~10 GB reference + the
#          genome-wide MAGMA/PoPS inputs, an HPC job, see flames_hpc/). It is NOT
#          a real FLAMES run, exactly like the SBayesRC mock in the fine-mapping
#          session. Numbers are plausible, not produced by the XGBoost model.
# Author:  course materials
# Date:    2026-06-17
# Usage:   Rscript make_flames_mock.R   (from session-flames/)
#
# Real schema (from FLAMES_scoring.py, v1.1.2):
#   .pred : locus filename symbol ensg FLAMES_scaled FLAMES_raw
#           estimated_cumulative_precision           (only FLAMES_causal genes)
#   .raw  : + XGB_score PoPS_Score rel_XGB_score rel_PoPS_Score
#           FLAMES_highest FLAMES_causal             (all scored genes)
# Gate (verified in the code): a gene is "causal"/reported only if FLAMES_raw > 0.134.
suppressMessages(library(data.table))

out_dir <- "../../input/flames/example_output/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Real genes in the SORCS3 locus (GRCh37 TSS from the FLAMES PoPS gene_annot).
# rs11596214 (~chr10:106.45 Mb) is INTRONIC in SORCS3, so SORCS3 is the nearest
# gene AND the gene Demontis et al. 2023 nominate (common + rare variant support).
genes <- data.table(
  symbol = c("SORCS3",          "CCDC147",        "ITPRIP",         "GSTO2",          "GSTO1",          "SFR1"),
  ensg   = c("ENSG00000156395", "ENSG00000120051","ENSG00000148841","ENSG00000065621","ENSG00000148834","ENSG00000156384"),
  tss    = c(106400859,         106113522,        106098162,        106028631,        105995114,        105881816)
)

# Per-gene evidence (teaching values). SORCS3 carries converging support across
# the XGBoost evidence ensemble (eQTL/chromatin/distance) AND the PoPS similarity
# score; neighbours have at most one weak stream.
genes[, XGB_score  := c(0.74, 0.18, 0.12, 0.09, 0.07, 0.04)]
genes[, PoPS_Score := c(1.93, 0.41, 0.22,-0.05, 0.10,-0.18)]

# Within-locus relative scores (FLAMES normalises each stream across the locus).
genes[, rel_XGB_score  := XGB_score / sum(XGB_score)]
pops_pos <- pmax(genes$PoPS_Score, 0)
genes[, rel_PoPS_Score := pops_pos / sum(pops_pos)]

# Combined FLAMES raw score: default weight 0.725 XGB / 0.275 PoPS.
w <- 0.725
genes[, FLAMES_raw    := round(w * rel_XGB_score + (1 - w) * rel_PoPS_Score, 4)]
genes[, FLAMES_scaled := round(FLAMES_raw / sum(FLAMES_raw), 4)]

# Gate + highest-gene flags (mirrors the real code: FLAMES_raw > 0.134).
genes[, FLAMES_highest := as.integer(FLAMES_raw == max(FLAMES_raw))]
genes[, FLAMES_causal  := as.integer(FLAMES_raw > 0.134 & FLAMES_highest == 1L)]
# Precision only for the reported (highest, above-gate) gene; 0 otherwise.
genes[, estimated_cumulative_precision := ifelse(FLAMES_causal == 1L, 0.86, 0)]

genes[, locus    := "1"]
genes[, filename := "locus_files/locus1_SORCS3.cred"]
setorder(genes, -FLAMES_scaled)

raw_cols <- c("locus","filename","symbol","ensg","XGB_score","PoPS_Score",
              "rel_XGB_score","rel_PoPS_Score","FLAMES_raw","FLAMES_scaled",
              "estimated_cumulative_precision","FLAMES_highest","FLAMES_causal")
fwrite(genes[, ..raw_cols], file.path(out_dir, "FLAMES_scores.raw"), sep = "\t")

pred_cols <- c("locus","filename","symbol","ensg","FLAMES_scaled","FLAMES_raw",
               "estimated_cumulative_precision")
fwrite(genes[FLAMES_causal == 1L, ..pred_cols],
       file.path(out_dir, "FLAMES_scores.pred"), sep = "\t")

cat("Wrote FLAMES_scores.pred (", sum(genes$FLAMES_causal), "reported gene) and",
    "FLAMES_scores.raw (", nrow(genes), "scored genes) to", out_dir, "\n")
