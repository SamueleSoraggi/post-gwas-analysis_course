#!/usr/bin/env Rscript
# Purpose: Take the genome-wide SBayesRC output (.snpRes) and write the per-locus
#          file the fine-mapping exercise's Block 4 reads, subset to exactly the
#          SNPs in the SORCS3 teaching locus and mapped to the column names the
#          .qmd expects (Name Chrom Position A1 A2 A1Freq b se p pip Var).
# Author:  course materials
# Date:    2026-06-17
# Usage:   Rscript subset_locus.R <genomewide.snpRes> <locus_sumstats.tsv> <out.snpRes>
#
# Real SBayesRC .snpRes columns:
#   Index Name Chrom Position A1 A2 A1Frq A1Effect SE VarExplained PEP PIP GelmanRubin_R
# (no per-SNP p-value column, and effect/SE/PIP are capitalised), so we rename
# A1Effect->b, SE->se, PIP->pip, A1Frq->A1Freq, VarExplained->Var and recover p
# from the locus Z-scores.
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L)
  stop("Usage: Rscript subset_locus.R <genomewide.snpRes> <locus_sumstats.tsv> <out.snpRes>")
gw_file    <- args[[1]]
locus_file <- args[[2]]
out_file   <- args[[3]]

gw    <- fread(gw_file)
locus <- fread(locus_file)   # data/locus1_clean/sumstats.tsv: CHR POS SNP A1 A2 BETA SE Z N

# Keep only the SNPs used in the locus exercise, matched by name, in locus order
# so Block 4's merge with the SuSiE PIPs lines up.
sub <- gw[match(locus$SNP, Name)]
keep <- !is.na(sub$Name)
sub  <- sub[keep]
locus <- locus[keep]
if (nrow(sub) == 0L) stop("No overlap between genome-wide .snpRes and the locus SNPs")

out <- data.table(
  Name     = sub$Name,
  Chrom    = sub$Chrom,
  Position = sub$Position,
  A1       = sub$A1,
  A2       = sub$A2,
  A1Freq   = sub$A1Frq,
  b        = sub$A1Effect,
  se       = sub$SE,
  p        = signif(2 * pnorm(-abs(locus$Z)), 3),   # recovered from the GWAS Z
  PIP      = sub$PIP,
  Var      = sub$VarExplained
)
fwrite(out, out_file, sep = "\t")
cat("Wrote", nrow(out), "of", nrow(locus), "locus SNPs to", out_file, "\n")
