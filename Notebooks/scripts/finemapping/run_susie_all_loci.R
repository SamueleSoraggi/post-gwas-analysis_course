#!/usr/bin/env Rscript
# Purpose:  Fine-map ALL genome-wide significant ADHD loci with SuSiE and export
#           each credible set in FLAMES input format, so a downstream FLAMES run
#           covers every locus (not just the SORCS3 teaching locus).
# Author:   course materials
# Date:     2026-06-22
# Usage:    export G1000_EUR=/path/to/g1000_eur     # genome-wide 1000G EUR plink
#           pixi run Rscript run_susie_all_loci.R   # from session-finemapping/
# Requires: plink, tabix, bgzip on PATH (pixi provides them); internet for 1000G.
# Inputs:   raw ADHD sumstats (env PGC_COURSE_GWAS or data/raw/ADHD2022_...meta)
#           genome-wide 1000G EUR plink fileset for clumping (env G1000_EUR; the
#           .bim variant IDs MUST be rsIDs matching the GWAS SNP column).
# Outputs:  results/all_loci/      per-locus sumstats + LD + plink clump log
#           results/flames_input/  one CS file per credible set + flames_index.tsv
#
# FLAMES contract (verified against tools/FLAMES v1.1.2):
#   * credset file: whitespace-delim WITH header; columns named by `-sc`/`-pc`
#     (we write `SNP`<TAB>`PIP`). SNP id = CHR:POS:A1:A2, CHR/POS bare integers.
#   * index file (`-id`): tab-delim WITH header; REQUIRES a `Filename` column
#     (path FLAMES opens) + optional `GenomicLocus`. Filename must be findable
#     from FLAMES' cwd, so we write absolute paths.
set.seed(42)
suppressMessages({library(data.table); library(susieR)})

BASE <- Sys.getenv("PGC_COURSE_BASE", unset = getwd())
source(file.path(BASE, "build_locus_1000g.R"))

GWASF     <- Sys.getenv("PGC_COURSE_GWAS",
                        unset = file.path(BASE, "data/raw/ADHD2022_iPSYCH_deCODE_PGC.meta"))
G1000_EUR <- Sys.getenv("G1000_EUR", unset = "")
N_EFF    <- 51568
WINDOW   <- 500000   # +/-500kb around each lead SNP
P_THRESH <- 5e-8
MAX_LOCI <- as.integer(Sys.getenv("MAX_LOCI", unset = "0"))  # >0 = test mode: top-N loci only

results_dir <- file.path(BASE, "results/all_loci")
flames_dir  <- file.path(BASE, "results/flames_input")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(flames_dir,  showWarnings = FALSE, recursive = TRUE)

if (G1000_EUR == "" || !file.exists(paste0(G1000_EUR, ".bed")))
  stop("Set G1000_EUR to a genome-wide 1000G EUR plink prefix (rsID-keyed). ",
       "Got '", G1000_EUR, "' (no .bed found).")

# -- 1. Read + harmonise GWAS once --------------------------------------------
message("Reading + harmonising GWAS ...")
gwas <- read_harmonised_gwas(GWASF)
gwas[, P := as.numeric(P)]

# -- 2. Clump to find independent lead SNPs -----------------------------------
clump_input <- file.path(results_dir, "gwas_for_clump.tsv")
fwrite(gwas[!is.na(P), .(SNP, P)], clump_input, sep = "\t")

# Catch the classic rsID-vs-chr:pos panel mismatch before it silently yields 0 clumps.
bim_ids <- fread(paste0(G1000_EUR, ".bim"), header = FALSE, select = 2L)[[1]]
n_overlap <- sum(gwas$SNP %in% bim_ids)
message(sprintf("  %d / %d GWAS SNP IDs match the panel .bim", n_overlap, nrow(gwas)))
if (n_overlap == 0)
  stop("No GWAS SNP IDs match the panel .bim. The panel is probably keyed on ",
       "chr:pos, not rsIDs. Clumping cannot match anything.")

message("Clumping ...")
clump_out <- file.path(results_dir, "clumped")
clump_log <- file.path(results_dir, "plink_clump.log")
rc <- system2("plink", c(
  "--bfile", G1000_EUR,
  "--clump", clump_input,
  "--clump-snp-field", "SNP", "--clump-field", "P",
  "--clump-p1", format(P_THRESH, scientific = FALSE),
  "--clump-r2", "0.1",
  "--clump-kb", "500",
  "--out", clump_out
), stdout = clump_log, stderr = clump_log)
if (rc != 0L)
  stop("plink --clump returned ", rc, ". See ", clump_log)
if (!file.exists(paste0(clump_out, ".clumped")))
  stop("plink wrote no .clumped file (no SNP passed P<", format(P_THRESH, scientific = FALSE),
       ", or none matched the panel). See ", clump_log)

# read.table collapses plink's space-padded columns robustly (fread can mis-split them)
clumped <- as.data.table(read.table(paste0(clump_out, ".clumped"),
                                    header = TRUE, stringsAsFactors = FALSE))
clumped <- clumped[CHR %in% 1:22]
if (MAX_LOCI > 0L) {
  clumped <- head(clumped[order(P)], MAX_LOCI)   # test mode: top-N most significant loci
  message(sprintf("[MAX_LOCI=%d] TEST MODE: fine-mapping only the top %d loci", MAX_LOCI, nrow(clumped)))
}
leads <- clumped[, .(SNP, CHR = as.integer(CHR), POS = as.integer(BP))][order(CHR, POS)]
message(sprintf("Found %d independent autosomal lead SNPs (P < %g)", nrow(leads), P_THRESH))

# -- 3. Fine-map each locus + export credible sets ----------------------------
index_rows <- list()
gl <- 0L   # GenomicLocus counter (one per exported credible set)

for (i in seq_len(nrow(leads))) {
  snp   <- leads$SNP[i]
  chr   <- leads$CHR[i]
  pos   <- leads$POS[i]
  start <- max(1, pos - WINDOW)
  end   <- pos + WINDOW
  label <- sprintf("locus%03d_chr%d_%s", i, chr, snp)
  locus_dir <- file.path(results_dir, label)

  message(sprintf("\n-- [%d/%d] %s  chr%d:%d-%d --", i, nrow(leads), snp, chr, start, end))

  tryCatch({
    build_locus(chr, start, end, label, locus_dir, n_eff = N_EFF, gwas_dt = gwas)

    ss <- fread(file.path(locus_dir, "sumstats.tsv"))
    ld <- readRDS(file.path(locus_dir, "ld_matched.RDS"))

    fit <- susie_rss(z = ss$Z, R = ld, n = ss$N[1], L = 10)

    if (length(fit$sets$cs) == 0) {
      message(sprintf("  No credible sets (s_hat = %.3f)",
                      estimate_s_rss(ss$Z, ld, n = ss$N[1])))
      next
    }

    pip <- susie_get_pip(fit)
    for (cs_i in seq_along(fit$sets$cs)) {
      idx     <- fit$sets$cs[[cs_i]]
      cs_snps <- ss[idx]
      cs_dt   <- data.table(
        SNP = paste0(cs_snps$CHR, ":", cs_snps$POS, ":", cs_snps$A1, ":", cs_snps$A2),
        PIP = pip[idx]
      )
      cs_file <- file.path(flames_dir, sprintf("%s_CS%d.tsv", label, cs_i))
      fwrite(cs_dt, cs_file, sep = "\t")   # WITH header -> FLAMES -sc SNP -pc PIP

      gl <- gl + 1L
      index_rows[[length(index_rows) + 1]] <- data.table(
        Filename     = normalizePath(cs_file),
        GenomicLocus = gl
      )
      message(sprintf("  CS%d: %d SNPs, purity = %.3f, topPIP = %.3f -> %s",
                      cs_i, length(idx), fit$sets$purity[cs_i, "min.abs.corr"],
                      max(cs_dt$PIP), basename(cs_file)))
    }
  }, error = function(e) {
    message(sprintf("  FAILED: %s", conditionMessage(e)))
  })
}

# -- 4. Write the FLAMES index file -------------------------------------------
if (length(index_rows) > 0) {
  index_dt   <- rbindlist(index_rows)
  index_file <- file.path(flames_dir, "flames_index.tsv")
  fwrite(index_dt, index_file, sep = "\t")
  message(sprintf("\nDone: %d credible sets from %d loci.\n  Credible sets: %s\n  FLAMES index:  %s\n  -> feed to: FLAMES.py annotate -id %s -sc SNP -pc PIP ...",
                  nrow(index_dt), nrow(leads), flames_dir, index_file, index_file))
} else {
  message("\nNo credible sets found at any locus.")
}
