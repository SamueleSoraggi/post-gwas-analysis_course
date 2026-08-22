#!/usr/bin/env Rscript
# Build one fine-mapping locus from real data: ADHD sumstats + 1000G LD.
#   matched panel    = 1000G EUR  -> ld_matched.RDS
#   mismatched panel = 1000G EAS  -> ld_mismatched.RDS
# z is harmonised to the PLINK A1 allele convention so z and R are consistent.
#
# Use as a function (source this file):   build_locus(chr, start, end, label, outdir)
# Or from the shell (under `pixi run`):    Rscript build_locus_1000g.R <chr> <start> <end> <label> <outdir>
# Requires on PATH (provided by pixi): tabix, bgzip, plink.
#
# Author: course materials   Date: 2026-06-10
suppressMessages({library(data.table); library(susieR)})

# BASE/GWASF are overridable by env var for HPC (default = this repo + ADHD sumstats).
BASE   <- Sys.getenv("PGC_COURSE_BASE",
                     unset = "/home/dipe/projects/projects/dimitris/post-gwas-course/session-finemapping")
GWASF  <- Sys.getenv("PGC_COURSE_GWAS",
                     unset = file.path(BASE, "data/raw/ADHD2022_iPSYCH_deCODE_PGC.meta"))
URLTPL <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL  <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

# Read + harmonise a GWAS file to the common schema CHR POS SNP A1 A2 BETA SE Z N
# (the source P column is retained when present, e.g. for clumping). Two formats:
#   - PGC daner case/control (has OR, Nca, Nco): BETA = log(OR); N = effective N.
#   - GIANT-style additive (10 positional columns): BETA used directly.
# Returns the FULL table (no region subset) so a caller fine-mapping many loci can
# read the GWAS once, clump it, and pass the result to build_locus() via gwas_dt.
read_harmonised_gwas <- function(path = GWASF) {
  gwas <- fread(path)
  if ("OR" %in% names(gwas)) {
    # PGC daner case/control format (ADHD2022): OR is on A1; SE is on the log-OR
    # scale; per-SNP case/control counts Nca/Nco give the effective sample size.
    setnames(gwas, "BP", "POS", skip_absent = TRUE)
    gwas[, BETA := log(OR)]
    gwas[, Z    := BETA / SE]
    gwas[, N    := 4 / (1 / Nca + 1 / Nco)]   # effective N for case/control
  } else {
    # GIANT-style additive quantitative format (10 positional columns)
    setnames(gwas, c("CHR","POS","SNP","A1","A2","FREQ","BETA","SE","P","N"))
    gwas[, Z := BETA / SE]
  }
  gwas[]
}

# n_eff: optional locus-level effective sample size to write into sumstats.tsv.
#   For a multi-cohort case/control meta the per-SNP 4/(1/Nca+1/Nco) on POOLED
#   counts overestimates the true effective N, so pass the study's reported value
#   (e.g. ADHD neff = 51,568). If NULL, the median per-SNP value is used.
# gwas_dt: optional pre-read harmonised GWAS (from read_harmonised_gwas) to avoid
#   re-reading the (large) sumstats file when building many loci in one run.
build_locus <- function(chr, start, end, label, outdir, n_eff = NULL, gwas_dt = NULL) {
  wd <- file.path(BASE, "data/genotypes", label); dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  sh <- function(cmd, args) stopifnot(system2(cmd, args, stdout = FALSE, stderr = FALSE) == 0L)
  p  <- function(...) file.path(wd, ...)

  # OFFLINE mode: if G1000_EUR points to a local 1000G EUR plink panel, build EUR
  # LD straight from it (no EBI streaming, no EAS). This is the genome-wide
  # research path (run_susie_all_loci.R exports G1000_EUR). When G1000_EUR is
  # unset the original streaming + EUR/EAS course path runs unchanged.
  LOCAL_PANEL <- Sys.getenv("G1000_EUR", unset = "")
  local_mode  <- nzchar(LOCAL_PANEL) && file.exists(paste0(LOCAL_PANEL, ".bed"))

  if (local_mode) {
    # region straight from the local EUR panel (already EUR -> no sample keep, no EAS)
    sh("plink", c("--bfile", LOCAL_PANEL, "--chr", chr, "--from-bp", start, "--to-bp", end,
                  "--keep-allele-order", "--biallelic-only", "strict", "--snps-only", "just-acgt",
                  "--set-missing-var-ids", "@:#:$1:$2", "--make-bed", "--out", p("region")))
    sh("plink", c("--bfile", p("region"), "--keep-allele-order", "--freq", "--out", p("freq_eur")))
  } else {
    # 1. sample lists (cached) --------------------------------------------------
    if (!file.exists(p("eur.keep"))) {
      panel <- fread(PANEL)
      fwrite(panel[super_pop == "EUR", .(sample, sample)], p("eur.keep"), sep = "\t", col.names = FALSE)
      fwrite(panel[super_pop == "EAS", .(sample, sample)], p("eas.keep"), sep = "\t", col.names = FALSE)
    }
    # 2. region VCF (cached) ----------------------------------------------------
    if (!file.exists(p("region.vcf.gz"))) {
      message(sprintf("[%s] streaming %d:%d-%d from 1000G ...", label, chr, start, end))
      url <- sprintf(URLTPL, chr)
      # run inside the (gitignored) cache dir so tabix's downloaded .tbi index lands there
      stopifnot(system(sprintf("cd '%s' && tabix -h '%s' %d:%d-%d | bgzip > region.vcf.gz",
                               wd, url, chr, start, end)) == 0L)
    }
    # 3. import + per-pop freq --------------------------------------------------
    sh("plink", c("--vcf", p("region.vcf.gz"), "--double-id", "--keep-allele-order",
                  "--biallelic-only", "strict", "--snps-only", "just-acgt",
                  "--set-missing-var-ids", "@:#:$1:$2", "--make-bed", "--out", p("region")))
    sh("plink", c("--bfile", p("region"), "--keep", p("eur.keep"), "--keep-allele-order", "--freq", "--out", p("freq_eur")))
    sh("plink", c("--bfile", p("region"), "--keep", p("eas.keep"), "--keep-allele-order", "--freq", "--out", p("freq_eas")))
  }

  # 4. GWAS window + SNP selection ---------------------------------------------
  # Harmonise to a common schema: CHR POS SNP A1 A2 BETA SE Z N(effective).
  # A1 is the effect allele in both supported formats. Re-use a pre-read table
  # (gwas_dt) when fine-mapping many loci, else read + harmonise the file here.
  gwas <- if (is.null(gwas_dt)) read_harmonised_gwas(GWASF) else gwas_dt
  g <- gwas[CHR == chr & POS >= start & POS <= end][order(POS)][!duplicated(POS)]

  bim <- fread(p("region.bim"), header = FALSE, col.names = c("chr","id","cm","pos","a1","a2"))
  fe  <- fread(p("freq_eur.frq"))[, .(id = SNP, maf_e = MAF)]
  bim <- merge(bim, fe, by = "id")
  if (!local_mode) {
    fa  <- fread(p("freq_eas.frq"))[, .(id = SNP, maf_a = MAF)]
    bim <- merge(bim, fa, by = "id")
  }

  m   <- merge(bim, g, by.x = "pos", by.y = "POS")
  amb <- function(x,y) (x=="A"&y=="T")|(x=="T"&y=="A")|(x=="C"&y=="G")|(x=="G"&y=="C")
  ok  <- function(A1,A2,a1,a2) (A1==a1&A2==a2)|(A1==a2&A2==a1)
  maf_ok <- if (local_mode) m$maf_e > 0.01 else (m$maf_e > 0.01 & m$maf_a > 0.01)
  m   <- m[ok(A1,A2,a1,a2) & !amb(A1,A2) & maf_ok][order(pos)]
  writeLines(m$id, p("snps.keep"))

  # 5. subset + signed LD (EUR; +EAS only in streaming course mode) ------------
  sh("plink", c("--bfile", p("region"), "--extract", p("snps.keep"), "--keep-allele-order", "--make-bed", "--out", p("region_sub")))
  if (local_mode) {
    sh("plink", c("--bfile", p("region_sub"), "--keep-allele-order", "--r", "square", "gz", "--out", p("ld_eur")))
  } else {
    sh("plink", c("--bfile", p("region_sub"), "--keep", p("eur.keep"), "--keep-allele-order", "--r", "square", "gz", "--out", p("ld_eur")))
    sh("plink", c("--bfile", p("region_sub"), "--keep", p("eas.keep"), "--keep-allele-order", "--r", "square", "gz", "--out", p("ld_eas")))
  }

  # 6. harmonise z to plink A1; write outputs ----------------------------------
  bs   <- fread(p("region_sub.bim"), header = FALSE, col.names = c("chr","id","cm","pos","a1","a2"))
  mm   <- m[match(bs$id, id)]
  stopifnot(nrow(mm) == nrow(bs), !anyNA(mm$Z))
  flip <- ifelse(mm$A1 == mm$a1, 1, -1)
  # one locus-level (effective) N so the exercise's `sumstats$N[1]` is stable
  N_locus <- if (!is.null(n_eff)) round(n_eff) else round(median(mm$N))
  sumstats <- data.table(CHR = chr, POS = mm$pos, SNP = mm$SNP, A1 = mm$a1, A2 = mm$a2,
                         BETA = mm$BETA * flip, SE = mm$SE, Z = mm$Z * flip, N = N_locus)
  fwrite(sumstats, file.path(outdir, "sumstats.tsv"), sep = "\t")

  read_ld <- function(f, ids) { R <- as.matrix(fread(f)); R[is.na(R)] <- 0; diag(R) <- 1; dimnames(R) <- list(ids, ids); R }
  R_eur <- read_ld(p("ld_eur.ld.gz"), bs$id); saveRDS(R_eur, file.path(outdir, "ld_matched.RDS"))
  if (!local_mode) { R_eas <- read_ld(p("ld_eas.ld.gz"), bs$id); saveRDS(R_eas, file.path(outdir, "ld_mismatched.RDS")) }

  # 7. verify -------------------------------------------------------------------
  n   <- sumstats$N[1]
  fe2 <- suppressWarnings(susie_rss(z = sumstats$Z, R = R_eur, n = n, L = 10))
  if (local_mode) {
    cat(sprintf("[%s] chr%d:%d-%d  SNPs:%d  lead|Z|:%.1f  s_hat EUR:%.3f  nCS EUR:%d  topPIP:%.2f\n",
                label, chr, start, end, nrow(sumstats), max(abs(sumstats$Z)),
                estimate_s_rss(sumstats$Z, R_eur, n = n),
                length(fe2$sets$cs), max(susie_get_pip(fe2))))
  } else {
    fa2 <- suppressWarnings(susie_rss(z = sumstats$Z, R = R_eas, n = n, L = 10))
    cat(sprintf("[%s] chr%d:%d-%d  SNPs:%d  lead|Z|:%.1f  s_hat EUR:%.3f EAS:%.3f  nCS EUR:%d EAS:%d  topPIP:%.2f\n",
                label, chr, start, end, nrow(sumstats), max(abs(sumstats$Z)),
                estimate_s_rss(sumstats$Z, R_eur, n = n), estimate_s_rss(sumstats$Z, R_eas, n = n),
                length(fe2$sets$cs), length(fa2$sets$cs), max(susie_get_pip(fe2))))
  }
  invisible(sumstats)
}

# CLI entry point (skipped when this file is sourced with no args) --------------
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 5)
  build_locus(as.integer(.args[1]), as.integer(.args[2]), as.integer(.args[3]), .args[4], .args[5])
