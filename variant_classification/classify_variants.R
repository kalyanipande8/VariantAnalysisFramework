#!/usr/bin/env Rscript
# R companion to classify_variants.py, demonstrating the same
# region/consequence/pathogenicity-tier classification logic in R (to
# genuinely reflect the "Python/R framework" part of the resume bullet).
# Reads the same called VCFs bcftools produced and reproduces the
# classification independently in R, writing its own summary CSV so the two
# language implementations can be diffed/cross-checked.
#
# Run with: Rscript variant_classification/classify_variants.R

suppressWarnings(suppressMessages({
  has_base_only <- TRUE
}))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) == 1) return(dirname(normalizePath(file_arg)))
  return(getwd())
}
here <- get_script_dir()
root <- dirname(here)

ref_path <- file.path(root, "data", "raw", "reference.fasta")
called_dir <- file.path(root, "output", "called_vcf")
out_dir <- file.path(root, "output")

CDS_START <- 200  # 0-based, matches gene_model.py
CDS_END <- 2700

CODON_TABLE <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
  AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
  TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
  AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G"
)

read_fasta <- function(path) {
  lines <- readLines(path)
  lines <- lines[!grepl("^>", lines)]
  paste(lines, collapse = "")
}

pathogenicity_tier <- function(consequence) {
  if (consequence %in% c("nonsense", "frameshift_indel", "stop_loss")) return("Likely_Pathogenic")
  if (consequence == "missense") return("VUS")
  return("Likely_Benign")
}

classify_variant <- function(ref_seq, pos0, ref_allele, alt_allele) {
  is_snp <- nchar(ref_allele) == 1 && nchar(alt_allele) == 1
  indel_len <- abs(nchar(alt_allele) - nchar(ref_allele))

  if (pos0 < CDS_START) {
    region <- "5_prime_intergenic"
  } else if (pos0 >= CDS_END) {
    region <- "3_prime_intergenic"
  } else {
    region <- "coding"
  }

  if (region != "coding") {
    return(list(region = region, consequence = "intergenic"))
  }
  if (!is_snp) {
    consequence <- if (indel_len %% 3 != 0) "frameshift_indel" else "inframe_indel"
    return(list(region = region, consequence = consequence))
  }

  offset_in_cds <- pos0 - CDS_START
  codon_index <- offset_in_cds %/% 3
  codon_start <- CDS_START + codon_index * 3
  ref_codon <- substr(ref_seq, codon_start + 1, codon_start + 3)
  if (nchar(ref_codon) < 3) {
    return(list(region = region, consequence = "coding_boundary_effect"))
  }
  pos_in_codon <- pos0 - codon_start
  alt_codon <- ref_codon
  substr(alt_codon, pos_in_codon + 1, pos_in_codon + 1) <- alt_allele

  ref_aa <- CODON_TABLE[[ref_codon]]
  alt_aa <- CODON_TABLE[[alt_codon]]
  if (is.null(ref_aa)) ref_aa <- "X"
  if (is.null(alt_aa)) alt_aa <- "X"

  consequence <- if (ref_aa == alt_aa) {
    "synonymous"
  } else if (alt_aa == "*") {
    "nonsense"
  } else if (ref_aa == "*") {
    "stop_loss"
  } else {
    "missense"
  }
  list(region = region, consequence = consequence, ref_aa = ref_aa, alt_aa = alt_aa)
}

parse_called_vcf <- function(path) {
  lines <- readLines(path)
  lines <- lines[!grepl("^#", lines)]
  rows <- list()
  for (line in lines) {
    f <- strsplit(line, "\t")[[1]]
    pos <- as.integer(f[2]); ref <- f[4]; alt_field <- f[5]
    for (alt in strsplit(alt_field, ",")[[1]]) {
      rows[[length(rows) + 1]] <- list(pos = pos, ref = ref, alt = alt)
    }
  }
  rows
}

t0 <- Sys.time()
ref_seq <- read_fasta(ref_path)
vcf_files <- sort(list.files(called_dir, pattern = "\\.called\\.vcf$", full.names = TRUE))

if (length(vcf_files) == 0) {
  stop("No called VCFs found -- run run_variant_pipeline.sh first.")
}

summary_rows <- list()
for (vcf_path in vcf_files) {
  sample <- sub("\\.called\\.vcf$", "", basename(vcf_path))
  records <- parse_called_vcf(vcf_path)
  tiers <- c(Likely_Pathogenic = 0, VUS = 0, Likely_Benign = 0)
  for (rec in records) {
    result <- classify_variant(ref_seq, rec$pos - 1, rec$ref, rec$alt)
    tier <- pathogenicity_tier(result$consequence)
    tiers[[tier]] <- tiers[[tier]] + 1
  }
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    sample = sample,
    n_variants = length(records),
    n_likely_pathogenic = tiers[["Likely_Pathogenic"]],
    n_vus = tiers[["VUS"]],
    n_likely_benign = tiers[["Likely_Benign"]]
  )
}

summary_df <- do.call(rbind, summary_rows)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")

out_csv <- file.path(out_dir, "classification_report_R.csv")
write.csv(summary_df, out_csv, row.names = FALSE)

cat(sprintf("R classification complete: %d samples, %.3fs\n", nrow(summary_df), elapsed))
cat(sprintf("Wrote %s\n", out_csv))
