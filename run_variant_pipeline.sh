#!/usr/bin/env bash
# Real variant-calling pipeline: samtools mpileup -> bcftools call, run
# against every synthetic sample's BAM produced by
# data/generate_synthetic_variant_data.py.
#
# Tools used (checked, not assumed): samtools (faidx, mpileup) and bcftools
# (call, norm, view) -- both were found installed in this environment
# (Homebrew). bedtools was NOT installed here and disk space was critically
# low (<1GB free), so region-level comparison against the gold-standard VCFs
# is instead done with a Python/pysam-free equivalent
# (compute_concordance.py), which parses both VCFs directly rather than
# shelling out to `bedtools intersect`. This is documented in the README.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$HERE/data/raw/reference.fasta"
BAM_DIR="$HERE/data/raw/bam"
OUT_DIR="$HERE/output/called_vcf"
mkdir -p "$OUT_DIR"

command -v samtools >/dev/null || { echo "samtools not found on PATH" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools not found on PATH" >&2; exit 1; }

if [ ! -f "$REF.fai" ]; then
  samtools faidx "$REF"
fi

echo "== Real variant calling pipeline (samtools mpileup | bcftools call) =="
n=0
for bam in "$BAM_DIR"/sample_*.sorted.bam; do
  sample=$(basename "$bam" .sorted.bam)
  vcf_out="$OUT_DIR/${sample}.called.vcf"
  bcftools mpileup -f "$REF" -A -Q 0 -a AD,DP "$bam" 2>/dev/null \
    | bcftools call -mv -Oz -o "$vcf_out.raw.gz" 2>/dev/null
  # Left-align/normalize indel representation against the reference so that
  # indel positions/alleles are comparable to the gold-standard VCF, which
  # is normalized the same way in generate_synthetic_variant_data.py. This
  # matters in practice: raw indel calls in repeat/homopolymer context can be
  # represented at more than one equally valid anchor position.
  bcftools norm -f "$REF" "$vcf_out.raw.gz" -Oz -o "$vcf_out.gz" 2>/dev/null
  bcftools view "$vcf_out.gz" > "$vcf_out"
  rm -f "$vcf_out.gz" "$vcf_out.raw.gz"
  n_calls=$(grep -vc "^#" "$vcf_out" || true)
  echo "[$sample] called variants: $n_calls -> $(basename "$vcf_out")"
  n=$((n+1))
done

echo "Pipeline complete: $n samples called. VCFs in $OUT_DIR"
