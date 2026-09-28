#!/usr/bin/env python3
"""
Compute real concordance / correlation metrics between the pipeline's called
variants (output/called_vcf/*.called.vcf, produced by run_variant_pipeline.sh
via samtools mpileup | bcftools call) and the gold-standard truth VCFs
(data/raw/gold_vcf/*.gold.vcf, produced by generate_synthetic_variant_data.py).

Metric definitions (stated explicitly, as required for an honest result):

1. Per-sample genotype concordance rate:
     concordance = |called_positions ∩ gold_positions with matching REF/ALT|
                   / |called_positions ∪ gold_positions|
   i.e. an (intersection-over-union) style rate over the set of variant
   sites, analogous to what `bedtools intersect` / `bedtools jaccard` would
   report at the position level. bedtools itself was not installed in this
   environment (and disk space was too tight to install it), so this
   Python-only equivalent parses both VCFs directly -- documented in README.

2. The headline "correlation with gold-standard calls" figure quoted in the
   resume bullet is defined here as the Pearson correlation coefficient,
   computed across every variant site observed in the union of gold-standard
   and called sites, pooled over all 25 samples, between:
     - true_vaf: the gold-standard's simulated variant allele fraction
       (0.0 for sites the caller found that are NOT in the gold standard --
       i.e. false positives)
     - called_vaf: the pipeline's called alt-allele fraction from the
       bcftools AD (allele depth) FORMAT field (0.0 for gold-standard sites
       the caller MISSED -- i.e. false negatives)
   This is a real, defensible quantitative concordance metric: it
   simultaneously penalizes false positives, false negatives, and VAF
   estimation error, rather than only counting exact position matches.
"""
import csv
import glob
import os
import re
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
GOLD_DIR = os.path.join(HERE, "data", "raw", "gold_vcf")
CALLED_DIR = os.path.join(HERE, "output", "called_vcf")


def parse_gold_vcf(path):
    """Return dict {(pos, ref, alt): true_vaf}."""
    variants = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos, ref, alt, info = fields[1], fields[3], fields[4], fields[7]
            vaf_match = re.search(r"VAF=([\d.]+)", info)
            vaf = float(vaf_match.group(1)) if vaf_match else 1.0
            variants[(int(pos), ref, alt)] = vaf
    return variants


def parse_called_vcf(path):
    """Return dict {(pos, ref, alt): called_vaf} using AD FORMAT field."""
    variants = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos, ref, alt_field = fields[1], fields[3], fields[4]
            fmt_keys = fields[8].split(":")
            fmt_vals = fields[9].split(":")
            fmt = dict(zip(fmt_keys, fmt_vals))
            for alt in alt_field.split(","):
                if alt == "<*>" or alt == ".":
                    continue
                vaf = None
                if "AD" in fmt:
                    ad_parts = fmt["AD"].split(",")
                    try:
                        ad_vals = [int(x) for x in ad_parts]
                        total = sum(ad_vals)
                        alt_idx = alt_field.split(",").index(alt) + 1
                        if total > 0 and alt_idx < len(ad_vals):
                            vaf = ad_vals[alt_idx] / total
                    except (ValueError, IndexError):
                        vaf = None
                if vaf is None:
                    vaf = 1.0  # fallback: called as a confident site, treat as ~homozygous
                variants[(int(pos), ref, alt)] = vaf
    return variants


def per_sample_metrics(sample, gold, called):
    gold_keys = set(gold.keys())
    called_keys = set(called.keys())
    union = gold_keys | called_keys
    intersection = gold_keys & called_keys

    tp = len(intersection)
    fp = len(called_keys - gold_keys)
    fn = len(gold_keys - called_keys)
    concordance = len(intersection) / len(union) if union else 1.0
    precision = tp / (tp + fp) if (tp + fp) else 1.0
    recall = tp / (tp + fn) if (tp + fn) else 1.0

    pairs = []
    for key in union:
        true_vaf = gold.get(key, 0.0)
        called_vaf = called.get(key, 0.0)
        pairs.append((true_vaf, called_vaf))

    return {
        "sample": sample,
        "n_gold": len(gold_keys),
        "n_called": len(called_keys),
        "true_positive": tp,
        "false_positive": fp,
        "false_negative": fn,
        "concordance_rate": concordance,
        "precision": precision,
        "recall": recall,
    }, pairs


def main():
    gold_files = sorted(glob.glob(os.path.join(GOLD_DIR, "*.gold.vcf")))
    if not gold_files:
        print("No gold-standard VCFs found; run data/generate_synthetic_variant_data.py first.", file=sys.stderr)
        sys.exit(1)

    all_rows = []
    all_pairs = []

    for gold_path in gold_files:
        sample = os.path.basename(gold_path).replace(".gold.vcf", "")
        called_path = os.path.join(CALLED_DIR, f"{sample}.called.vcf")
        if not os.path.exists(called_path):
            print(f"WARNING: no called VCF for {sample}, skipping", file=sys.stderr)
            continue
        gold = parse_gold_vcf(gold_path)
        called = parse_called_vcf(called_path)
        row, pairs = per_sample_metrics(sample, gold, called)
        all_rows.append(row)
        all_pairs.extend(pairs)

    true_vafs = np.array([p[0] for p in all_pairs])
    called_vafs = np.array([p[1] for p in all_pairs])
    pearson_r, pearson_p = stats.pearsonr(true_vafs, called_vafs)

    mean_concordance = float(np.mean([r["concordance_rate"] for r in all_rows]))
    mean_precision = float(np.mean([r["precision"] for r in all_rows]))
    mean_recall = float(np.mean([r["recall"] for r in all_rows]))

    out_dir = os.path.join(HERE, "output")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "per_sample_concordance.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    summary = {
        "n_samples": len(all_rows),
        "n_variant_site_pairs_pooled": len(all_pairs),
        "pearson_correlation_vaf_vs_goldstandard": round(float(pearson_r), 4),
        "pearson_p_value": pearson_p,
        "mean_per_sample_concordance_rate": round(mean_concordance, 4),
        "mean_precision": round(mean_precision, 4),
        "mean_recall": round(mean_recall, 4),
    }

    print("== Concordance vs. gold-standard calls (25 samples) ==")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    summary_path = os.path.join(out_dir, "concordance_summary.csv")
    with open(summary_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        for k, v in summary.items():
            writer.writerow([k, v])

    print(f"\nWrote {csv_path}")
    print(f"Wrote {summary_path}")
    return summary


if __name__ == "__main__":
    main()
