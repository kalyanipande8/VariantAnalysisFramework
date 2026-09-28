#!/usr/bin/env python3
"""
Python variant classification framework.

Reads each sample's real bcftools-called VCF (output/called_vcf/*.called.vcf)
and classifies every called variant using the synthetic gene/codon model in
gene_model.py (region + coding consequence + an illustrative pathogenicity
tier). Writes one per-sample-variant CSV and a per-sample summary report
(classification_report.csv / classification_report.md) across all 25
samples.

This is the automated stand-in for what the resume bullet calls "4-day
faster report generation": see the README for the honest framing of that
claim (there is no real historical multi-day manual baseline measured here --
what IS measured, and reported below, is the real wall-clock time this
automated classification + report step takes to run across all 25 samples).
"""
import csv
import glob
import os
import re
import time

from gene_model import CDS_START, CDS_END, classify_variant, load_reference, pathogenicity_tier

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REF_PATH = os.path.join(ROOT, "data", "raw", "reference.fasta")
CALLED_DIR = os.path.join(ROOT, "output", "called_vcf")
OUT_DIR = os.path.join(ROOT, "output")


def parse_called_vcf_records(path):
    records = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos, ref, alt_field, info = fields[1], fields[3], fields[4], fields[7]
            fmt_keys = fields[8].split(":")
            fmt_vals = fields[9].split(":")
            fmt = dict(zip(fmt_keys, fmt_vals))
            dp_match = re.search(r"(?:^|;)DP=(\d+)", info)
            dp = int(dp_match.group(1)) if dp_match else None
            for alt in alt_field.split(","):
                records.append({
                    "pos": int(pos), "ref": ref, "alt": alt,
                    "gt": fmt.get("GT", "."), "dp": dp,
                })
    return records


def classify_all_samples():
    t0 = time.time()
    ref_seq = load_reference(REF_PATH)
    vcf_files = sorted(glob.glob(os.path.join(CALLED_DIR, "*.called.vcf")))
    if not vcf_files:
        raise SystemExit("No called VCFs found -- run run_variant_pipeline.sh first.")

    all_variant_rows = []
    sample_summaries = []

    for vcf_path in vcf_files:
        sample = os.path.basename(vcf_path).replace(".called.vcf", "")
        records = parse_called_vcf_records(vcf_path)
        tier_counts = {"Likely_Pathogenic": 0, "VUS": 0, "Likely_Benign": 0}
        consequence_counts = {}

        for rec in records:
            pos0 = rec["pos"] - 1
            result = classify_variant(ref_seq, pos0, rec["ref"], rec["alt"])
            tier = pathogenicity_tier(result["consequence"])
            tier_counts[tier] += 1
            consequence_counts[result["consequence"]] = consequence_counts.get(result["consequence"], 0) + 1

            all_variant_rows.append({
                "sample": sample,
                "pos": rec["pos"],
                "ref": rec["ref"],
                "alt": rec["alt"],
                "genotype": rec["gt"],
                "depth": rec["dp"],
                "region": result["region"],
                "consequence": result["consequence"],
                "ref_aa": result["ref_aa"],
                "alt_aa": result["alt_aa"],
                "pathogenicity_tier": tier,
            })

        sample_summaries.append({
            "sample": sample,
            "n_variants": len(records),
            "n_likely_pathogenic": tier_counts["Likely_Pathogenic"],
            "n_vus": tier_counts["VUS"],
            "n_likely_benign": tier_counts["Likely_Benign"],
        })

    os.makedirs(OUT_DIR, exist_ok=True)
    variants_csv = os.path.join(OUT_DIR, "classified_variants.csv")
    with open(variants_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_variant_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_variant_rows)

    summary_csv = os.path.join(OUT_DIR, "classification_report.csv")
    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(sample_summaries[0].keys()))
        writer.writeheader()
        writer.writerows(sample_summaries)

    elapsed = time.time() - t0

    total_variants = len(all_variant_rows)
    total_path = sum(s["n_likely_pathogenic"] for s in sample_summaries)
    total_vus = sum(s["n_vus"] for s in sample_summaries)
    total_benign = sum(s["n_likely_benign"] for s in sample_summaries)

    md_path = os.path.join(OUT_DIR, "classification_report.md")
    with open(md_path, "w") as fh:
        fh.write("# Variant Classification Report\n\n")
        fh.write(f"Automated classification of called variants across **{len(sample_summaries)} samples** "
                 f"using the synthetic gene/codon model in `variant_classification/gene_model.py`.\n\n")
        fh.write(f"Generated automatically in **{elapsed:.3f} seconds** "
                 f"(measured wall-clock time for this run, see README for the honest 'report generation speed' framing).\n\n")
        fh.write(f"- Total classified variants: {total_variants}\n")
        fh.write(f"- Likely Pathogenic: {total_path}\n")
        fh.write(f"- VUS: {total_vus}\n")
        fh.write(f"- Likely Benign: {total_benign}\n\n")
        fh.write("## Per-sample summary\n\n")
        fh.write("| Sample | Variants | Likely Pathogenic | VUS | Likely Benign |\n")
        fh.write("|---|---|---|---|---|\n")
        for s in sample_summaries:
            fh.write(f"| {s['sample']} | {s['n_variants']} | {s['n_likely_pathogenic']} | "
                      f"{s['n_vus']} | {s['n_likely_benign']} |\n")

    print(f"Classified {total_variants} variants across {len(sample_summaries)} samples in {elapsed:.3f}s")
    print(f"Wrote {variants_csv}")
    print(f"Wrote {summary_csv}")
    print(f"Wrote {md_path}")
    return {
        "elapsed_seconds": elapsed,
        "total_variants": total_variants,
        "n_samples": len(sample_summaries),
    }


if __name__ == "__main__":
    classify_all_samples()
