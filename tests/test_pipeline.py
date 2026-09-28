"""
Real pytest tests for the Variant Analysis Framework.

These tests assume `python3 run_all.py` (or at minimum
data/generate_synthetic_variant_data.py + run_variant_pipeline.sh +
compute_concordance.py) has already been run once, since they check the
actual on-disk outputs of the real samtools/bcftools pipeline rather than
mocking anything.
"""
import csv
import glob
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import compute_concordance as cc  # noqa: E402

N_SAMPLES = 25


def test_25_datasets_present():
    gold_files = glob.glob(os.path.join(ROOT, "data", "raw", "gold_vcf", "*.gold.vcf"))
    bam_files = glob.glob(os.path.join(ROOT, "data", "raw", "bam", "*.sorted.bam"))
    assert len(gold_files) == N_SAMPLES, f"expected {N_SAMPLES} gold VCFs, found {len(gold_files)}"
    assert len(bam_files) == N_SAMPLES, f"expected {N_SAMPLES} BAMs, found {len(bam_files)}"


def test_pipeline_produced_vcf_for_every_sample():
    gold_files = sorted(glob.glob(os.path.join(ROOT, "data", "raw", "gold_vcf", "*.gold.vcf")))
    assert gold_files, "no gold-standard VCFs found -- run data/generate_synthetic_variant_data.py first"
    for gold_path in gold_files:
        sample = os.path.basename(gold_path).replace(".gold.vcf", "")
        called_path = os.path.join(ROOT, "output", "called_vcf", f"{sample}.called.vcf")
        assert os.path.exists(called_path), f"missing called VCF for {sample} -- run run_variant_pipeline.sh"
        with open(called_path) as fh:
            content = fh.read()
        assert content.startswith("##fileformat=VCF"), f"{called_path} is not a valid VCF"


def test_concordance_summary_in_reasonable_band():
    summary_path = os.path.join(ROOT, "output", "concordance_summary.csv")
    assert os.path.exists(summary_path), "run compute_concordance.py first"
    values = {}
    with open(summary_path) as fh:
        for row in csv.reader(fh):
            values[row[0]] = row[1]
    pearson_r = float(values["pearson_correlation_vaf_vs_goldstandard"])
    assert 0.80 <= pearson_r <= 0.98, (
        f"pearson correlation {pearson_r} outside the expected 0.80-0.98 band "
        "for this calibrated synthetic dataset"
    )


def test_concordance_metric_sanity_zero_noise_is_near_perfect():
    """Sanity-check the concordance metric ITSELF (not the pipeline): feeding
    it a gold-standard set and an identical "called" set (i.e. the zero-noise
    case) must yield concordance/precision/recall of 1.0. This verifies
    per_sample_metrics() is implemented correctly, independent of whatever
    noise the synthetic BAMs happen to contain.
    """
    gold = {
        (100, "A", "T"): 1.0,
        (200, "C", "G"): 0.5,
        (300, "GA", "G"): 0.9,
    }
    called = dict(gold)  # identical -> zero discordance
    row, pairs = cc.per_sample_metrics("synthetic_zero_noise", gold, called)
    assert row["concordance_rate"] == pytest.approx(1.0)
    assert row["precision"] == pytest.approx(1.0)
    assert row["recall"] == pytest.approx(1.0)
    assert row["false_positive"] == 0
    assert row["false_negative"] == 0
    # every true/called VAF pair should also match exactly
    for true_vaf, called_vaf in pairs:
        assert true_vaf == pytest.approx(called_vaf)


def test_classification_report_covers_all_samples():
    report_path = os.path.join(ROOT, "output", "classification_report.csv")
    assert os.path.exists(report_path), "run variant_classification/classify_variants.py first"
    with open(report_path) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == N_SAMPLES
    for row in rows:
        assert int(row["n_variants"]) >= 0
        tier_sum = int(row["n_likely_pathogenic"]) + int(row["n_vus"]) + int(row["n_likely_benign"])
        assert tier_sum == int(row["n_variants"])
