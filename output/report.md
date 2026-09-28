# Variant Analysis Framework -- Run Report

Real pipeline run: `samtools`/`bcftools` variant calling on 25 synthetic BAM samples, compared against per-sample gold-standard VCFs.

## Concordance vs. gold-standard

- **n_samples**: 25
- **n_variant_site_pairs_pooled**: 521
- **pearson_correlation_vaf_vs_goldstandard**: 0.9135
- **pearson_p_value**: 8.214133056998314e-205
- **mean_per_sample_concordance_rate**: 0.7692
- **mean_precision**: 0.7706
- **mean_recall**: 0.9969

See `output/figures/per_sample_concordance.png` and `output/figures/precision_recall.png`.

## Classification report

# Variant Classification Report

Automated classification of called variants across **25 samples** using the synthetic gene/codon model in `variant_classification/gene_model.py`.

Generated automatically in **0.004 seconds** (measured wall-clock time for this run, see README for the honest 'report generation speed' framing).

- Total classified variants: 520
- Likely Pathogenic: 81
- VUS: 227
- Likely Benign: 212

## Per-sample summary

| Sample | Variants | Likely Pathogenic | VUS | Likely Benign |
|---|---|---|---|---|
| sample_01 | 19 | 3 | 8 | 8 |
| sample_02 | 24 | 8 | 10 | 6 |
| sample_03 | 20 | 4 | 10 | 6 |
| sample_04 | 14 | 3 | 6 | 5 |
| sample_05 | 23 | 3 | 9 | 11 |
| sample_06 | 25 | 3 | 9 | 13 |
| sample_07 | 18 | 5 | 8 | 5 |
| sample_08 | 22 | 1 | 11 | 10 |
| sample_09 | 21 | 4 | 10 | 7 |
| sample_10 | 22 | 6 | 5 | 11 |
| sample_11 | 20 | 1 | 12 | 7 |
| sample_12 | 17 | 1 | 10 | 6 |
| sample_13 | 18 | 3 | 4 | 11 |
| sample_14 | 26 | 4 | 15 | 7 |
| sample_15 | 19 | 3 | 7 | 9 |
| sample_16 | 23 | 0 | 10 | 13 |
| sample_17 | 19 | 5 | 7 | 7 |
| sample_18 | 18 | 1 | 8 | 9 |
| sample_19 | 23 | 4 | 10 | 9 |
| sample_20 | 19 | 2 | 9 | 8 |
| sample_21 | 22 | 5 | 11 | 6 |
| sample_22 | 24 | 3 | 12 | 9 |
| sample_23 | 25 | 3 | 10 | 12 |
| sample_24 | 18 | 5 | 5 | 8 |
| sample_25 | 21 | 1 | 11 | 9 |
