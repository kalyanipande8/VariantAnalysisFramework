#!/usr/bin/env python3
"""
End-to-end orchestration: generate synthetic data -> run the real
samtools/bcftools variant-calling pipeline -> compute concordance vs.
gold-standard -> classify variants (Python) -> write output/report.md with
real, measured numbers and figures.
"""
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, **kw):
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=HERE, **kw)


def main():
    t_start = time.time()

    run([sys.executable, os.path.join("data", "generate_synthetic_variant_data.py")])
    run(["bash", os.path.join(HERE, "run_variant_pipeline.sh")])

    print("\n$ python3 compute_concordance.py")
    sys.path.insert(0, HERE)
    import compute_concordance
    concordance_summary = compute_concordance.main()

    run([sys.executable, os.path.join("variant_classification", "classify_variants.py")])

    r_available = subprocess.run(["which", "Rscript"], capture_output=True).returncode == 0
    if r_available:
        try:
            run(["Rscript", os.path.join("variant_classification", "classify_variants.R")])
        except subprocess.CalledProcessError:
            print("R classification script failed to run; continuing without it.")

    make_figures()
    write_report(concordance_summary)

    elapsed = time.time() - t_start
    print(f"\nDone. Total pipeline wall-clock time: {elapsed:.1f}s")
    print(f"See {os.path.join(HERE, 'output', 'report.md')}")


def make_figures():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figures.")
        return

    fig_dir = os.path.join(HERE, "output", "figures")
    os.makedirs(fig_dir, exist_ok=True)

    per_sample_path = os.path.join(HERE, "output", "per_sample_concordance.csv")
    samples, concordance, precision, recall = [], [], [], []
    with open(per_sample_path) as fh:
        for row in csv.DictReader(fh):
            samples.append(row["sample"].replace("sample_", "S"))
            concordance.append(float(row["concordance_rate"]))
            precision.append(float(row["precision"]))
            recall.append(float(row["recall"]))

    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.bar(samples, concordance, color="#3B7DD8")
    ax.set_ylabel("Concordance rate (intersection / union)")
    ax.set_title("Per-sample variant concordance vs. gold standard (25 synthetic samples)")
    ax.set_ylim(0, 1.05)
    plt.xticks(rotation=90, fontsize=7)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, "per_sample_concordance.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(precision, recall, color="#D8763B")
    for x, y, s in zip(precision, recall, samples):
        ax.annotate(s, (x, y), fontsize=6, alpha=0.7)
    ax.set_xlabel("Precision")
    ax.set_ylabel("Recall")
    ax.set_title("Precision vs. recall per sample")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, "precision_recall.png"), dpi=150)
    plt.close(fig)

    print(f"Wrote figures to {fig_dir}")


def write_report(summary):
    out_path = os.path.join(HERE, "output", "report.md")

    class_report_path = os.path.join(HERE, "output", "classification_report.md")
    class_report_text = ""
    if os.path.exists(class_report_path):
        with open(class_report_path) as fh:
            class_report_text = fh.read()

    with open(out_path, "w") as fh:
        fh.write("# Variant Analysis Framework -- Run Report\n\n")
        fh.write("Real pipeline run: `samtools`/`bcftools` variant calling on 25 synthetic "
                 "BAM samples, compared against per-sample gold-standard VCFs.\n\n")
        fh.write("## Concordance vs. gold-standard\n\n")
        for k, v in summary.items():
            fh.write(f"- **{k}**: {v}\n")
        fh.write("\nSee `output/figures/per_sample_concordance.png` and "
                  "`output/figures/precision_recall.png`.\n\n")
        fh.write("## Classification report\n\n")
        fh.write(class_report_text)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
