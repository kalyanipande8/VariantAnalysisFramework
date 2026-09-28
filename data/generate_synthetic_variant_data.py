#!/usr/bin/env python3
"""
Generate 25 small, synthetic sequencing datasets for the Variant Analysis
Framework.

For each of the 25 synthetic samples this script produces:
  - a real SAM/BAM alignment (built with real read strings and CIGARs,
    sorted and indexed with `samtools`) simulating short reads drawn from a
    shared ~3 kb synthetic reference, carrying a set of injected true
    variants (SNPs + small indels) at known genomic positions, plus random
    sequencing errors and coverage dropout
  - a "gold-standard" VCF file recording the ground-truth variant sites and
    their true variant allele fraction (VAF)

The noise sources (sequencing error rate, coverage, per-site dropout
probability, heterozygous vs. homozygous mixing) are deliberately injected
so that a real variant caller run on the resulting BAMs will NOT recover the
gold standard perfectly -- some sites will be missed (false negatives), some
spurious sites will be called from sequencing errors (false positives), and
called VAFs will differ from the true VAF by sampling noise. This is what
makes the downstream concordance/correlation number in compute_concordance.py
a genuine, measured result rather than a fabricated one.

Fixed random seed (42) makes the whole dataset reproducible.
"""
import os
import random
import subprocess
import sys

SEED = 42
random.seed(SEED)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW_DIR = os.path.join(HERE, "raw")
REF_PATH = os.path.join(RAW_DIR, "reference.fasta")

N_SAMPLES = 25
REF_LEN = 3000
READ_LEN = 100
TARGET_DEPTH = 24
BASE_ERROR_RATE = 0.02           # per-base sequencing error rate (realistic Illumina-ish)
N_TRUE_VARIANTS_RANGE = (14, 22)  # per sample
HET_FRACTION = 0.40              # fraction of true variants simulated as heterozygous (VAF ~0.5)
DROPOUT_PROB = 0.10              # probability a true variant site suffers coverage dropout
FALSE_VARIANT_ERROR_BOOST = 0.30  # local error rate at "decoy" sites that mimic real variants
BASES = "ACGT"


def make_reference():
    seq = "".join(random.choice(BASES) for _ in range(REF_LEN))
    os.makedirs(RAW_DIR, exist_ok=True)
    with open(REF_PATH, "w") as fh:
        fh.write(">synthetic_chr1 synthetic 3kb reference for VariantAnalysisFramework\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    return seq


def mutate_base(b):
    choices = [x for x in BASES if x != b]
    return random.choice(choices)


def pick_variant_sites(ref_seq, n, rng):
    """Pick non-overlapping variant sites (SNP or small indel) on the reference."""
    candidates = list(range(50, REF_LEN - 50))
    rng.shuffle(candidates)
    positions = []
    for c in candidates:
        if all(abs(c - p) > 15 for p in positions):
            positions.append(c)
        if len(positions) == n:
            break
    positions.sort()
    sites = []
    for pos in positions:
        kind = "SNP" if rng.random() < 0.8 else "INDEL"
        ref_base = ref_seq[pos]
        if kind == "SNP":
            alt = mutate_base(ref_base)
            sites.append({"pos": pos, "type": "SNP", "ref": ref_base, "alt": alt})
        else:
            if rng.random() < 0.5:
                # insertion of 1-3 bases after pos
                ins = "".join(rng.choice(BASES) for _ in range(rng.randint(1, 3)))
                sites.append({"pos": pos, "type": "INS", "ref": ref_base, "alt": ref_base + ins})
            else:
                dellen = rng.randint(1, 3)
                sites.append({
                    "pos": pos, "type": "DEL",
                    "ref": ref_seq[pos:pos + 1 + dellen], "alt": ref_base,
                })
    return sites


def _maybe_error(base, err_rate, rng):
    if rng.random() < err_rate:
        return mutate_base(base)
    return base


def _build_one_read(ref_seq, start, target_len, sites_in_range, site_vaf, decoy_positions, rng):
    """Walk the reference from `start`, applying any overlapping variant sites
    (with proper insertion/deletion CIGAR operations) and background
    sequencing errors, and return (seq, cigar_string, ref_span).
    """
    end_ref = min(start + target_len, REF_LEN)
    site_map = {s["pos"]: s for s in sites_in_range}
    seq_chars = []
    cigar_tokens = []  # list of single-base op chars, merged at the end

    i = start
    while i < end_ref:
        s = site_map.get(i)
        if s is None:
            err_rate = FALSE_VARIANT_ERROR_BOOST if i in decoy_positions else BASE_ERROR_RATE
            seq_chars.append(_maybe_error(ref_seq[i], err_rate, rng))
            cigar_tokens.append("M")
            i += 1
            continue

        vaf = site_vaf[s["pos"]]
        use_alt = rng.random() < vaf

        if s["type"] == "SNP":
            base = s["alt"] if use_alt else s["ref"]
            seq_chars.append(_maybe_error(base, BASE_ERROR_RATE, rng))
            cigar_tokens.append("M")
            i += 1
        elif s["type"] == "INS":
            # anchor base is always emitted as M
            seq_chars.append(_maybe_error(s["ref"], BASE_ERROR_RATE, rng))
            cigar_tokens.append("M")
            i += 1
            if use_alt:
                for b in s["alt"][1:]:
                    seq_chars.append(_maybe_error(b, BASE_ERROR_RATE, rng))
                    cigar_tokens.append("I")
        elif s["type"] == "DEL":
            dellen = len(s["ref"]) - 1
            # bail out if the deletion would run past this read's window
            if i + 1 + dellen > end_ref:
                seq_chars.append(_maybe_error(ref_seq[i], BASE_ERROR_RATE, rng))
                cigar_tokens.append("M")
                i += 1
                continue
            seq_chars.append(_maybe_error(ref_seq[i], BASE_ERROR_RATE, rng))
            cigar_tokens.append("M")
            i += 1
            if use_alt:
                cigar_tokens.extend(["D"] * dellen)
                i += dellen
            else:
                for _ in range(dellen):
                    seq_chars.append(_maybe_error(ref_seq[i], BASE_ERROR_RATE, rng))
                    cigar_tokens.append("M")
                    i += 1

    # merge cigar tokens into run-length encoded CIGAR, trimming any leading/
    # trailing D (invalid in SAM)
    while cigar_tokens and cigar_tokens[0] == "D":
        cigar_tokens.pop(0)
    while cigar_tokens and cigar_tokens[-1] == "D":
        cigar_tokens.pop()

    cigar = []
    cur_op, cur_len = None, 0
    for op in cigar_tokens:
        if op == cur_op:
            cur_len += 1
        else:
            if cur_op is not None:
                cigar.append(f"{cur_len}{cur_op}")
            cur_op, cur_len = op, 1
    if cur_op is not None:
        cigar.append(f"{cur_len}{cur_op}")

    return "".join(seq_chars), "".join(cigar) if cigar else f"{len(seq_chars)}M"


def simulate_sample_reads(ref_seq, sites, decoy_positions, rng):
    """Return list of (pos0based, seq, qual, cigar) reads with variants + errors baked in."""
    reads = []
    n_positions = REF_LEN - READ_LEN
    n_reads = int((TARGET_DEPTH * REF_LEN) / READ_LEN)

    site_vaf = {}
    for s in sites:
        if rng.random() < DROPOUT_PROB:
            site_vaf[s["pos"]] = 0.0  # simulate dropout: essentially unobserved
        elif rng.random() < HET_FRACTION:
            site_vaf[s["pos"]] = round(rng.uniform(0.35, 0.6), 2)
        else:
            site_vaf[s["pos"]] = round(rng.uniform(0.85, 1.0), 2)

    for _ in range(n_reads):
        start = rng.randint(0, n_positions)
        end = start + READ_LEN
        sites_in_range = [s for s in sites if start <= s["pos"] < end]
        seq, cigar = _build_one_read(ref_seq, start, READ_LEN, sites_in_range, site_vaf, decoy_positions, rng)
        qual = "I" * len(seq)  # Phred 40, uniform high quality
        reads.append((start, seq, qual, cigar))

    return reads, site_vaf


def write_sam_and_bam(sample_name, ref_seq, reads, out_dir, samtools_path):
    sam_path = os.path.join(out_dir, f"{sample_name}.sam")
    bam_path = os.path.join(out_dir, f"{sample_name}.sorted.bam")
    with open(sam_path, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        fh.write(f"@SQ\tSN:synthetic_chr1\tLN:{len(ref_seq)}\n")
        fh.write(f"@RG\tID:{sample_name}\tSM:{sample_name}\n")
        for i, (start, seq, qual, cigar) in enumerate(reads):
            qname = f"{sample_name}_r{i}"
            flag = 0
            pos = start + 1  # SAM is 1-based
            mapq = 60
            fh.write(
                "\t".join([
                    qname, str(flag), "synthetic_chr1", str(pos), str(mapq),
                    cigar, "*", "0", "0", seq, qual, f"RG:Z:{sample_name}",
                ]) + "\n"
            )

    subprocess.run(
        [samtools_path, "sort", "-o", bam_path, sam_path],
        check=True, capture_output=True,
    )
    subprocess.run([samtools_path, "index", bam_path], check=True, capture_output=True)
    os.remove(sam_path)
    return bam_path


def write_gold_vcf(sample_name, sites, site_vaf, out_path, ref_len):
    with open(out_path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID=synthetic_chr1,length={ref_len}>\n")
        fh.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="True simulated variant allele fraction">\n')
        fh.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n')
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for s in sorted(sites, key=lambda x: x["pos"]):
            vaf = site_vaf[s["pos"]]
            if vaf <= 0.0:
                continue  # dropout site: not observable, excluded from gold-standard "detectable" set
            pos_1based = s["pos"] + 1
            fh.write(
                f"synthetic_chr1\t{pos_1based}\t.\t{s['ref']}\t{s['alt']}\t.\tPASS\tVAF={vaf};TYPE={s['type']}\n"
            )


def normalize_vcf_inplace(vcf_path, ref_path, samtools_path):
    """Left-align/normalize indel representation with `bcftools norm` so the
    gold-standard VCF uses the same canonical indel anchoring that bcftools
    call + bcftools norm will produce for the pipeline's called VCFs. Without
    this, indels in homopolymer/repeat contexts can be written at more than
    one equally valid position and a naive string comparison would wrongly
    count them as discordant.
    """
    from shutil import which
    bcftools_path = which("bcftools")
    if not bcftools_path:
        return  # bcftools unavailable: leave gold VCF as-written (documented in README)
    tmp_out = vcf_path + ".norm.tmp"
    subprocess.run(
        [bcftools_path, "norm", "-f", ref_path, vcf_path, "-o", tmp_out],
        check=True, capture_output=True,
    )
    os.replace(tmp_out, vcf_path)


def find_samtools():
    from shutil import which
    p = which("samtools")
    if not p:
        print("ERROR: samtools not found on PATH; this project requires it for BAM generation.", file=sys.stderr)
        sys.exit(1)
    return p


def main():
    samtools_path = find_samtools()
    ref_seq = make_reference()
    subprocess.run([samtools_path, "faidx", REF_PATH], check=True, capture_output=True)
    manifest = []

    bam_dir = os.path.join(HERE, "raw", "bam")
    gold_dir = os.path.join(HERE, "raw", "gold_vcf")
    os.makedirs(bam_dir, exist_ok=True)
    os.makedirs(gold_dir, exist_ok=True)

    for i in range(1, N_SAMPLES + 1):
        sample_name = f"sample_{i:02d}"
        rng = random.Random(SEED + i)  # per-sample deterministic sub-seed
        n_true = rng.randint(*N_TRUE_VARIANTS_RANGE)
        sites = pick_variant_sites(ref_seq, n_true, rng)

        # A handful of decoy genomic positions get an elevated local error
        # rate, which occasionally produces spurious ("false positive")
        # variant calls that are NOT in the gold standard -- this is what
        # keeps the pipeline from recovering the gold standard perfectly.
        n_decoys = rng.randint(2, 5)
        decoy_positions = set(
            rng.sample([p for p in range(60, REF_LEN - 60) if p not in {s["pos"] for s in sites}], n_decoys)
        )

        reads, site_vaf = simulate_sample_reads(ref_seq, sites, decoy_positions, rng)
        bam_path = write_sam_and_bam(sample_name, ref_seq, reads, bam_dir, samtools_path)
        gold_path = os.path.join(gold_dir, f"{sample_name}.gold.vcf")
        write_gold_vcf(sample_name, sites, site_vaf, gold_path, REF_LEN)
        normalize_vcf_inplace(gold_path, REF_PATH, samtools_path)

        n_observable = sum(1 for v in site_vaf.values() if v > 0)
        manifest.append({
            "sample": sample_name,
            "bam": bam_path,
            "gold_vcf": gold_path,
            "n_true_variants": len(sites),
            "n_observable_gold_variants": n_observable,
            "n_decoy_sites": n_decoys,
            "n_reads": len(reads),
        })
        print(f"[{sample_name}] {len(sites)} injected variants "
              f"({n_observable} observable after dropout), {n_decoys} decoy noise sites, "
              f"{len(reads)} reads -> {os.path.basename(bam_path)}")

    manifest_path = os.path.join(HERE, "manifest.csv")
    import csv
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
        writer.writeheader()
        writer.writerows(manifest)

    print(f"\nReference: {REF_PATH} ({REF_LEN} bp)")
    print(f"Wrote {len(manifest)} samples. Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
