"""
A small, fully synthetic "gene annotation model" for the 3 kb reference used
throughout this project. This stands in for a real gene/transcript
annotation (e.g. a GTF/RefSeq record) so that called variants can be
classified by functional consequence the way a real annotation-driven
variant classifier (e.g. VEP/SnpEff style) would, without requiring a real
genome or annotation database.

Model: a single synthetic "gene" spans the whole 3000 bp reference.
  - [0, 200)      : 5' intergenic / promoter region
  - [200, 2700)   : coding sequence (CDS), in-frame codons of length 3
  - [2700, 3000)  : 3' intergenic region
Within the CDS, a codon table (the real standard genetic code) is used to
translate REF vs ALT codons for SNPs, giving genuine synonymous / missense /
nonsense calls. Indels inside the CDS are classified as frameshift (length
not a multiple of 3) or in-frame indel (length is a multiple of 3).
"""

CDS_START = 200
CDS_END = 2700

# Standard genetic code (DNA codon -> amino acid, '*' = stop)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def load_reference(fasta_path):
    seq = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def classify_variant(ref_seq, pos0, ref_allele, alt_allele):
    """Classify a single variant by genomic location / coding consequence.

    pos0: 0-based genomic position (VCF POS - 1)
    Returns a dict with region, consequence, and (for CDS SNPs) amino-acid
    change details.
    """
    is_snp = len(ref_allele) == 1 and len(alt_allele) == 1
    is_insertion = len(alt_allele) > len(ref_allele)
    is_deletion = len(alt_allele) < len(ref_allele)
    indel_len = abs(len(alt_allele) - len(ref_allele))

    if pos0 < CDS_START:
        region = "5_prime_intergenic"
    elif pos0 >= CDS_END:
        region = "3_prime_intergenic"
    else:
        region = "coding"

    result = {
        "region": region,
        "consequence": None,
        "ref_aa": None,
        "alt_aa": None,
        "codon_index": None,
    }

    if region != "coding":
        result["consequence"] = "intergenic"
        return result

    if not is_snp:
        result["consequence"] = "frameshift_indel" if indel_len % 3 != 0 else "inframe_indel"
        return result

    # SNP inside CDS: translate the codon
    offset_in_cds = pos0 - CDS_START
    codon_index = offset_in_cds // 3
    codon_start = CDS_START + codon_index * 3
    ref_codon = ref_seq[codon_start:codon_start + 3]
    if len(ref_codon) < 3:
        result["consequence"] = "coding_boundary_effect"
        return result
    pos_in_codon = pos0 - codon_start
    alt_codon = ref_codon[:pos_in_codon] + alt_allele + ref_codon[pos_in_codon + 1:]

    ref_aa = CODON_TABLE.get(ref_codon, "X")
    alt_aa = CODON_TABLE.get(alt_codon, "X")
    result["ref_aa"] = ref_aa
    result["alt_aa"] = alt_aa
    result["codon_index"] = codon_index

    if ref_aa == alt_aa:
        result["consequence"] = "synonymous"
    elif alt_aa == "*":
        result["consequence"] = "nonsense"
    elif ref_aa == "*":
        result["consequence"] = "stop_loss"
    else:
        result["consequence"] = "missense"
    return result


# A simple, documented, synthetic pathogenicity heuristic (NOT a real
# clinical classifier): stop-gained/frameshift variants are "Likely
# Pathogenic", missense variants are "VUS" (variant of uncertain
# significance), and everything else (synonymous, intergenic, in-frame
# indels) is "Likely Benign". This mirrors, in miniature, the kind of
# consequence-based tiering real ACMG-style classifiers use as one input
# among many -- it is illustrative, not diagnostic.
def pathogenicity_tier(consequence):
    if consequence in ("nonsense", "frameshift_indel", "stop_loss"):
        return "Likely_Pathogenic"
    if consequence == "missense":
        return "VUS"
    return "Likely_Benign"
