"""
Design two-fragment Golden Gate assemblies from DNA sequences.

Takes sequences in the 250-500bp range and fragments them into two oligos
bridged by an optimal shared Golden Gate site. Fidelity prediction is adapted
from OMEGA (Freschlin et al., 2025 -- https://github.com/RomeroLab/omega).

Oligo structure:
  Fragment A: fwd_primer - GGTCTCA - 5'_bbsite - seq[0:split+4] - AGAGACC - rc(rev_primer)
  Fragment B: fwd_primer - GGTCTCA - seq[split:] - 3'_bbsite  - AGAGACC - rc(rev_primer)
"""

import os
import sys
import argparse
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# BsaI info
ENZYME_PREFIX   = 'GGTCTCA' # BsaI site 'GGTCTC' + 'A' padding base
ENZYME_SUFFIX   = 'AGAGACC' # Reverse complement of prefix
ENZYME_OVERHEAD = len(ENZYME_PREFIX) + len(ENZYME_SUFFIX)  # 14bp per oligo
OVERHANG_SIZE   = 4

# Seq info
SEQ_MIN_LEN  = 300
SEQ_MAX_LEN  = 500
OLIGO_MAX_LEN = 300

# Remove palindromic sites (per OMEGA Oligos)
PALINDROMIC_SITES = {
    'GCGC', 'CGCG', 'ATAT', 'TATA', 'GGCC', 'CCGG',
    'AATT', 'TTAA', 'TGCA', 'AGCT', 'TCGA', 'ACGT',
    'GATC', 'GTAC', 'CATG', 'CTAG'
}

# File locations
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_LIGATION_DATA = os.path.join(SCRIPT_DIR, '..', 'ligation_data', 'T4_37C.csv')
DEFAULT_PRIMERS = os.path.join(SCRIPT_DIR, '..', 'primers', 'primers.csv')

# Fidelity prediction adapted from OMEGA (Freschlin et al., 2025)

def load_ligation_data(filepath: str) -> pd.DataFrame:
    """Load Potapov ligation frequency matrix from CSV."""
    return pd.read_csv(filepath, index_col=0)


def correct_ligations(site: str, wc_site: str, data: pd.DataFrame) -> float:
    """Correct ligation count for a site with its Watson-Crick pair."""
    return data.loc[site, wc_site]


def total_ligations(site: str, sites: list, data: pd.DataFrame) -> float:
    """Total ligation events observed for a site against all sites in the set."""
    wc_sites = [str(Seq(s).reverse_complement()) for s in sites]
    cols = [c for c in list(sites) + wc_sites if c in data.columns]
    return data.loc[site, cols].sum()


def site_probability(site: str, all_sites: list, data: pd.DataFrame) -> float:
    """Orthogonality probability for a single GG site in context of a full site set."""
    wc_site = str(Seq(site).reverse_complement())

    correct = (
        correct_ligations(site, wc_site, data) +
        correct_ligations(wc_site, site, data)
    )
    total = (
        total_ligations(site, all_sites, data) +
        total_ligations(wc_site, all_sites, data)
    )

    if total == 0:
        return 0.0
    return correct / total


def predict_fidelity(sites: list, data: pd.DataFrame) -> float:
    """Predicted assembly fidelity for a set of GG sites (product of per-site probabilities)."""
    return float(np.prod([site_probability(s, sites, data) for s in sites]))

# Input parsing
def parse_fasta(filepath: str) -> list:
    sequences = []
    for record in SeqIO.parse(filepath, 'fasta'):
        sequences.append((record.id, str(record.seq).upper()))
    return sequences


def parse_csv(filepath: str, id_col: str, seq_col: str) -> list:
    df = pd.read_csv(filepath)
    for col in [id_col, seq_col]:
        if col not in df.columns:
            print(f"[ERROR] Column '{col}' not found in CSV.")
            print(f"  Available columns: {list(df.columns)}")
            print("  Use --csv_id_column and --csv_sequence_column to specify correct names.")
            sys.exit(1)
    return list(zip(df[id_col].astype(str), df[seq_col].str.upper()))


def load_primers(filepath: str) -> list:
    """Load primers from CSV with 'name' and 'sequence' columns."""
    df = pd.read_csv(filepath)
    for col in ['name', 'sequence']:
        if col not in df.columns:
            print(f"[ERROR] Column '{col}' not found in primer CSV.")
            print(f"  Available columns: {list(df.columns)}")
            sys.exit(1)
    return list(zip(df['name'], df['sequence'].str.upper()))


# Seq filtering
def filter_sequences(sequences: list) -> list:
    """Discard sequences outside the 250-500bp range."""
    kept, skipped = [], []
    for seq_id, seq in sequences:
        if SEQ_MIN_LEN <= len(seq) <= SEQ_MAX_LEN:
            kept.append((seq_id, seq))
        else:
            skipped.append((seq_id, len(seq)))

    for seq_id, length in skipped:
        print(f"  [SKIPPED] {seq_id}: {length} bp -- outside {SEQ_MIN_LEN}-{SEQ_MAX_LEN} bp range.")

    return kept


# Backbone gg validation
def validate_backbone_sites(upstream: str, downstream: str) -> None:
    """Warn if backbone sites may cause assembly issues."""
    downstream_rc = str(Seq(downstream).reverse_complement())

    if upstream in PALINDROMIC_SITES:
        warnings.warn(
            f"Upstream site '{upstream}' is palindromic and may reduce fidelity.",
            stacklevel=2
        )
    if downstream in PALINDROMIC_SITES:
        warnings.warn(
            f"Downstream site '{downstream}' is palindromic and may reduce fidelity.",
            stacklevel=2
        )
    if upstream == downstream:
        warnings.warn(
            f"Upstream and downstream sites are identical ('{upstream}'). "
            "This will prevent directional cloning.",
            stacklevel=2
        )
    if upstream == downstream_rc:
        warnings.warn(
            f"Upstream site is the reverse complement of the downstream site. "
            "This may cause mis-assembly.",
            stacklevel=2
        )


# gg finder 
def get_disallowed_sites(upstream: str, downstream: str) -> set:
    """Build set of overhangs that cannot be used as the shared internal site."""
    disallowed = set(PALINDROMIC_SITES)
    for site in [upstream, downstream]:
        disallowed.add(site)
        disallowed.add(str(Seq(site).reverse_complement()))
    return disallowed


def find_best_split(
    seq: str,
    upstream: str,
    downstream: str,
    fwd_primer: str,
    rev_primer: str,
    ligation_data: pd.DataFrame,
    disallowed_sites: set
) -> tuple | None:
    """Find the optimal position for the shared internal GG site.

    Enumerates all valid positions, scores each using Potapov ligation data
    (3-site system: upstream, shared, downstream), and returns the highest-fidelity
    position. Ties are broken by proximity to the sequence centre.

    Returns (split_position, shared_ggsite, fidelity) or None if no valid position.
    """
    n = len(seq)
    center = (n - OVERHANG_SIZE) / 2.0

    # Oligo length constraints:
    # Oligo A: fwd(L) + prefix(7) + upstream(4) + seq[0:p+4] + suffix(7) + rev(L)
    #   total = len(fwd) + len(rev) + p + 22  ≤  OLIGO_MAX_LEN
    #   → p ≤ OLIGO_MAX_LEN - len(fwd) - len(rev) - 22
    #
    # Oligo B: fwd(L) + prefix(7) + seq[p:] + downstream(4) + suffix(7) + rev(L)
    #   total = len(fwd) + len(rev) + (n - p) + 18  ≤  OLIGO_MAX_LEN
    #   → p ≥ len(fwd) + len(rev) + n - OLIGO_MAX_LEN + 18

    primer_len = len(fwd_primer) + len(rev_primer)
    p_max = OLIGO_MAX_LEN - primer_len - 22
    p_min = primer_len + n - OLIGO_MAX_LEN + 18

    p_min = max(0, p_min)
    p_max = min(n - OVERHANG_SIZE, p_max)

    if p_min > p_max:
        return None

    scored = []
    for p in range(p_min, p_max + 1):
        ggsite = seq[p:p + OVERHANG_SIZE]
        if ggsite in disallowed_sites:
            continue
        fidelity = predict_fidelity([upstream, ggsite, downstream], ligation_data)
        dist_from_center = abs(p - center)
        scored.append((p, ggsite, fidelity, dist_from_center))

    if not scored:
        return None

    # Best fidelity first w/ ties broken by closest to center
    scored.sort(key=lambda x: (-x[2], x[3]))
    p, ggsite, fidelity, _ = scored[0]
    return p, ggsite, fidelity


# Oligo builder
def build_oligos(
    seq: str,
    upstream: str,
    downstream: str,
    split_pos: int,
    fwd_primer: str,
    rev_primer: str
) -> tuple:
    """Build the two orderable oligo sequences.

    Fragment A: fwd - GGTCTCA - upstream_bbsite - seq[0:split+4] - AGAGACC - rc(rev)
    Fragment B: fwd - GGTCTCA - seq[split:] - downstream_bbsite - AGAGACC - rc(rev)

    The shared 4bp overhang (seq[split:split+4]) is present in both fragments
    and becomes the scarless ligation junction.
    """
    p = split_pos
    rev_rc = str(Seq(rev_primer).reverse_complement())

    coding_a = upstream + seq[:p + OVERHANG_SIZE]
    coding_b = seq[p:] + downstream

    oligo_a = fwd_primer + ENZYME_PREFIX + coding_a + ENZYME_SUFFIX + rev_rc
    oligo_b = fwd_primer + ENZYME_PREFIX + coding_b + ENZYME_SUFFIX + rev_rc

    return oligo_a, oligo_b


# Primer pair 
def assign_primer_pairs(sequences: list, primers: list) -> list:
    """Assign a unique primer pair to each sequence.

    Each gene gets a (fwd, rev) pair where fwd ≠ rev. Pairs are unordered:
    (p1, p2) and (p2, p1) are treated as the same pair and will not both be used.

    Returns list of (fwd_name, fwd_seq, rev_name, rev_seq) per sequence.
    """
    all_pairs = list(combinations(range(len(primers)), 2))
    n_seqs = len(sequences)

    if n_seqs > len(all_pairs):
        print(
            f"[ERROR] Not enough unique primer pairs. "
            f"{len(primers)} primers yield {len(all_pairs)} unique pairs, "
            f"but {n_seqs} sequences need assignment."
        )
        sys.exit(1)

    assignments = []
    for i in range(n_seqs):
        idx1, idx2 = all_pairs[i]
        name1, seq1 = primers[idx1]
        name2, seq2 = primers[idx2]
        assignments.append((name1, seq1, name2, seq2))

    return assignments


# Output
def write_outputs(results: list, output_dir: str) -> None:
    """Write order.csv and fidelity_report.csv."""
    os.makedirs(output_dir, exist_ok=True)

    # Order
    order_rows = []
    for r in results:
        order_rows.append({'name': f"{r['id']}_A", 'sequence': r['oligo_a']})
        order_rows.append({'name': f"{r['id']}_B", 'sequence': r['oligo_b']})

    order_path = os.path.join(output_dir, 'order.csv')
    pd.DataFrame(order_rows).to_csv(order_path, index=False)
    print(f"  Orderable oligos written to:  {order_path}")

    # Report
    fidelity_rows = [{
        'id':                r['id'],
        'sequence_length':   r['sequence_length'],
        'split_position':    r['split_position'],
        'shared_ggsite':     r['shared_ggsite'],
        'upstream_bbsite':   r['upstream_bbsite'],
        'downstream_bbsite': r['downstream_bbsite'],
        'fidelity':          round(r['fidelity'], 6),
        'oligo_a_length':    len(r['oligo_a']),
        'oligo_b_length':    len(r['oligo_b']),
        'fwd_primer':        r['fwd_primer_name'],
        'rev_primer':        r['rev_primer_name'],
    } for r in results]

    fidelity_path = os.path.join(output_dir, 'fidelity_report.csv')
    pd.DataFrame(fidelity_rows).to_csv(fidelity_path, index=False)
    print(f"  Fidelity report written to:   {fidelity_path}")


# Main
def main(args: argparse.Namespace) -> None:

    # Load ligation data
    print(f"\n[Setup] Loading ligation data: {args.ligation_data}")
    ligation_data = load_ligation_data(args.ligation_data)

    # Validate and prepare backbone sites
    upstream  = args.upstream_site.upper()
    downstream = args.downstream_site.upper()
    print(f"[Setup] Backbone sites: 5' = {upstream}, 3' = {downstream}")
    validate_backbone_sites(upstream, downstream)
    disallowed_sites = get_disallowed_sites(upstream, downstream)

    # Load sequences
    print(f"\n[Input] Loading sequences: {args.input}")
    ext = os.path.splitext(args.input)[1].lower()
    if ext in ('.fasta', '.fa'):
        raw_sequences = parse_fasta(args.input)
    elif ext == '.csv':
        raw_sequences = parse_csv(args.input, args.csv_id_column, args.csv_sequence_column)
    else:
        print(f"[ERROR] Unrecognized file extension '{ext}'. Use .fasta, .fa, or .csv.")
        sys.exit(1)
    print(f"  Loaded {len(raw_sequences)} sequences.")

    # Filter by length
    print(f"[Filter] Keeping sequences between {SEQ_MIN_LEN}-{SEQ_MAX_LEN} bp...")
    sequences = filter_sequences(raw_sequences)
    print(f"  {len(sequences)} sequences passed.")

    if not sequences:
        print("[ERROR] No sequences passed the length filter. Exiting.")
        sys.exit(1)

    # Load and assign primers
    print(f"\n[Primers] Loading: {args.primers}")
    primers = load_primers(args.primers)
    n_pairs = len(primers) * (len(primers) - 1) // 2
    print(f"  Loaded {len(primers)} primers ({n_pairs} unique pairs available).")

    print("[Primers] Assigning unique pairs to sequences...")
    primer_assignments = assign_primer_pairs(sequences, primers)

    # Design fragments
    print("\n[Design] Finding optimal shared GG sites...")
    results, failed = [], []

    for (seq_id, seq), (fwd_name, fwd_seq, rev_name, rev_seq) in zip(sequences, primer_assignments):

        result = find_best_split(
            seq=seq,
            upstream=upstream,
            downstream=downstream,
            fwd_primer=fwd_seq,
            rev_primer=rev_seq,
            ligation_data=ligation_data,
            disallowed_sites=disallowed_sites
        )

        if result is None:
            print(f"  [FAILED] {seq_id}: no valid GG site found within oligo length constraints.")
            failed.append(seq_id)
            continue

        split_pos, shared_ggsite, fidelity = result

        oligo_a, oligo_b = build_oligos(
            seq=seq,
            upstream=upstream,
            downstream=downstream,
            split_pos=split_pos,
            fwd_primer=fwd_seq,
            rev_primer=rev_seq
        )

        print(f"  {seq_id}: split @ {split_pos}, shared site = {shared_ggsite}, "
              f"fidelity = {fidelity:.4f}, oligos = {len(oligo_a)}/{len(oligo_b)} bp")

        results.append({
            'id':                seq_id,
            'sequence_length':   len(seq),
            'split_position':    split_pos,
            'shared_ggsite':     shared_ggsite,
            'upstream_bbsite':   upstream,
            'downstream_bbsite': downstream,
            'fidelity':          fidelity,
            'oligo_a':           oligo_a,
            'oligo_b':           oligo_b,
            'fwd_primer_name':   fwd_name,
            'rev_primer_name':   rev_name,
        })

    # Write outputs
    if results:
        print(f"\n[Output] Writing to: {args.output_dir}")
        write_outputs(results, args.output_dir)

    print(f"\n[Done] {len(results)}/{len(sequences)} sequences designed successfully.")
    if failed:
        print(f"  Failed ({len(failed)}): {failed}")


# Entry point
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "Design two-fragment Golden Gate assemblies from DNA sequences. "
            "Fidelity prediction adapted from OMEGA (Freschlin et al., 2025)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python design_fragments.py sequences.fasta --upstream_site AATG --downstream_site TTAG --primers primers.csv --ligation_data FileS04_T4_18h_37C.csv
  python design_fragments.py sequences.csv --upstream_site AATG --downstream_site TTAG --primers primers.csv --ligation_data FileS04_T4_18h_37C.csv --output_dir ./results
        """
    )

    # Required
    parser.add_argument('input', type=str,
                        help='Input sequences (.fasta, .fa, or .csv)')
    parser.add_argument('--upstream_site', type=str, required=True,
                        help="4bp 5' backbone Golden Gate overhang (e.g. AATG)")
    parser.add_argument('--downstream_site', type=str, required=True,
                        help="4bp 3' backbone Golden Gate overhang (e.g. TTAG)")
    parser.add_argument('--primers', type=str, required=False, default=DEFAULT_PRIMERS,
                        help='Primer CSV file (columns: name, sequence)')
    parser.add_argument('--ligation_data', type=str, required=False, default=DEFAULT_LIGATION_DATA,
                        help='Potapov ligation frequency CSV (from OMEGA data/ligation_data/, recommend FileS04_T4_18h_37C.csv)')

    # Optional
    parser.add_argument('--output_dir', type=str, default='./results',
                        help='Output directory (default: ./results)')
    parser.add_argument('--csv_id_column', type=str, default='id',
                        help='ID column name for CSV input (default: id)')
    parser.add_argument('--csv_sequence_column', type=str, default='sequence',
                        help='Sequence column name for CSV input (default: sequence)')

    args = parser.parse_args()
    main(args)
