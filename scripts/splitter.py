"""
Design Golden Gate assemblies from DNA sequences.

Takes sequences and fragments them into oligos bridged by optimal shared Golden Gate sites. 
Fidelity prediction is adapted from OMEGA (Freschlin et al., 2025 -- https://github.com/RomeroLab/omega).
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
from math import ceil, exp
import string

# Nuc info
ENZYME_SITES = {
    'BsaI':  {'prefix': 'GGTCTCA', 'suffix': 'AGAGACC'},   # GGTCTC + A padding
    'BsmBI': {'prefix': 'CGTCTCA', 'suffix': 'TGAGACG'},   # CGTCTC + A padding
}
ENZYME_OVERHEAD = 14  # 7bp prefix + 7bp suffix, same for both enzymes
OVERHANG_SIZE   = 4

# Seq info
OLIGO_MAX_LEN = 350
PRIMER_LEN = 20
FIXED_OVERHEAD = PRIMER_LEN * 2 + ENZYME_OVERHEAD
BACKBONE_OVERHANG = 4
MAX_END_CODING = OLIGO_MAX_LEN - FIXED_OVERHEAD - BACKBONE_OVERHANG
MAX_INT_CODING = OLIGO_MAX_LEN - FIXED_OVERHEAD
SEQ_MIN_LEN = 50

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
    """Skip sequences that are shorter than minimum length."""
    kept, skipped = [], []
    for seq_id, seq in sequences:
        if len(seq) >= SEQ_MIN_LEN:
            kept.append((seq_id, seq))
        else:
            skipped.append((seq_id, len(seq)))
    
    for seq_id, length in skipped:
        print(f"    [SKIPPED] {seq_id}: {length} bp -- shorter than cutoff (minimum {SEQ_MIN_LEN} bp).")

    return kept


def estimate_n_fragments(seq_len: int) -> int:
    """Calculate the minimum number of fragments needed to encode a sequence."""
    if seq_len <= MAX_END_CODING:
        return 1
    return ceil(seq_len / MAX_END_CODING)


# Backbone gg validation
def validate_backbone_sites(upstream: str, downstream: str) -> None:
    """Warn if backbone sites may cause assembly issues."""
    for name, site in [("Upstream", upstream), ("Downstream", downstream)]:
        if len(site) != 4:
            print(f"[ERROR] {name} site '{site}' is not 4bp.")
            sys.exit(1)
        if not all(c in 'ATCG' for c in site):
            print(f"[ERROR] {name} site '{site}' contains invalid characters.")
            sys.exit(1)
    
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
            "This will prevent directional cloning into a linearized backbone."
            "If you intend to circularize the fragment(s), this is fine.",
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


def find_split_position(
    seq: str,
    upstream: str,
    downstream: str,
    ligation_data: pd.DataFrame,
    disallowed_sites: set
) -> tuple | None:
    """Finds good positions to split large fragments into subfragments."""
    n_frags = estimate_n_fragments(len(seq))

    if n_frags == 1:
        fidelity = predict_fidelity(list(dict.fromkeys([upstream, downstream])), ligation_data)
        
        return ([], [], fidelity)
    
    if n_frags > 1:
        ideal_positions = [round(len(seq) * i / n_frags) for i in range(1, n_frags)]

        n = len(seq)
        placed_positions = []
        placed_sites = []

        for j, ideal_pos in enumerate(ideal_positions):
            ideal_gap = n / n_frags

            # Window bounds from oligo length constraints
            left_bound = placed_positions[-1] + OVERHANG_SIZE if placed_positions else 0
            right_bound = (placed_positions[-1] + MAX_INT_CODING if placed_positions else MAX_END_CODING) - OVERHANG_SIZE
            
            half_gap = int(ideal_gap // 2)
            p_min = max(left_bound, ideal_pos - half_gap)
            p_max = min(right_bound, ideal_pos + half_gap)
            
            if j == len(ideal_positions) - 1:
                p_min = max(p_min, n - MAX_END_CODING)
            
            scored = []
            for p in range(p_min, p_max + 1):
                ggsite = seq[p:p + OVERHANG_SIZE]
                if ggsite in disallowed_sites:
                    continue
                all_sites = list(dict.fromkeys([upstream] + placed_sites + [ggsite] + [downstream]))
                fidelity = predict_fidelity(all_sites, ligation_data)
                gap = p if j == 0 else p - placed_positions[-1]
                deviation = abs(gap - ideal_gap) / ideal_gap
                score = fidelity * exp(-5 * deviation ** 2)
                scored.append((p, ggsite, score, fidelity))

            if not scored:
                return None

            scored.sort(key=lambda x: -x[2])
            best_p, best_site, _, _ = scored[0]
            placed_positions.append(best_p)
            placed_sites.append(best_site)

        all_sites = list(dict.fromkeys([upstream] + placed_sites + [downstream]))
        final_fidelity = predict_fidelity(all_sites, ligation_data)
        return placed_positions, placed_sites, final_fidelity


# Oligo builder
def build_oligos(
    seq: str,
    upstream: str,
    downstream: str,
    positions: list,
    ggsites: list,
    fwd_primer: str,
    rev_primer: str,
    enzyme_prefix: str,
    enzyme_suffix: str,
) -> list:
    """Build orderable oligos from a fragmented sequence. """
    
    rev_rc = str(Seq(rev_primer).reverse_complement())

    if not positions:
        single_fragment = upstream + seq + downstream
        oligos = [fwd_primer + enzyme_prefix + single_fragment + enzyme_suffix + rev_rc]
    
    else:
        first_frag = upstream + seq[0:positions[0] + OVERHANG_SIZE]
        middle_fragments = [seq[positions[j-1]:positions[j] + OVERHANG_SIZE] for j in range(1, len(positions))]
        last_fragment = seq[positions[-1]:] + downstream

        all_fragments = [first_frag] + middle_fragments + [last_fragment]
        oligos = [fwd_primer + enzyme_prefix + fragment + enzyme_suffix + rev_rc for fragment in all_fragments]
    
    return oligos


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


def format_fragment_ranges(seq_len: int, positions: list) -> str:
    """Format fragment sequence ranges as a human-readable string, e.g. '1-198 | 195-387'."""
    boundaries = [0] + positions + [seq_len]
    ranges = []
    for i in range(len(boundaries) - 1):
        start = boundaries[i] + 1 if i == 0 else boundaries[i]  # 1-indexed
        end = boundaries[i+1] + (OVERHANG_SIZE if i < len(positions) else 0)
        ranges.append(f"{start}-{end}")
    return ' | '.join(ranges)


# Output
def write_outputs(results: list, output_dir: str) -> None:
    """Write order.csv and fidelity_report.csv."""
    os.makedirs(output_dir, exist_ok=True)

    # Order
    order_rows = []
    for r in results:
        for oligo, letter in zip(r['oligos'], string.ascii_uppercase):
            order_rows.append({'name': f"{r['id']}_{letter}", 'sequence': oligo})

    order_path = os.path.join(output_dir, 'order.csv')
    pd.DataFrame(order_rows).to_csv(order_path, index=False)
    print(f"  Orderable oligos written to:  {order_path}")

    # Report
    fidelity_rows = [{
        'id':                r['id'],
        'sequence_length':   r['sequence_length'],
        'n_fragments':       len(r['oligos']),
        'fragment_ranges':   format_fragment_ranges(r['sequence_length'], r['positions']),
        'shared_ggsites':    '|'.join(r['ggsites']),
        'oligo_lengths':     '|'.join(str(len(o)) for o in r['oligos']),
        'upstream_bbsite':   r['upstream_bbsite'],
        'downstream_bbsite': r['downstream_bbsite'],
        'fidelity':          round(r['fidelity'], 6),
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

    enzyme = args.enzyme
    if enzyme not in ENZYME_SITES:
        print(f"[ERROR] Unknown enzyme '{enzyme}'. Choose from: {list(ENZYME_SITES.keys())}")
        sys.exit(1)
    enzyme_prefix = ENZYME_SITES[enzyme]['prefix']
    enzyme_suffix = ENZYME_SITES[enzyme]['suffix']
    print(f"[Setup] Enzyme: {enzyme} (prefix: {enzyme_prefix}, suffix: {enzyme_suffix})")

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
    print(f"[Filter] Skipping sequences shorter than {SEQ_MIN_LEN} bp...")
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

        result = find_split_position(
            seq=seq,
            upstream=upstream,
            downstream=downstream,
            ligation_data=ligation_data,
            disallowed_sites=disallowed_sites
        )

        if result is None:
            print(f"  [FAILED] {seq_id}: no valid GG site found within oligo length constraints.")
            failed.append(seq_id)
            continue

        positions, ggsites, fidelity = result

        oligos = build_oligos(
            seq=seq,
            upstream=upstream,
            downstream=downstream,
            positions=positions,
            ggsites=ggsites,
            fwd_primer=fwd_seq,
            rev_primer=rev_seq,
            enzyme_prefix=enzyme_prefix,
            enzyme_suffix=enzyme_suffix
        )

        print(f"  {seq_id}: {len(oligos)} fragment(s), sites = {ggsites}, fidelity = {fidelity:.4f}")

        results.append({
            'id':                seq_id,
            'sequence_length':   len(seq),
            'positions':         positions,
            'ggsites':           ggsites,
            'upstream_bbsite':   upstream,
            'downstream_bbsite': downstream,
            'fidelity':          fidelity,
            'oligos':            oligos,
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
            "Design n-fragment Golden Gate assemblies from DNA sequences. "
            "Fidelity prediction adapted from OMEGA (Freschlin et al., 2025)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python splitter.py sequences.fasta --upstream_site AATG --downstream_site TTAG
  python splitter.py sequences.csv --upstream_site AATG --downstream_site TTAG --enzyme BsaI --primers primers.csv --ligation_data FileS04_T4_18h_37C.csv --output_dir ./results
        """
    )

    # Required
    parser.add_argument('input', type=str,
                        help='Input sequences (.fasta, .fa, or .csv)')
    parser.add_argument('--upstream_site', type=str, required=True,
                        help="4bp 5' backbone Golden Gate overhang (e.g. AATG)")
    parser.add_argument('--downstream_site', type=str, required=True,
                        help="4bp 3' backbone Golden Gate overhang (e.g. TTAG)")
    
    # Optional
    parser.add_argument('--enzyme', type=str, required=False, default='BsaI',
                    choices=['BsaI', 'BsmBI'],
                    help="Type IIS restriction enzyme to use (default: BsaI)")
    parser.add_argument('--primers', type=str, required=False, default=DEFAULT_PRIMERS,
                        help='Primer CSV file (columns: name, sequence)')
    parser.add_argument('--ligation_data', type=str, required=False, default=DEFAULT_LIGATION_DATA,
                        help='Potapov ligation frequency CSV (from OMEGA data/ligation_data/, recommend FileS04_T4_18h_37C.csv)')
    parser.add_argument('--output_dir', type=str, default='./results',
                        help='Output directory (default: ./results)')
    parser.add_argument('--csv_id_column', type=str, default='id',
                        help='ID column name for CSV input (default: id)')
    parser.add_argument('--csv_sequence_column', type=str, default='sequence',
                        help='Sequence column name for CSV input (default: sequence)')

    args = parser.parse_args()
    main(args)
