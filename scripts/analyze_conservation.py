import os
import argparse
import numpy as np
from Bio import SeqIO
from collections import Counter


def parse_a3m(a3m_path: str) -> tuple[str, list[str]]:
    """
    Parse an .a3m file, returning the query sequence and aligned sequences
    with insertion columns (lowercase) stripped.
    """
    records = list(SeqIO.parse(a3m_path, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {a3m_path}")

    query_seq = str(records[0].seq)

    # Identify which columns are NOT insertions based on the query sequence
    # (query has no lowercase letters, so we keep columns that are uppercase or gap in query)
    keep_cols = [i for i, c in enumerate(query_seq) if not c.islower()]

    # Strip insertion columns from all sequences
    aligned = []
    for record in records:
        seq = str(record.seq)
        stripped = ''.join(seq[i] for i in keep_cols if i < len(seq))
        aligned.append(stripped.upper())

    return aligned[0], aligned


def compute_conservation(query: str, aligned: list[str]) -> list[dict]:
    """
    For each position in the query, compute amino acid frequencies
    and a conservation score (frequency of the most common AA).
    Gaps are excluded from frequency calculations.
    """
    results = []
    n_seqs = len(aligned)

    for pos in range(len(query)):
        col = [seq[pos] for seq in aligned if pos < len(seq)]
        # Exclude gaps
        col_no_gaps = [aa for aa in col if aa != '-']

        if not col_no_gaps:
            results.append({
                'position': pos + 1,  # 1-indexed
                'query_aa': query[pos],
                'top_aa': '-',
                'top_freq': 0.0,
                'counts': {},
                'n_seqs': 0
            })
            continue

        counts = Counter(col_no_gaps)
        top_aa, top_count = counts.most_common(1)[0]
        top_freq = top_count / len(col_no_gaps)

        results.append({
            'position': pos + 1,
            'query_aa': query[pos],
            'top_aa': top_aa,
            'top_freq': top_freq,
            'counts': dict(counts),
            'n_seqs': len(col_no_gaps)
        })

    return results


def select_fixed_positions(results: list[dict], percentiles: list[int] = [30, 50, 70]) -> dict:
    """
    Rank positions by conservation score and return the top X% for each
    specified percentile threshold, matching the paper's approach.
    """
    # Sort by conservation score descending
    ranked = sorted(results, key=lambda x: x['top_freq'], reverse=True)
    n = len(ranked)

    fixed = {}
    for pct in percentiles:
        n_fixed = max(1, int(np.ceil(n * pct / 100)))
        fixed[pct] = [r['position'] for r in ranked[:n_fixed]]

    return fixed, ranked


def write_output(results: list[dict], fixed: dict, ranked: list[dict], out_path: str, splits: bool, csv: bool):
    with open(out_path, 'w') as f:
        f.write("# Conservation Analysis\n")
        f.write(f"# Total positions: {len(results)}\n")
        f.write(f"# Total sequences in MSA: {results[0]['n_seqs'] if results else 0}\n\n")

        # Fixed position summary
        for pct, positions in fixed.items():
            f.write(f"## Top {pct}% most conserved positions ({len(positions)} residues):\n")
            f.write(', '.join(str(p) for p in sorted(positions)) + '\n\n')

        # Full ranked table
        f.write("## Full conservation ranking (position, query_aa, top_aa, conservation_score):\n")
        f.write(f"{'Rank':<6}{'Position':<10}{'Query_AA':<12}{'Top_AA':<10}{'Score':<10}{'N_seqs':<10}\n")
        f.write("-" * 58 + "\n")
        for rank, r in enumerate(ranked, 1):
            f.write(f"{rank:<6}{r['position']:<10}{r['query_aa']:<12}{r['top_aa']:<10}"
                    f"{r['top_freq']:.4f}    {r['n_seqs']:<10}\n")
    
    if splits:
        base = os.path.splitext(out_path)[0]
        for pct in [30, 50, 70]:
            with open(f'{base}_{pct}.txt', 'w') as f:
                f.write(', '.join(str(p) for p in sorted(fixed[pct])))
    if csv:
        csv_path = os.path.splitext(out_path)[0] + "_conservation.csv"
        with open(csv_path, 'w') as f:
            f.write("position,query_aa,top_aa,conservation_score,n_seqs\n")
            for r in sorted(ranked, key=lambda x: x['position']):
                f.write(f"{r['position']},{r['query_aa']},{r['top_aa']},{r['top_freq']:.4f},{r['n_seqs']}\n")


def main(args: argparse.Namespace):
    print(f"Parsing {args.a3m_file}...")
    query, aligned = parse_a3m(args.a3m_file)
    print(f"  Query length: {len(query)} residues")
    print(f"  Sequences in MSA: {len(aligned)}")

    print("Computing conservation scores...")
    results = compute_conservation(query, aligned)

    percentiles = [30, 50, 70]
    fixed, ranked = select_fixed_positions(results, percentiles)

    for pct, positions in fixed.items():
        print(f"  Top {pct}%: {len(positions)} fixed positions")

    out_path = args.output or os.path.splitext(args.a3m_file)[0] + "_conservation.txt"
    splits = args.splits
    csv = args.csv
    write_output(results, fixed, ranked, out_path, splits, csv)
    print(f"\nResults written to: {out_path}")

    if splits:
        splits_example = os.path.splitext(out_path)[0] + "_##.txt"
        print(f"\nSplits results written to: {splits_example}")
    
    if csv:
        csv_path = os.path.splitext(out_path)[0] + "_conservation.csv"
        print(f"\nCSV results written to: {csv_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze conservation of an MSA in .a3m format, "
                    "following the fixed residue selection approach from the TEV protease paper.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Examples:
    python analyze_conservation.py mmseqs_msa.a3m
    python analyze_conservation.py mmseqs_msa.a3m --output results/conservation.txt --splits True. --csv True
        """
    )
    parser.add_argument(
        "a3m_file", type=str,
        help="Path to the .a3m MSA file (e.g. alignments/myprotein/mmseqs_msa.a3m)"
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Path for output file. Defaults to <a3m_file>_conservation.txt"
    )
    parser.add_argument(
        "--splits", type=bool, default=False,
        help="Set True to split conserved positions out to unique files"
    )
    parser.add_argument(
        "--csv", type=bool, default=False,
        help="Set True to get a csv of conservation of each position across the protein sequence"
    )
    args = parser.parse_args()
    main(args)