import argparse
import os
import re
import sys


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, seq_lines = None, []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]  # strip leading '>'
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def extract_id(header):
    """Return the integer id from 'id=N' in the header, or None if absent."""
    match = re.search(r'\bid=(\d+)\b', header)
    return int(match.group(1)) if match else None


def main(args):
    os.makedirs(args.outdir, exist_ok=True)

    written = 0
    skipped = 0

    for header, seq in parse_fasta(args.input):
        seq_id = extract_id(header)

        if seq_id is None:
            # Reference / scaffold sequence
            if args.include_ref:
                out_path = os.path.join(args.outdir, f"{args.basename}_ref.fasta")
                with open(out_path, "w") as f:
                    f.write(f">{header}\n{seq}\n")
                print(f"  Wrote reference → {out_path}")
                written += 1
            else:
                skipped += 1
            continue

        out_path = os.path.join(args.outdir, f"{args.basename}_{seq_id}.fasta")
        with open(out_path, "w") as f:
            f.write(f">{header}\n{seq.split(':')[0]}\n")
        written += 1

    print(f"\nDone. {written} file(s) written to '{args.outdir}/'  |  {skipped} record(s) skipped.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a LigandMPNN output FASTA into individual files per design (id=N).\n"
                    "The reference sequence (no id= tag) is skipped unless --include-ref is set.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python split_lmpnn_fasta.py -i lmpnn_model.fa -o ./designs -b lmpnn_model
  python split_lmpnn_fasta.py -i lmpnn_model.fa -o ./designs -b lmpnn_model --include-ref
        """
    )
    parser.add_argument("-i", "--input",    required=True,
                        help="Input FASTA file (LigandMPNN format)")
    parser.add_argument("-o", "--outdir",   required=True,
                        help="Output directory (created if it doesn't exist)")
    parser.add_argument("-b", "--basename", required=True,
                        help="Base name for output files (e.g. 'lmpnn_model')")
    parser.add_argument("--include-ref",    action="store_true",
                        help="Also write the reference sequence (no id=) as <basename>_ref.fasta")
    args = parser.parse_args()
    main(args)