import argparse
import os


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, seq_lines = None, []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def main(args):
    os.makedirs(args.outdir, exist_ok=True)

    written = 0
    for header, seq in parse_fasta(args.input):
        # Strip protein|name= prefix if already present to get the bare name
        if header.startswith("protein|name="):
            bare_name = header[len("protein|name="):]
        else:
            bare_name = header

        # Sanitize bare name for use as a filename
        filename = bare_name.replace("|", "_").replace(" ", "_").replace("/", "_")
        out_path = os.path.join(args.outdir, f"{filename}.fasta")
        with open(out_path, "w") as f:
            f.write(f">protein|name={bare_name}\n{seq}\n")
        written += 1

    print(f"\nDone. {written} file(s) written to '{args.outdir}/'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a multi-sequence FASTA into individual files named after each sequence header.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python split_fasta.py -i master.fasta -o ./sequences
        """
    )
    parser.add_argument("-i", "--input",  required=True, help="Input FASTA file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory (created if it doesn't exist)")
    args = parser.parse_args()
    main(args)