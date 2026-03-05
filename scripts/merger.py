import argparse


def main(args):
    numbers = set()

    for path in args.inputs:
        with open(path, 'r') as f:
            for token in f.read().split():
                numbers.add(int(token))

    sorted_numbers = sorted(numbers)

    prefix = args.prefix or ""
    result = ' '.join(f"{prefix}{n}" for n in sorted_numbers)

    with open(args.output, 'w') as f:
        f.write(result)

    print(f"Written {len(sorted_numbers)} numbers to: {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge, deduplicate, and sort number files, with an optional prefix character.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python merge_numbers.py file1.txt file2.txt -o merged.txt
  python merge_numbers.py file1.txt file2.txt -o merged.txt --prefix A
        """
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more paths to input .txt files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path for the output file"
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Optional character to prepend to each number (e.g. A)"
    )
    args = parser.parse_args()
    main(args)