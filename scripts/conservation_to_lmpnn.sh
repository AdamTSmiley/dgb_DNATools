#!/bin/bash
# conservation_to_lmpnn.sh
# Parses a _conservation.txt file and outputs LigandMPNN --fixed_residues format.
#
# Examples:
#   bash conservation_to_lmpnn.sh -i myprotein_conservation.txt -p 50 -c A
#   bash conservation_to_lmpnn.sh -i myprotein_conservation.txt -p 30 -c A -o fixed.txt

# Default values
PERCENTILE=50
CHAIN="A"
OUTPUT=""
EXTRA_RESIDUES=""

usage() {
    echo "Usage: $0 -i <conservation.txt> [-p <30|50|70>] [-c <chain>] [-o <output.txt>]"
    echo ""
    echo "  -i  Path to _conservation.txt file (required)"
    echo "  -p  Conservation percentile: 30, 50, or 70 (default: 50)"
    echo "  -c  Chain ID letter for LigandMPNN format (default: A)"
    echo "  -o  Output file path (default: print to terminal)"
    echo "  -r  Extra residue numbers to fix, space-separated (e.g. '5 18 42')"
    exit 1
}

# Flags
while getopts "i:p:c:o:r:" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        p) PERCENTILE="$OPTARG" ;;
        c) CHAIN="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) EXTRA_RESIDUES="$OPTARG" ;;
        *) usage ;;
    esac
done

# Validation
if [ -z "$INPUT" ]; then
    echo "[ERROR] Input file is required (-i)."
    usage
fi

if [ ! -f "$INPUT" ]; then
    echo "[ERROR] File not found: $INPUT"
    exit 1
fi

# [[ ]] is bash's extended test syntax -- better for string comparisons
if [[ "$PERCENTILE" != "30" && "$PERCENTILE" != "50" && "$PERCENTILE" != "70" ]]; then
    echo "[ERROR] Percentile must be 30, 50, or 70. Got: $PERCENTILE"
    exit 1
fi

# Parse conservation file
POSITIONS_LINE=$(grep -A 1 "^## Top ${PERCENTILE}%" "$INPUT" | tail -1)

if [ -z "$POSITIONS_LINE" ]; then
    echo "[ERROR] Could not find a 'Top ${PERCENTILE}%' section in: $INPUT"
    exit 1
fi

# Extract position numbers from conservation file as a newline-separated list
CONSERVED_NUMS=$(echo "$POSITIONS_LINE" \
    | tr ',' '\n' \
    | tr -d ' ' \
    | grep -v '^$')

# Combine with any extra residues the user passed in
ALL_NUMS=$(printf "%s\n%s" "$CONSERVED_NUMS" "$(echo "$EXTRA_RESIDUES" | tr ' ' '\n')" \
    | grep -v '^$')

# Then prepend chain letter and collapse back to a space-separated string
FIXED_RESIDUES=$(echo "$ALL_NUMS" \
    | sort -un \
    | sed "s/^/${CHAIN}/" \
    | tr '\n' ' ' \
    | sed 's/ $//')

# Write output
if [ -n "$OUTPUT" ]; then
    echo "$FIXED_RESIDUES" > "$OUTPUT"
    echo "[Done] Written to: $OUTPUT"
    echo "  Residues: $FIXED_RESIDUES"
else
    echo "$FIXED_RESIDUES"
fi
