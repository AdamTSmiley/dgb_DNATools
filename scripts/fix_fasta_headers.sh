#!/bin/bash
# fix_fasta_headers.sh
# Usage: bash fix_fasta_headers.sh <directory>

DIR="${1:-.}"

find "$DIR" -name "*.fasta" | while read -r fasta; do
    if grep -q '^>"' "$fasta" 2>/dev/null || grep -q "^>[^|]" "$fasta" 2>/dev/null; then
        sed -i 's/^>\(.*\)/>protein|name=\1/' "$fasta"
        echo "Fixed: $fasta"
    fi
done

echo "Done."
