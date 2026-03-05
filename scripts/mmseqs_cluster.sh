#!/usr/bin/env bash
# cluster_top_reps.sh
# Self-search to compute mean pairwise identity, cluster at that threshold,
# then output a FASTA of representatives from the top 1-3 largest clusters.
#
# Usage: bash cluster_top_reps.sh input.fasta output.fasta

set -euo pipefail

INPUT="${1:-}"
OUTPUT="${2:-}"

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Usage: bash $(basename "$0") input.fasta output.fasta"
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: file not found: $INPUT"
    exit 1
fi

WORKDIR="cluster_work"
TMP="$WORKDIR/tmp"
M8="$WORKDIR/allvsall.m8"
CLUSTER_PREFIX="$WORKDIR/clusterRes"
TSV="${CLUSTER_PREFIX}_cluster.tsv"
REP_SEQS="${CLUSTER_PREFIX}_rep_seq.fasta"

mkdir -p "$WORKDIR" "$TMP"

# ── Step 1: Self-search ────────────────────────────────────────────────────────
echo "[1/4] Running all-vs-all self-search..."
mmseqs easy-search "$INPUT" "$INPUT" "$M8" "$TMP" \
    --min-seq-id 0.3 \
    -c 0.8 \
    --cov-mode 0 \
    --threads 8 \
    --format-output "query,target,fident" \
    -v 0

# ── Step 2: Compute mean pairwise identity (excluding self-hits) ───────────────
echo "[2/4] Computing mean pairwise identity..."
MEAN_ID=$(awk '$1 != $2 { sum += $3; n++ } END { if (n > 0) printf "%.2f", sum/n; else print "0.80" }' "$M8")
echo "      Mean pairwise identity: $MEAN_ID  →  using as --min-seq-id"

# ── Step 3: Cluster at mean identity ──────────────────────────────────────────
echo "[3/4] Clustering at --min-seq-id $MEAN_ID..."
mmseqs easy-cluster "$INPUT" "$CLUSTER_PREFIX" "$TMP" \
    --min-seq-id "$MEAN_ID" \
    -c 0.8 \
    --cov-mode 0 \
    --threads 8 \
    -v 0

# ── Step 4: Extract representatives from top clusters ─────────────────────────
echo "[4/4] Extracting top representatives..."

# Rank clusters by member count, write top 3 representative IDs to a temp file
TOP_REPS_FILE="$WORKDIR/top_reps.txt"
cut -f1 "$TSV" | sort | uniq -c | sort -rn | head -3 | awk '{print $2}' > "$TOP_REPS_FILE"
N_REPS=$(wc -l < "$TOP_REPS_FILE" | tr -d ' ')
echo "      Found $N_REPS cluster(s) — extracting representative(s)..."

# Pull matching sequences from rep_seq.fasta using the IDs file
awk '
NR==FNR { wanted[$1]=1; next }
/^>/ { id=substr($0,2); sub(/ .*/,"",id); capture=(id in wanted) }
{ if (capture) print }
' "$TOP_REPS_FILE" "$REP_SEQS" > "$OUTPUT"

# Confirm output
NSEQ=$(grep -c "^>" "$OUTPUT" || true)
echo ""
echo "Done. $NSEQ representative sequence(s) written to: $OUTPUT"
echo ""
echo "Cluster summary (rank, size, representative ID):"
cut -f1 "$TSV" | sort | uniq -c | sort -rn | head -3 | \
    awk '{printf "  Rank %d  |  %d members  |  %s\n", NR, $1, $2}'