import os
import sys
import argparse
import yaml
import pandas as pd
from dnachisel import *
from Bio import SeqIO

# Input parsing functions

def parse_fasta(filepath):
    
    sequences = []
    
    for record in SeqIO.parse(filepath, "fasta"):
        current_id = record.id
        current_seq = str(record.seq)
        
        sequences.append((current_id, current_seq))
    
    return sequences
        
def parse_csv(filepath, config):

    id_col = config.get("csv_id_column", "id")
    seq_col = config.get("csv_sequence_column", "sequence")
    
    df = pd.read_csv(filepath)
    
    for column in [id_col, seq_col]:
        if column not in df.columns:
            print(f"[ERROR] Column '{column}' not found in CSV.")
            print(f"Available columns: {list(df.columns)}")
            print("Set 'csv_id_col' and 'csv_sequence_column' in your config YAML.")
            sys.exit(1)
    
    sequences = list(zip(df[id_col].astype(str), df[seq_col].str.upper()))
    
    return sequences

# Sequence type ID and helper

def is_protein_sequence(seq):

    dna_characters = set("ATCG")
    
    return bool(set(seq.upper()) - dna_characters)

def prepare_dna_sequence(seq_id, seq, config):

    if is_protein_sequence(seq):
        print(f"[{seq_id}] Protein sequence -- reverse translating ...")
        dna_seq = reverse_translate(seq)
    else:
        dna_seq = seq
        if len(dna_seq) % 3 != 0:
            print(f"  [{seq_id}] WARNING: DNA sequence length {len(dna_seq)} is not divisible by 3. Skipping.")
            return None
    
    return dna_seq
            
# Constraints and objectives

def build_constraints(config):
    
    constraints = []

    constraints.append(EnforceTranslation())

    # GC content enforcement
    gc_config = config.get("enforce_gc_content", {})
    if gc_config.get("enabled", False):
        constraints.append(
            EnforceGCContent(
                mini=gc_config.get("min", 0.4),
                maxi=gc_config.get("max", 0.6),
                window=gc_config.get("window", 100),
            )
        )
        print(f"Enforcing GC content {gc_config.get('min', 0.4)}-{gc_config.get('max', 0.6)} in a {gc_config.get('window', 100)}bp window.")

    # Restriction site avoidance
    pattern_config = config.get("avoid_patterns", {})
    if pattern_config.get("enabled", False):
        for site in pattern_config.get("sites", []):
            constraints.append(AvoidPattern(site))
            print(f"Avoiding {site} restriction site.")

    # Hairpin avoidance
    hairpin_config = config.get("avoid_hairpins", {})
    if hairpin_config.get("enabled", False):
        constraints.append(
            AvoidHairpins(
                stem_size=hairpin_config.get("stem_size", 20),
                hairpin_window=hairpin_config.get("window", 200),
            )
        )
        print(f"Enforcing hairpin content {hairpin_config.get('stem_size', 20)}bp stem sizw in a {hairpin_config.get('window', 200)}bp window.")
        
    return constraints

def build_objectives(config):
  
    objectives = []

    # Codon optimization objective
    codon_config = config.get("codon_optimize", {})
    if codon_config.get("enabled", True):
        species = codon_config.get("species", "e_coli")
        objectives.append(CodonOptimize(species=species))
        print(f"Codon optimizing for '{species}'.")

    kmer_config = config.get("uniquify_kmers", {})
    if kmer_config.get("enabled", False):
        k = kmer_config.get("k", 9)
        include_rc = kmer_config.get("include_reverse_complement", True)
        objectives.append(
            UniquifyAllKmers(
                k=k,
                include_reverse_complement=include_rc,
            )
        )
        print(f"Making all kmers of {k}bp unique.")

    return objectives

# Sequence optimization

def optimize_sequence(seq_id, dna_seq, config):

    print(f"Optimizing: {seq_id} ({len(dna_seq)} bp)")

    constraints = build_constraints(config)
    objectives = build_objectives(config)
    
    add_stop = config.get("stop_codon", {})
    if add_stop.get("enabled", False):
        stop = add_stop.get("codon", "TAA")
        print(f"Adding {stop} stop codon.")
    else:
        stop = ""

    if not objectives:
        print(f"[{seq_id}] WARNING: No objectives defined. Check your config.")
        return None

    try:
        problem = DnaOptimizationProblem(
            sequence=dna_seq,
            constraints=constraints,
            objectives=objectives,
        )

        # Step 1: Resolve hard constraints first
        problem.resolve_constraints()

        # Step 2: Optimize towards the objectives
        problem.optimize()

        # Check that all constraints still pass after optimization
        if not problem.all_constraints_pass():
            print(f"[{seq_id}] WARNING: Some constraints failed after optimization.")
            print(problem.constraints_text_summary())

        return {
            "id": seq_id,
            "original_sequence": dna_seq,
            "optimized_sequence": problem.sequence + stop
        }

    except Exception as e:
        print(f"[{seq_id}] ERROR during optimization: {e}.")
        return None

# Writing output

def write_output(results, output_dir, output_name):

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{output_name}.csv")

    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)
    print(f"\n[Done] Results written to: {output_path}.")
    print(f"{len(results)} sequences optimized successfully.")
    return output_path

# Main function

def main(args):
    # load config
    print(f"\n[Config] Loading: {args.config}.")
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    # load seqs
    print(f"[Input] Loading: {args.input}.")
    ext = os.path.splitext(args.input)[1].lower()

    if ext in (".fasta", ".fa", ".faa"):
        raw_sequences = parse_fasta(args.input)
    elif ext == ".csv":
        raw_sequences = parse_csv(args.input, config)
    else:
        print(f"[ERROR] Unrecognized file extension: '{ext}'.")
        print("Supported formats: .fasta, .fa, .faa, .csv.")
        sys.exit(1)

    print(f"Loaded {len(raw_sequences)} sequences.\n")

    # seq optimize
    print("[Optimization] Starting...")
    results = []

    for seq_id, seq in raw_sequences:
        dna_seq = prepare_dna_sequence(seq_id, seq, config)
        if dna_seq is None:
            continue

        result = optimize_sequence(seq_id, dna_seq, config)
        if result is not None:
            results.append(result)

    if not results:
        print("\n[ERROR] No sequences were successfully optimized. Check your input and config.")
        sys.exit(1)

    # write out
    write_output(results, args.output_dir, args.output_name)

# Entry point

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Codon optimization pipeline using DNA Chisel.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python codon_optimize.py sequences.fasta config.yaml
  python codon_optimize.py sequences.csv config.yaml --output_dir ./results
  python codon_optimize.py sequences.fasta config.yaml --output_dir ./results --output_name my_run
        """,
    )

    # required arguments
    parser.add_argument(
        "input",
        help="Path to input file (.fasta, .fa, .faa, or .csv)",
    )
    parser.add_argument(
        "config",
        help="Path to YAML config file",
    )

    # optional arguments
    parser.add_argument(
        "--output_dir",
        default="./results",
        help="Directory to write output files (default: ./results)",
    )
    parser.add_argument(
        "--output_name",
        default="optimized_sequences",
        help="Base name for the output CSV file (default: optimized_sequences)",
    )

    args = parser.parse_args()
    main(args)








