import os
import argparse
import requests
import time
import random
import tarfile
from typing import Sequence, Dict, Tuple, Union
from Bio import SeqIO


def list_files_with_extensions(dir, extensions):
    return [f for f in os.listdir(dir) if f.endswith(extensions)]


def runMMseqs2(
        out_path: str,
        sequences: Union[Sequence[str], str],
        use_env: bool = True,
        use_filter: bool = True,
        use_pairing: bool = False,
        host_url: str = 'https://api.colabfold.com'
        ) -> Sequence[str]:
    """ Computes MSAs by querying MMseqs2 API. """

    submission_endpoint = 'ticket/pair' if use_pairing else 'ticket/msa'

    def submit(seqs: Sequence[str], mode: str, N=101) -> Dict[str, str]:
        """ Submits a query of sequences to MMseqs2 API. """
        n, query = N, ''
        for seq in seqs:
            query += f'>{n}\n{seq}\n'
            n += 1

        res = requests.post(f'{host_url}/{submission_endpoint}',
                            data={'q': query, 'mode': mode}, verify=False)
        try:
            out = res.json()
        except ValueError:
            out = {'status': 'ERROR'}

        return out

    def status(ID: int) -> Dict[str, str]:
        """ Obtains the status of a submitted query. """
        res = requests.get(f'{host_url}/ticket/{ID}', verify=False)
        try:
            out = res.json()
        except ValueError:
            out = {'status': 'ERROR'}

        return out

    def download(ID: int, path: str) -> None:
        """ Downloads the completed MMseqs2 query. """
        res = requests.get(f'{host_url}/result/download/{ID}', verify=False)
        with open(path, 'wb') as out:
            out.write(res.content)

    # Make input sequence a list if not already.
    sequences = [sequences] if isinstance(sequences, str) else sequences

    # Set the mode for MMseqs2.
    if use_filter:
        mode = 'env' if use_env else 'all'
    else:
        mode = 'env-nofilter' if use_env else 'nofilter'

    if use_pairing:
        mode = "pairgreedy"
        if use_env:
            mode = mode + "-env"

    # Set up output path.
    os.makedirs(out_path, exist_ok=True)
    tar_gz_file = os.path.join(out_path, 'out.tar.gz')
    N, REDO = 101, True

    # Deduplicate and keep track of order.
    unique_seqs = []
    [unique_seqs.append(seq) for seq in sequences if seq not in unique_seqs]
    Ms = [N + unique_seqs.index(seq) for seq in sequences]

    # Call MMseqs2 API.
    if not os.path.isfile(tar_gz_file):
        while REDO:
            # Resubmit job until it goes through
            out = submit(seqs=unique_seqs, mode=mode, N=N)

            while out['status'] in ['UNKNOWN', 'RATELIMIT']:
                # Resubmit
                time.sleep(10 + 5 * (len(unique_seqs) // 50) + random.randint(0, 5))
                out = submit(seqs=unique_seqs, mode=mode, N=N)

            if out['status'] == 'ERROR':
                raise Exception('MMseqs2 API is giving errors. Please confirm '
                                'your input is a valid protein sequence. If '
                                'error persists, please try again in an hour.')

            if out['status'] == 'MAINTENANCE':
                raise Exception('MMseqs2 API is undergoing maintenance. Please '
                                'try again in a few minutes.')

            # Wait for job to finish
            ID = out['id']
            while out['status'] in ['UNKNOWN', 'RUNNING', 'PENDING']:
                time.sleep(5 + 5 * (len(unique_seqs) // 50) + random.randint(0, 5))
                out = status(ID)

            if out['status'] == 'COMPLETE':
                REDO = False

            if out['status'] == 'ERROR':
                REDO = False
                raise Exception('MMseqs2 API is giving errors. Please confirm '
                                'your input is a valid protein sequence. If '
                                'error persists, please try again in an hour.')

        # Download results
        download(ID, tar_gz_file)

    # Get and extract a list of .a3m files.
    if use_pairing:
        a3m_files = [os.path.join(out_path, 'pair.a3m')]
    else:
        a3m_files = [os.path.join(out_path, 'uniref.a3m')]
    if use_env:
        a3m_files.append(
            os.path.join(out_path, 'bfd.mgnify30.metaeuk30.smag30.a3m'))
    if not os.path.isfile(a3m_files[0]):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(out_path)

    # Gather .a3m lines.
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        with open(a3m_file, 'r') as f:
            for line in f:
                if len(line) > 0:
                    # Replace NULL values
                    if '\x00' in line:
                        line = line.replace('\x00', '')
                        update_M = True
                    if line.startswith('>') and update_M:
                        M = int(line[1:].rstrip())
                        update_M = False
                        if M not in a3m_lines:
                            a3m_lines[M] = []
                    a3m_lines[M].append(line)

    # Return results.
    return [''.join(a3m_lines[n]) for n in Ms]


def main(args: argparse.Namespace):
    os.makedirs(args.alignment_dir, exist_ok=True)

    tag_list = []
    seq_list = []
    for fasta_file in list_files_with_extensions(args.fasta_dir, (".fasta", ".fa")):
        fasta_path = os.path.join(args.fasta_dir, fasta_file)
        for record in SeqIO.parse(fasta_path, "fasta"):
            tag_list.append(record.id)
            seq_list.append(str(record.seq))

    for tag, seq in zip(tag_list, seq_list):
        tag_dir = os.path.join(args.alignment_dir, tag)
        a3m_lines = runMMseqs2(
            out_path=tag_dir,
            sequences=seq,
            use_env=args.use_env,
            use_pairing=args.use_pairing
        )

        for f in list_files_with_extensions(tag_dir, ".a3m"):
            fp = os.path.join(tag_dir, f)
            fp_new = fp + "2"  # Rename to .a3m2 so OpenFold won't detect
            os.rename(fp, fp_new)

        with open(os.path.join(tag_dir, 'mmseqs_msa.a3m'), 'w') as f:
            f.writelines(a3m_lines)

    print(f'Computed MSAs for tags: {tag_list}')
    print(f'Use "{args.alignment_dir}" for the "use_precomputed_alignments" flag.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch MSAs for protein sequences using the MMseqs2/ColabFold API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python get_mmseqs2_msa.py ./fastas
  python get_mmseqs2_msa.py ./fastas --alignment_dir ./alignments
  python get_mmseqs2_msa.py ./fastas --alignment_dir ./alignments --use_pairing
        """,
    )
    parser.add_argument(
        "fasta_dir", type=str,
        help="Path to directory containing FASTA files"
    )
    parser.add_argument(
        "--alignment_dir", type=str, default=os.path.join(os.getcwd(), "alignments"),
        help="Name of the directory in which to output the computed MSAs. Defaults to './alignments'"
    )
    parser.add_argument(
        "--use_env", action='store_true', default=True,
        help="Include the environmental sequence databases in search. Defaults to True"
    )
    parser.add_argument(
        "--use_pairing", action='store_true', default=False,
        help="Attempt to pair sequences from same species. Defaults to False."
    )
    args = parser.parse_args()

    main(args)