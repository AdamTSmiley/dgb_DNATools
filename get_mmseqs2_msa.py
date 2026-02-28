import os
import argparse
import requests
import time
import random
import tarfile
from typing import Dict

from Bio import SeqIO


def runMMseqs2(
        out_path: str,
        sequence: str,
        use_filter: bool = True,
        host_url: str = 'https://api.colabfold.com'
        ) -> str:
    """ Computes an MSA for a single sequence by querying MMseqs2 API.
    Always searches environmental databases (env). """

    def submit(seq: str, mode: str, N: int = 101) -> Dict[str, str]:
        """ Submits a query sequence to MMseqs2 API. """
        query = f'>{N}\n{seq}\n'
        res = requests.post(f'{host_url}/ticket/msa',
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

    # Set the mode for MMseqs2 (always use env databases).
    mode = 'env' if use_filter else 'env-nofilter'

    # Set up output path.
    os.makedirs(out_path, exist_ok=True)
    tar_gz_file = os.path.join(out_path, 'out.tar.gz')
    N, REDO = 101, True

    # Call MMseqs2 API.
    if not os.path.isfile(tar_gz_file):
        while REDO:
            # Resubmit job until it goes through
            out = submit(seq=sequence, mode=mode, N=N)

            while out['status'] in ['UNKNOWN', 'RATELIMIT']:
                # Resubmit
                time.sleep(10 + random.randint(0, 5))
                out = submit(seq=sequence, mode=mode, N=N)

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
                time.sleep(5 + random.randint(0, 5))
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

    # Selectively extract only the two .a3m files we want.
    a3m_names = ['uniref.a3m', 'bfd.mgnify30.metaeuk30.smag30.a3m']
    a3m_files = [os.path.join(out_path, name) for name in a3m_names]

    if not os.path.isfile(a3m_files[0]):
        with tarfile.open(tar_gz_file) as tar_gz:
            for name in a3m_names:
                member = tar_gz.getmember(name)
                member.name = os.path.basename(member.name)
                tar_gz.extract(member, path=out_path)

    # Gather .a3m lines, skipping the query header in subsequent files
    # to avoid duplicate >101 headers in the merged output.
    a3m_lines = []
    for i, a3m_file in enumerate(a3m_files):
        skipped_header = (i == 0)  # first file keeps its header; subsequent files skip theirs
        with open(a3m_file, 'r') as f:
            for line in f:
                if len(line) > 0:
                    # Replace NULL values
                    if '\x00' in line:
                        line = line.replace('\x00', '')
                    if not skipped_header and line.startswith('>'):
                        skipped_header = True
                        continue
                    a3m_lines.append(line)

    return ''.join(a3m_lines)


def main(args: argparse.Namespace):
    os.makedirs(args.alignment_dir, exist_ok=True)

    records = list(SeqIO.parse(args.fasta_file, "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected exactly 1 sequence in {args.fasta_file}, found {len(records)}.")

    tag = records[0].id
    seq = str(records[0].seq)

    tag_dir = os.path.join(args.alignment_dir, tag)
    os.makedirs(tag_dir, exist_ok=True)
    a3m = runMMseqs2(
        out_path=tag_dir,
        sequence=seq,
    )

    with open(os.path.join(tag_dir, 'mmseqs_msa.a3m'), 'w') as f:
        f.write(a3m)

    print(f'Computed MSA for: {tag}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch an MSA for a protein sequence using the MMseqs2/ColabFold API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python get_mmseqs2_msa.py my_protein.fasta
  python get_mmseqs2_msa.py my_protein.fasta --alignment_dir ./alignments
        """,
    )
    parser.add_argument(
        "fasta_file", type=str,
        help="Path to a FASTA file containing a single protein sequence"
    )
    parser.add_argument(
        "--alignment_dir", type=str, default=os.path.join(os.getcwd(), "alignments"),
        help="Directory to write the computed MSA. Defaults to './alignments'"
    )
    args = parser.parse_args()

    main(args)