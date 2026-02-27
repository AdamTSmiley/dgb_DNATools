import os
import argparse
import requests
import subprocess
import time
import random
import tarfile
import io
from typing import Sequence, Optional, Dict, Tuple, Union

from Bio import SeqIO


def list_files_with_extensions(dir, extensions):
    return [f for f in os.listdir(dir) if f.endswith(extensions)]


class HHSearch:
    def __init__(self, binary_path: str, databases: list):
        self.binary_path = binary_path
        self.databases = databases

    def query(self, a3m: str) -> str:
        db_args = []
        for db in self.databases:
            db_args += ["-d", db]
        cmd = [self.binary_path] + db_args + ["-i", "stdin", "-o", "stdout"]
        result = subprocess.run(
            cmd,
            input=a3m,
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout


def runMMseqs2(
        out_path: str,
        sequences: Union[Sequence[str], str],
        use_env: bool = True,
        use_filter: bool = True,
        use_templates: bool = False,
        num_templates: int = 20,
        use_pairing: bool = False,
        host_url: str = 'https://api.colabfold.com'
        ) -> Tuple[Sequence[str], Sequence[Optional[str]]]:
    """ Computes MSAs and templates by querying MMseqs2 API. """
    
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

    # Get templates if necessary. 
    if use_templates:
        templates = {}
        
        # Read MMseqs2 template outputs and sort templates based on query seq.
        with open(os.path.join(out_path, 'pdb70.m8'), 'r') as f:
            for line in f:
                p = line.rstrip().split()
                M, pdb = p[0], p[1]
                M = int(M)
                if M not in templates:
                    templates[M] = []
                templates[M].append(pdb)

        # Obtain template structures and data files
        template_paths = {}
        for k, TMPL in templates.items():
            TMPL_PATH = os.path.join(out_path, f'templates_{k}')
            if not os.path.isdir(TMPL_PATH):
                os.mkdir(TMPL_PATH)
                TMPL_LINE = ','.join(TMPL[:num_templates])
                # Obtain the .cif and data files for the templates
                os.system(
                    f'curl -s -L {host_url}/template/{TMPL_LINE} '
                    f'| tar xzf - -C {TMPL_PATH}/')
                # Rename data files
                os.system(
                    f'cp {TMPL_PATH}/pdb70_a3m.ffindex '
                    f'{TMPL_PATH}/pdb70_cs219.ffindex')
                os.system(f'touch {TMPL_PATH}/pdb70_cs219.ffdata')
            template_paths[k] = TMPL_PATH

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
    a3m_lines = [''.join(a3m_lines[n]) for n in Ms]

    if use_templates:
        template_paths_ = []
        for n in Ms:
            if n not in template_paths:
                template_paths_.append(None)
            else:
                template_paths_.append(template_paths[n])
        template_paths = template_paths_
    else:
        template_paths = []
        for n in Ms:
            template_paths.append(None)

    if isinstance(sequences, str):
        return (a3m_lines, template_paths[0])
    else:
        return (a3m_lines, template_paths)


def main(args: argparse.Namespace):
    os.makedirs(args.alignment_dir, exist_ok=True)
    os.makedirs(args.mmcif_dir, exist_ok=True)
    
    tag_list = []
    seq_list = []
    for fasta_file in list_files_with_extensions(args.fasta_dir, (".fasta", ".fa")):
        fasta_path = os.path.join(args.fasta_dir, fasta_file)
        for record in SeqIO.parse(fasta_path, "fasta"):
            tag_list.append(record.id)
            seq_list.append(str(record.seq))
    
    no_templates = []
    for tag, seqs in zip(tag_list, seq_list):
        tag_dir = os.path.join(args.alignment_dir, tag)
        a3m_lines, template_paths = runMMseqs2(
            out_path=tag_dir,
            sequences=seqs,
            use_env=args.use_env,
            use_templates=True,
            use_pairing=args.use_pairing
        )
        
        for f in list_files_with_extensions(tag_dir, ".a3m"):
            fp = os.path.join(tag_dir, f)
            fp_new = fp + "2" # Rename to .a3m2 file so OpenFold won't detect
            os.rename(fp, fp_new)
            
        with open(os.path.join(tag_dir, 'mmseqs_msa.a3m'), 'w') as f:
            f.writelines(a3m_lines)
            
        mmcif_dir = template_paths[0]
        if mmcif_dir is not None:
            hhsearch_pdb70_runner = HHSearch(
                binary_path="hhsearch", databases=[os.path.join(mmcif_dir, "pdb70")]
            )

            hhsearch_result = hhsearch_pdb70_runner.query(''.join(a3m_lines))
            pdb70_out_path = os.path.join(tag_dir, "pdb70_hits.hhr")
            with open(pdb70_out_path, "w") as f:
                f.write(hhsearch_result)
            
            for f in list_files_with_extensions(mmcif_dir, ".cif"):
                fp = os.path.join(mmcif_dir, f)
                fp_new = os.path.join(args.mmcif_dir, f)
                os.rename(fp, fp_new)
        else:
            no_templates.append(tag)
                    
    print(f'Computed MSAs and templates for tags: {tag_list}')
    print(f'Use "{args.alignment_dir}" for the "use_precomputed_alignments" flag.')
    if len(no_templates) != len(tag_list):
        print(f'Use "{args.mmcif_dir}" for the "mmcif_dir" flag.')
    if len(no_templates) > 0:
        print(f'No templates were found for tags: {no_templates}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta_dir", type=str,
        help="Path to directory containing FASTA files"
    )
    parser.add_argument(
        "--alignment_dir", type=str, default=os.path.join(os.getcwd(), "alignments"),
        help="Name of the directory in which to output the computed MSAs. Defaults to './alignments'"
    )
    parser.add_argument(
        "--mmcif_dir", type=str, default=os.path.join(os.getcwd(), "mmcif_dir"),
        help="Name of the directory in which to output the template mmCIFs. Defaults to './mmcif_dir'"
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