import logging
import os
from subprocess import Popen, PIPE
from typing import List

from ..utils import sample_name_from_fasta_path, run_command, sample_name_from_fastq_paths


def sketch_fasta(fasta_path, mash_bin="mash", tmp_dir="/tmp", sample_name=None, k=16, s=400):
    """Create Mash sketch file

    Args:
        fasta_path (str):
        mash_bin (str): Mash binary path
        tmp_dir (str): Temp user directory
        sample_name (str): Genome name
        fasta_path (str): Genome fasta file path
        k (int): kmer length
        s (int): number of sketches

    Returns:
        str: Mash sketch file path for genome fasta file
    """
    logging.info('Creating Mash sketch file for %s', fasta_path)
    if sample_name is None:
        sample_name = sample_name_from_fasta_path(fasta_path=fasta_path)

    msh_path = os.path.join(tmp_dir, sample_name + '.msh')
    cmd_list = [mash_bin,
                'sketch',
                '-k', str(k),
                '-s', str(s),
                '-o', msh_path,
                fasta_path]
    exit_code, stdout, stderr = run_command(cmd_list)
    if exit_code != 0:
        raise Exception(
            'Could not create Mash sketch. EXITCODE={} STDERR="{}" STDOUT="{}"'.format(exit_code, stderr, stdout))
    assert os.path.exists(msh_path), 'Mash sketch file does not exist at "{}"'.format(msh_path)
    logging.info('Created Mash sketch file at "%s"', msh_path)
    return msh_path


def sketch_fastqs(fastqs: List[str],
                  mash_bin: str = 'mash',
                  sample_name: str = None,
                  tmp_dir: str = '/tmp',
                  k: int = 16,
                  s: int = 400,
                  m: int = 8) -> str:
    """Create Mash sketch database from one or more FASTQ files

    Args:
        fastqs: list of FASTQ files (may be gzipped)
        mash_bin: Mash binary path
        sample_name: Sample name
        tmp_dir: Temporary working directory
        k: Mash kmer size
        s: Mash number of min-hashes
        m: Mash number of times a k-mer needs to be observed in order to be considered for Mash sketch DB

    Returns:
        (str): path to Mash sketch database for input FASTQs
    """
    if sample_name is None:
        sample_name = sample_name_from_fastq_paths(fastqs)
    p = Popen(['cat', *fastqs], stdout=PIPE)
    msh_path = os.path.join(tmp_dir, sample_name + '.msh')
    cmd_list = [mash_bin,
                'sketch',
                '-k', str(k),  # kmer size
                '-s', str(s),  # number of sketches
                '-m', str(m),  # min times a kmer needs to be observed to add to sketch DB
                '-o', msh_path,
                '-']
    logging.info('Creating Mash sketch file at "%s" from "%s"', msh_path, fastqs)
    exit_code, stdout, stderr = run_command(cmd_list, stdin=p.stdout, stderr=None)
    if exit_code != 0:
        raise Exception(
            'Could not create Mash sketch. EXITCODE={} STDERR="{}" STDOUT="{}"'.format(exit_code, stderr, stdout))
    assert os.path.exists(msh_path), 'Mash sketch file does not exist at "{}"'.format(msh_path)
    logging.info('Created Mash sketch file at "%s"', msh_path)
    return msh_path