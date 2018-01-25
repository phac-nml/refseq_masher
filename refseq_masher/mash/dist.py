import logging
import os
from io import StringIO
from subprocess import Popen, PIPE
from typing import Optional, List, Union

import pandas as pd

from refseq_masher.parser import parse_refseq_info
from refseq_masher.utils import sample_name_from_fasta_path, run_command, sample_name_from_fastq_paths
from .const import MASH_REFSEQ_MSH

MASH_DIST_4_COLUMNS = """
match_id
distance
pvalue
matching
""".strip().split('\n')

MASH_DIST_5_COLUMNS = """
match_id
query_id
distance
pvalue
matching
""".strip().split('\n')

MASH_SCREEN_COLUMNS = """
identity
shared_hashes
median_multiplicity
pvalue
match_id
match_comment
""".strip().split('\n')


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


def query_refseq(sketch_path: str, mash_bin: str = "mash") -> str:
    """Compute Mash distances of sketch file of genome fasta to RefSeq sketch DB.

    Args:
        mash_bin (str): Mash binary path
        sketch_path (str): Mash sketch file path or genome fasta file path

    Returns:
        (str): Mash STDOUT string
    """
    assert os.path.exists(sketch_path)
    cmd_list = [mash_bin,
                'dist',
                MASH_REFSEQ_MSH,
                sketch_path]
    exit_code, stdout, stderr = run_command(cmd_list)
    if exit_code != 0:
        raise Exception(
            'Could not run Mash dist. EXITCODE="{}" STDERR="{}" STDOUT="{}"'.format(exit_code, stderr, stdout))

    return stdout


def mash_dist_output_to_dataframe(mash_out: str) -> pd.DataFrame:
    """Mash dist stdout to Pandas DataFrame

    Args:
        mash_out (str): Mash dist stdout

    Returns:
        (pd.DataFrame): Mash dist table ordered by ascending distance
    """
    df = pd.read_table(StringIO(mash_out))
    ncols = df.shape[1]
    if ncols == 5:
        df.columns = MASH_DIST_5_COLUMNS
        df = df[MASH_DIST_4_COLUMNS]
    if ncols == 4:
        df.columns = MASH_DIST_4_COLUMNS
    df.sort_values(by='distance', ascending=True, inplace=True)
    match_ids = df.match_id
    refseq_matches = [parse_refseq_info(match_id=match_id) for match_id in match_ids]
    dfmatch = pd.DataFrame(refseq_matches)
    dfmerge = pd.merge(dfmatch, df, on='match_id')
    return dfmerge

def mash_screen_output_to_dataframe(mash_out: str) -> pd.DataFrame:
    df = pd.read_table(StringIO(mash_out))
    ncols = df.shape[1]
    df.columns = MASH_SCREEN_COLUMNS[:ncols]
    df.sort_values(by=['identity', 'median_multiplicity'], ascending=[False, False], inplace=True)
    match_ids = df.match_id
    refseq_matches = [parse_refseq_info(match_id=match_id) for match_id in match_ids]
    dfmatch = pd.DataFrame(refseq_matches)
    dfmerge = pd.merge(dfmatch, df, on='match_id')
    return dfmerge


def mash_dist_fasta_against_refseq(fasta_path: str,
                                   mash_bin: str = "mash",
                                   sample_name: Optional[str] = None,
                                   tmp_dir: str = "/tmp",
                                   k: int = 16,
                                   s: int = 400) -> pd.DataFrame:
    """

    :param fasta_path:
    :param mash_bin:
    :param sample_name:
    :param tmp_dir:
    :param k:
    :param s:
    :return:
    """
    sketch_path = None
    try:
        sketch_path = sketch_fasta(fasta_path,
                                   mash_bin=mash_bin,
                                   tmp_dir=tmp_dir,
                                   sample_name=sample_name,
                                   k=k,
                                   s=s)
        mashout = query_refseq(sketch_path, mash_bin=mash_bin)
        logging.info('Ran Mash dist successfully (output length=%s). Parsing Mash dist output', len(mashout))
        df_mash = mash_dist_output_to_dataframe(mashout)
        df_mash['sample'] = sample_name
        logging.info('Parsed Mash dist output into Pandas DataFrame with %s rows', df_mash.shape[0])
        logging.debug('df_mash: %s', df_mash.head(5))
        return df_mash
    finally:
        if sketch_path and os.path.exists(sketch_path):
            logging.info('Deleting temporary sketch file "%s"', sketch_path)
            os.remove(sketch_path)
            logging.info('Sketch file "%s" deleted!', sketch_path)


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


def mash_dist_fastq_against_refseq(fastqs: List[str],
                                   mash_bin: str = 'mash',
                                   sample_name: str = None,
                                   tmp_dir: str = '/tmp',
                                   k: int = 16,
                                   s: int = 400,
                                   m: int = 8) -> pd.DataFrame:
    """Compute Mash distances between input reads against all RefSeq genomes

    Args:
        fastqs: FASTQ paths
        mash_bin: Mash binary path
        sample_name: Sample name
        tmp_dir: Temporary working directory
        k: Mash kmer size
        s: Mash number of min-hashes
        m: Mash number of times a k-mer needs to be observed in order to be considered for Mash sketch DB

    Returns:
        (pd.DataFrame): Mash genomic distance results ordered by ascending distance
    """

    assert len(fastqs) > 0, "Must supply one or more FASTQ paths"
    sketch_path = None
    try:
        sketch_path = sketch_fastqs(fastqs,
                                    mash_bin=mash_bin,
                                    tmp_dir=tmp_dir,
                                    sample_name=sample_name,
                                    k=k,
                                    s=s,
                                    m=m)
        logging.info('Mash sketch database created for "%s" at "%s"', fastqs, sketch_path)
        logging.info('Querying Mash sketches "%s" against RefSeq sketch database', sketch_path)
        mashout = query_refseq(sketch_path, mash_bin=mash_bin)
        logging.info('Queried "%s" against RefSeq sketch database. Parsing into Pandas DataFrame', sketch_path)
        df_mash = mash_dist_output_to_dataframe(mashout)
        df_mash['sample'] = sample_name
        logging.info('Parsed Mash distance results into DataFrame with %s entries', df_mash.shape[0])
        logging.debug('df_mash %s', df_mash.head(5))
        return df_mash
    finally:
        if sketch_path and os.path.exists(sketch_path):
            logging.info('Deleting temporary sketch file "%s"', sketch_path)
            os.remove(sketch_path)
            logging.info('Sketch file "%s" deleted!', sketch_path)


def mash_screen_against_refseq(inputs: Union[str, List[str]],
                               mash_bin: str = 'mash',
                               sample_name: str = None,
                               max_pvalue: float = 0.01,
                               min_identity: float = 0.9,
                               parallelism: int = 1) -> pd.DataFrame:

    cmd_list = [mash_bin, 'screen',
                '-v', str(max_pvalue),
                '-p', str(parallelism),
                '-i', str(min_identity),
                MASH_REFSEQ_MSH]
    if isinstance(inputs, list):
        cmd_list += inputs
    elif isinstance(inputs, str):
        cmd_list.append(inputs)
    else:
        raise TypeError('Unexpected type "{}" for "inputs": {}'.format(type(inputs), inputs))
    logging.info('Running Mash Screen with NCBI RefSeq sketch database '
                 'against sample "%s" with inputs: %s', sample_name, inputs)
    exit_code, stdout, stderr = run_command(cmd_list, stderr=None)
    df = mash_screen_output_to_dataframe(stdout)
    df['sample'] = sample_name
    return df
