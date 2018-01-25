import logging
import os
from typing import Optional, List

import pandas as pd

from .sketch import sketch_fasta, sketch_fastqs
from .parser import mash_dist_output_to_dataframe
from ..utils import run_command
from ..const import MASH_REFSEQ_MSH


def mash_dist_refseq(sketch_path: str, mash_bin: str = "mash") -> str:
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


def fasta_vs_refseq(fasta_path: str,
                    mash_bin: str = "mash",
                    sample_name: Optional[str] = None,
                    tmp_dir: str = "/tmp",
                    k: int = 16,
                    s: int = 400) -> pd.DataFrame:
    """Compute Mash distances between input FASTA against all RefSeq genomes

    Args:
        fasta_path: FASTA file path
        mash_bin: Mash binary path
        sample_name: Sample name
        tmp_dir: Temporary working directory
        k: Mash kmer size
        s: Mash number of min-hashes

    Returns:
        (pd.DataFrame): Mash genomic distance results ordered by ascending distance
    """
    sketch_path = None
    try:
        sketch_path = sketch_fasta(fasta_path,
                                   mash_bin=mash_bin,
                                   tmp_dir=tmp_dir,
                                   sample_name=sample_name,
                                   k=k,
                                   s=s)
        mashout = mash_dist_refseq(sketch_path, mash_bin=mash_bin)
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


def fastq_vs_refseq(fastqs: List[str],
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
        mashout = mash_dist_refseq(sketch_path, mash_bin=mash_bin)
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
