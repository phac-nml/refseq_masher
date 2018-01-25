import logging
import os
import re
from collections import defaultdict
from subprocess import Popen, PIPE
from typing import List, Tuple, Union, Optional, Any

import pandas as pd

from refseq_masher.const import REGEX_FASTA, REGEX_FASTQ
from .const import REGEX_FASTQ, REGEX_FASTA

NT_SUB = {x: y for x, y in zip('acgtrymkswhbvdnxACGTRYMKSWHBVDNX', 'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')}


def run_command(cmdlist: List[str], stdin: Optional[Any] = None, stderr: Optional[Any] = PIPE) -> (int, str, str):
    p = Popen(cmdlist,
              stdout=PIPE,
              stderr=stderr,
              stdin=stdin)
    stdout, stderr = p.communicate()
    exit_code = p.returncode
    if isinstance(stdout, bytes):
        stdout = stdout.decode()
    if isinstance(stderr, bytes):
        stderr = stderr.decode()
    return exit_code, stdout, stderr


def exc_exists(exc_name: str) -> bool:
    """Check if an executable exists

    Args:
        exc_name (str): Executable name or path (e.g. "blastn")

    Returns:
        bool: Does the executable exists in the user's $PATH?
    """
    cmd = ['which', exc_name]
    exit_code, stdout, stderr = run_command(cmd)
    if exit_code == 0:
        return True
    else:
        logging.warning('which exited with non-zero code {} with command "{}"'.format(exit_code, ' '.join(cmd)))
        logging.warning(stderr)
        return False


def sample_name_from_fasta_path(fasta_path: str) -> str:
    """Extract genome name from fasta filename

    Get the filename without directory and remove the file extension.

    Example:
        With fasta file path ``/path/to/genome_1.fasta``::

            fasta_path = '/path/to/genome_1.fasta'
            genome_name = genome_name_from_fasta_path(fasta_path)
            print(genome_name)
            # => "genome_1"

    Args:
        fasta_path (str): fasta file path

    Returns:
        str: genome name
    """
    filename = os.path.basename(fasta_path)
    filename = re.sub(r'\.gz$', '', filename)
    return re.sub(r'\.(fa|fas|fasta|fna|\w{1,})(\.gz)?$', '', filename)


def sample_name_from_fastq_paths(fastqs: List[str]) -> str:
    """Get the sample name from FASTQ file paths

    Group FASTQs based on base filename without file extensions or expected FASTQ file delimiter (e.g. `_1`/`_2`)

    Args:
        fastqs: FASTQ paths

    Returns:
        (str): sample name
    """
    grouped_fastqs = group_fastqs(fastqs)
    for fastq_paths, sample_name in grouped_fastqs:
        return sample_name


def group_fastqs(fastqs: List[str]) -> List[Tuple[List[str], str]]:
    """Group FASTQs based on common base filename

    For example, if you have 2 FASTQs:

    - reads_1.fastq
    - reads_2.fastq

    The common name would be `reads` and the files would be grouped based on that common name.

    Args:
        fastqs: FASTQ file paths

    Returns:
        list of grouped FASTQs grouped by common base filename
    """
    genome_fastqs = defaultdict(list)
    for fastq in fastqs:
        filename = os.path.basename(fastq)
        basefilename = re.sub(r'_\d', '', REGEX_FASTQ.sub(r'\1', filename))
        genome_fastqs[basefilename].append(fastq)
    return [(fastq_paths, sample_name) for sample_name, fastq_paths in genome_fastqs.items()]


def collect_fasta_from_dir(input_directory: str) -> List[Tuple[str, str]]:
    fastas = []
    for x in os.listdir(input_directory):
        full_file_path = os.path.abspath(os.path.join(input_directory, x))
        if os.path.isfile(full_file_path) and REGEX_FASTA.match(x):
            sample_name = sample_name_from_fasta_path(full_file_path)
            fastas.append((full_file_path, sample_name))
    return fastas


def collect_fastq_from_dir(input_directory):
    fastqs = []
    for x in os.listdir(input_directory):
        full_file_path = os.path.abspath(os.path.join(input_directory, x))
        if os.path.isfile(full_file_path) and REGEX_FASTQ.match(x):
            fastqs.append(full_file_path)
    if len(fastqs) > 0:
        logging.info('Found %s FASTQ files in %s',
                     len(fastqs),
                     input_directory)
        reads_from_dir = group_fastqs(fastqs)
        logging.info('Collected %s read sets from %s FASTQ files in %s',
                     len(reads_from_dir),
                     len(fastqs),
                     input_directory)
        return reads_from_dir
    return []


def collect_inputs(inputs: List[str]) -> Tuple[List[Tuple[str, str]], List[Tuple[List[str], str]]]:
    """Collect all input files for analysis

    Sample names are derived from the base filename with no extensions.
    Sequencing reads are paired if they share a common filename name without "_\d".
    Filepaths for contigs and reads files are collected from an input directory if provided.

    Args:
        inputs: paths to FASTA/FASTQ files

    Returns:
        List of (contig filename, sample name)
        List of ([reads filepaths], sample name)
    """
    contigs = []
    reads = []

    fastas = [x for x in inputs if REGEX_FASTA.match(x)]
    fastqs = [x for x in inputs if REGEX_FASTQ.match(x)]
    dirs = [x for x in inputs if os.path.isdir(x)]
    if len(fastas) > 0:
        for fasta_path in fastas:
            fasta_path = os.path.abspath(fasta_path)
            if os.path.exists(fasta_path):
                genome_name = sample_name_from_fasta_path(fasta_path)
                contigs.append((fasta_path, genome_name))
            else:
                logging.error('Input fasta "%s" does not exist!', fasta_path)
    if len(fastqs) > 0:
        grouped_fastqs = group_fastqs(fastqs)
        logging.info('Grouped %s fastqs into %s groups',
                     len(fastqs),
                     len(grouped_fastqs))
        reads += grouped_fastqs
    for d in dirs:
        fasta_from_dir = collect_fasta_from_dir(d)
        if len(fasta_from_dir) > 0:
            logging.info('Collected %s FASTA from dir "%s"', len(fasta_from_dir), d)
            contigs = contigs + fasta_from_dir
        fastq_from_dir = collect_fastq_from_dir(d)
        if len(fastq_from_dir) > 0:
            logging.info('Collected %s FASTQ from dir "%s"', len(fastq_from_dir), d)
            reads += fastq_from_dir
    logging.info('Collected %s FASTA inputs and %s read sets', len(contigs), len(reads))
    return contigs, reads


LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'


def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)
    return lvl


def order_output_columns(dfout: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    set_columns = set(dfout.columns)
    present_columns = [x for x in cols if x in set_columns]
    rest_columns = list(set_columns - set(present_columns))
    return dfout[present_columns + rest_columns]