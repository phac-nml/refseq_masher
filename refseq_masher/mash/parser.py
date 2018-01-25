import logging
from io import StringIO
from typing import Optional

import pandas as pd

#: Sometimes Mash dist outputs 4 columns other times it outputs 5 columns
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
#: Mash screen output columns
MASH_SCREEN_COLUMNS = """
identity
shared_hashes
median_multiplicity
pvalue
match_id
match_comment
""".strip().split('\n')


def _no_periods(s: str) -> Optional[str]:
    return s if s != '.' else None


def parse_refseq_info(match_id: str) -> dict:
    """Parse a RefSeq Mash match_id

    For example from the following `match_id`:

    ./rcn/refseq-NZ-1147754-PRJNA224116-.-GCF_000313715.1-.-Salmonella_enterica_subsp._enterica_serovar_Enteritidis_str._LA5.fna

    If you split on '-' and ignoring the first two elements, you can extract, in order, the NCBI:

    - Taxonomy UID = 1147754
    - BioProject accession = PRJNA224116
    - BioSample accession = None
    - Genome accession = GCF_000313715.1
    - plasmid name = None
    - FNA filename (Salmonella_enterica_subsp._enterica_serovar_Enteritidis_str._LA5.fna)

    If "Salmonella" is found in the FNA filename, then serovar and subspecies will be parsed if present.
    For the example above, the subspecies would be "enterica" and the serovar would be "Enteritidis".

    Values with periods ('.') will be treated as None (null).

    Args:
        match_id (str): Mash RefSeq match_id with taxid, bioproject, full strain name, etc delimited by '-'

    Returns:
        (dict): parsed NCBI accession and other info
    """
    logging.debug('Parsing RefSeq info from "%s"', match_id)
    sp = match_id.split('-')
    _, prefix, taxid_str, bioproject, biosample, assembly_acc, plasmid, fullname = sp
    taxid = int(taxid_str)
    fullname = fullname.replace('.fna', '')
    serovar = None
    subsp = None
    if 'Salmonella' in fullname:
        if '_serovar_' in fullname:
            serovar = fullname.split('_serovar_')[-1].split('_str.')[0]
        if '_subsp._' in fullname:
            subsp = fullname.split('_subsp._')[-1].split('_')[0]

    return dict(match_id=match_id,
                taxid=taxid,
                biosample=_no_periods(biosample),
                bioproject=_no_periods(bioproject),
                assembly_accession=_no_periods(assembly_acc),
                plasmid=_no_periods(plasmid),
                serovar=serovar,
                subspecies=subsp)


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
    dfmatch = pd.DataFrame([parse_refseq_info(match_id=match_id) for match_id in match_ids])
    return pd.merge(dfmatch, df, on='match_id')


def mash_screen_output_to_dataframe(mash_out: str) -> pd.DataFrame:
    """Mash screen stdout to Pandas DataFrame

    Args:
        mash_out: Mash screen stdout

    Returns:
        (pd.DataFrame): Mash screen output table ordered by `identity` and `median_multiplicity` columns in descending order
    """
    df = pd.read_table(StringIO(mash_out))
    ncols = df.shape[1]
    df.columns = MASH_SCREEN_COLUMNS[:ncols]
    df.sort_values(by=['identity', 'median_multiplicity'], ascending=[False, False], inplace=True)
    match_ids = df.match_id
    refseq_matches = [parse_refseq_info(match_id=match_id) for match_id in match_ids]
    dfmatch = pd.DataFrame(refseq_matches)
    dfmerge = pd.merge(dfmatch, df, on='match_id')
    return dfmerge