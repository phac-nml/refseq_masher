# -*- coding: utf-8 -*-

"""NCBI Taxonomy information assignment

All taxonomic information for all unique NCBI Taxonomy UIDs of RefSeq genomes
in the Mash RefSeq sketch database is available in `NCBI_TAXID_INFO`. This info
is merged with Mash results on the `taxid` column.

"""

import logging
from pkg_resources import resource_filename

import pandas as pd

from . import program_name

#: NCBI taxonomy info table package resource path
NCBI_TAXID_INFO_CSV = resource_filename(program_name, 'data/ncbi_refseq_taxonomy_summary.csv')
#: Eagerly read the NCBI taxonomy info table for merging with Mash results
NCBI_TAXID_INFO = pd.read_csv(NCBI_TAXID_INFO_CSV, low_memory=False)


def merge_ncbi_taxonomy_info(dfmash: pd.DataFrame) -> pd.DataFrame:
    """Merge/join NCBI Taxonomy info with Mash results table

    Merge/join on `taxid` (NCBI taxonomy UID)

    Args:
        dfmash: Mash results dataframe

    Returns:
        (pd.DataFrame): dataframe with Mash results and taxonomy information
    """
    logging.info('Fetching all taxonomy info for %s unique NCBI Taxonomy UIDs', dfmash.taxid.unique().size)
    df_tax_info = NCBI_TAXID_INFO.loc[NCBI_TAXID_INFO.taxid.isin(dfmash.taxid), :]
    if df_tax_info.shape[0] > 0:
        logging.info('Dropping columns with all NA values (ncol=%s)', df_tax_info.shape[1])
        df_tax_info = df_tax_info.dropna(axis=1, how='all')
        logging.info('Columns with all NA values dropped (ncol=%s)', df_tax_info.shape[1])
        logging.info('Merging Mash results with relevant taxonomic information')
        dfmerge = pd.merge(dfmash, df_tax_info, how='left', on='taxid')
        logging.info('Merged Mash results with taxonomy info')
        return dfmerge
    else:
        logging.warning('No taxonomy info merged with Mash results!')

    return dfmash
