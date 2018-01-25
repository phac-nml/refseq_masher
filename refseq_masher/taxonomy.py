"""NCBI Taxonomy information assignment

All taxonomic information for all unique NCBI Taxonomy UIDs of RefSeq genomes 
in the Mash RefSeq sketch database

"""
import logging
from typing import List

import pandas as pd

from refseq_masher.const import NCBI_TAXID_INFO


def get_full_taxonomy(l: List[str]) -> str:
    """From an ordered list of taxonomic classifications, return a string of concatenated taxonomic info 
    
    Remove some of the redundancy between classifications, e.g.
    
    ["Salmonella enterica", "Salmonella enterica subsp. enterica"]
    
    is concatenated to 
    
    "Salmonella enterica; subsp. enterica" rather than "Salmonella enterica; Salmonella enterica subsp. enterica"
    
    Args:
        l: list of taxonomic classifications
    
    Returns:
        (str): concatenated taxonomic classifications with some redundancy removed
    """
    out = [l[0]]
    for i in range(1, len(l)):
        prev = l[i - 1]
        curr = l[i]
        out.append(curr.replace(prev, '').strip())
    return '; '.join(out)


def merge_ncbi_taxonomy_info(dfmash):
    """Merge/join NCBI Taxonomy info with Mash results table 
    
    Merge/join on `taxid` (NCBI taxonomy UID)
    
    Args:
        dfmash: Mash results dataframe 
    
    Returns:
        (pd.DataFrame): dataframe with Mash results and taxonomy information
    """
    df = NCBI_TAXID_INFO[NCBI_TAXID_INFO.query_taxid.isin(dfmash.taxid)]
    if df.shape[0] > 0:
        logging.debug(df)
        df_has_rank = df[~pd.isnull(df['rank'])]
        dicts = []
        for taxid in dfmash.taxid.unique():
            d = {'taxid': taxid}
            df_has_rank_taxid = df_has_rank[df_has_rank['query_taxid'] == taxid]
            for _, r in df_has_rank_taxid.iterrows():
                d['taxonomic_{}'.format(r['rank'])] = r['name']
            df_matching_taxid = df[df['query_taxid'] == taxid]
            # highest resolution taxonomic classification is the last entry for that taxid so get last row
            _, row = df_matching_taxid[~df_matching_taxid['query_taxid'].duplicated(keep='last')].iterrows().__next__()
            d['top_taxonomy_name'] = row['name']
            # build string of concatenated all taxonomic info with some redundancy between successive terms 
            d['full_taxonomy'] = get_full_taxonomy(list(df_matching_taxid['name']))
            dicts.append(d)
        df_taxid_info = pd.DataFrame(dicts)
        dfmerge = pd.merge(dfmash, df_taxid_info, how='left', on='taxid')
        return dfmerge

    return dfmash