# -*- coding: utf-8 -*-
# import logging
import re
from pkg_resources import resource_filename

import pandas as pd

from . import program_name

MASH_REFSEQ_MSH = resource_filename(program_name, 'data/RefSeqSketches.msh')
NCBI_TAXID_INFO_CSV = resource_filename(program_name, 'data/refseq-taxid-info.csv')


def read_taxid_info_csv():
    # logging.error('Reading NCBI Taxonomy ID info from %s', NCBI_TAXID_INFO_CSV)

    df = pd.read_csv(NCBI_TAXID_INFO_CSV)
    # logging.error('Read NCBI Taxonomy ID info table with %s rows, %s cols', df.shape[0], df.shape[1])
    return df

NCBI_TAXID_INFO = read_taxid_info_csv()

REGEX_FASTQ = re.compile(r'^(.+)\.(fastq|fq)(\.gz)?$')
REGEX_FASTA = re.compile(r'^.+\.(fasta|fa|fna|fas)(\.gz)?$')
MASH_DIST_ORDERED_COLUMNS = '''
sample
top_taxonomy_name
distance
pvalue
matching
full_taxonomy
taxonomic_subspecies
taxonomic_species
taxonomic_genus
taxonomic_family
taxonomic_order
taxonomic_class
taxonomic_phylum
taxonomic_kingdom
taxonomic_superkingdom
subspecies
serovar
plasmid
bioproject
biosample
taxid
assembly_accession
match_id
'''.strip().split('\n')
MASH_SCREEN_ORDERED_COLUMNS = '''
sample
top_taxonomy_name
identity
shared_hashes
median_multiplicity
pvalue
full_taxonomy
taxonomic_subspecies
taxonomic_species
taxonomic_genus
taxonomic_family
taxonomic_order
taxonomic_class
taxonomic_phylum
taxonomic_kingdom
taxonomic_superkingdom
subspecies
serovar
plasmid
bioproject
biosample
taxid
assembly_accession
match_id
'''.strip().split('\n')