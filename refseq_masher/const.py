# -*- coding: utf-8 -*-

import re
from pkg_resources import resource_filename

from . import program_name

#: Mash sketch database with sketches from 54,925 RefSeq genomes package resource path
MASH_REFSEQ_MSH = resource_filename(program_name, 'data/RefSeqSketches.msh')
#: Regex for matching FASTQ filenames with optional .gz
REGEX_FASTQ = re.compile(r'^(.+)\.(fastq|fq)(\.gz)?$')
#: Regex for matching FASTA filenames with optional .gz
REGEX_FASTA = re.compile(r'^.+\.(fasta|fa|fna|fas)(\.gz)?$')
#: Ordered Mash dist and select taxonomy columns
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
#: Ordered Mash screen and select taxonomy columns
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
