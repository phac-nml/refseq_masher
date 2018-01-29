# RefSeq Masher

[![travis-ci](https://travis-ci.org/phac-nml/refseq_masher.svg?branch=master)](https://travis-ci.org/phac-nml/refseq_masher) 
[![pypi](https://badge.fury.io/py/refseq-masher.svg)](https://pypi.python.org/pypi/refseq_masher/)


Find what NCBI RefSeq genomes match or are contained within your sequence data using [Mash MinHash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x) with a Mash sketch database of 54,925 NCBI RefSeq Genomes.

## Installation

Easiest way to install `refseq_masher` and all its dependencies is with [Conda](https://conda.io/docs/) through the [BioConda](https://bioconda.github.io/) channel:

```
conda install -c bioconda refseq_masher
```

Otherwise you can install `refseq_masher` from [PyPI](https://pypi.python.org/pypi/refseq-masher) with `pip install refseq_masher`, but you would need to manually install [Mash v2.0+](https://github.com/marbl/Mash/releases).


### Dependencies

Other than Python 3.5/3.6, the only external dependency of `refseq_masher` is [Mash v2.0+](https://github.com/marbl/Mash/releases).


### Python dependencies

- Pandas
- NumPy
- Click
- pytest for running tests


## Usage

If you run `refseq_masher` without any arguments, you should see the following usage info:

```
Usage: refseq_masher [OPTIONS] COMMAND [ARGS]...

  Find the closest matching NCBI RefSeq genomes or the genomes contained in
  your contigs or reads.

Options:
  --version      Show the version and exit.
  -v, --verbose  Logging verbosity (-v for logging warnings; -vvv for logging
				 debug info)
  -h, --help     Show this message and exit.

Commands:
  contains  Find the NCBI RefSeq genomes contained in...
  matches   Find NCBI RefSeq genome matches for an input...
```

`refseq_masher` has 2 commands:

- `matches` for finding the closest NCBI RefSeq genome matches to your input sequences


- `contains` for finding what RefSeq genomes are contained within your input sequences
	- useful for finding what genomes may be contained within your metagenomic sample


### `matches` - find the closest matching NCBI RefSeq Genomes in your input sequences

```
Usage: refseq_masher matches [OPTIONS] INPUT...

  Find NCBI RefSeq genome matches for an input genome fasta file

  Input is expected to be one or more FASTA/FASTQ files or one or more
  directories containing FASTA/FASTQ files. Files can be Gzipped.

Options:
  --mash-bin TEXT                 Mash binary path (default="mash")
  -o, --output PATH               Output file path (default="-"/stdout)
  --output-type [tab|csv]         Output file type (tab|csv)
  -n, --top-n-results INTEGER     Output top N results sorted by distance in
								  ascending order (default=5)
  -m, --min-kmer-threshold INTEGER
								  Mash sketch of reads: "Minimum copies of
								  each k-mer required to pass noise filter for
								  reads" (default=8)
  -h, --help                      Show this message and exit.
```

#### Example

With the [FNA.GZ](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/329/025/GCF_000329025.1_ASM32902v1/GCF_000329025.1_ASM32902v1_genomic.fna.gz) file for Salmonella enterica subsp. enterica serovar Enteritidis str. [CHS44](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/329/025/GCF_000329025.1_ASM32902v1/):

```
# download sequence file
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/329/025/GCF_000329025.1_ASM32902v1/GCF_000329025.1_ASM32902v1_genomic.fna.gz

# find RefSeq matches
refseq_masher -vv matches GCF_000329025.1_ASM32902v1_genomic.fna.gz
```

Log:

```
2018-01-29 11:02:13,786 INFO: Collected 1 FASTA inputs and 0 read sets [in ...refseq_masher/refseq_masher/utils.py:185]
2018-01-29 11:02:13,786 INFO: Creating Mash sketch file for ...refseq_masher/GCF_000329025.1_ASM32902v1_genomic.fna.gz [in ...refseq_masher/refseq_masher/mash/sketch.py:24]
2018-01-29 11:02:14,055 INFO: Created Mash sketch file at "/tmp/GCF_000329025.1_ASM32902v1_genomic.msh" [in ...refseq_masher/refseq_masher/mash/sketch.py:40]
2018-01-29 11:02:14,613 INFO: Ran Mash dist successfully (output length=11647035). Parsing Mash dist output [in ...refseq_masher/refseq_masher/mash/dist.py:64]
2018-01-29 11:02:15,320 INFO: Parsed Mash dist output into Pandas DataFrame with 54924 rows [in ...refseq_masher/refseq_masher/mash/dist.py:67]
2018-01-29 11:02:15,321 INFO: Deleting temporary sketch file "/tmp/GCF_000329025.1_ASM32902v1_genomic.msh" [in ...refseq_masher/refseq_masher/mash/dist.py:72]
2018-01-29 11:02:15,321 INFO: Sketch file "/tmp/GCF_000329025.1_ASM32902v1_genomic.msh" deleted! [in ...refseq_masher/refseq_masher/mash/dist.py:74]
2018-01-29 11:02:15,322 INFO: Ran Mash dist on all input. Merging NCBI taxonomic information into results output. [in ...refseq_masher/refseq_masher/cli.py:88]
2018-01-29 11:02:15,323 INFO: Fetching all taxonomy info for 5 unique NCBI Taxonomy UIDs [in ...refseq_masher/refseq_masher/taxonomy.py:35]
2018-01-29 11:02:15,325 INFO: Dropping columns with all NA values (ncol=32) [in ...refseq_masher/refseq_masher/taxonomy.py:38]
2018-01-29 11:02:15,327 INFO: Columns with all NA values dropped (ncol=11) [in ...refseq_masher/refseq_masher/taxonomy.py:40]
2018-01-29 11:02:15,327 INFO: Merging Mash results with relevant taxonomic information [in ...refseq_masher/refseq_masher/taxonomy.py:41]
2018-01-29 11:02:15,329 INFO: Merged Mash results with taxonomy info [in ...refseq_masher/refseq_masher/taxonomy.py:43]
2018-01-29 11:02:15,329 INFO: Merged taxonomic info into results output [in ...refseq_masher/refseq_masher/cli.py:90]
2018-01-29 11:02:15,329 INFO: Reordering output columns [in ...refseq_masher/refseq_masher/cli.py:91]
2018-01-29 11:02:15,331 INFO: Writing output to stdout [in ...refseq_masher/refseq_masher/writers.py:16]
```


Output:

Table output to standard output:

```
sample  top_taxonomy_name       distance        pvalue  matching        full_taxonomy   taxonomic_subspecies    taxonomic_species       taxonomic_genus taxonomic_family     taxonomic_order taxonomic_class taxonomic_phylum        taxonomic_superkingdom  subspecies      serovar plasmid bioproject      biosample   taxid    assembly_accession      match_id
GCF_000329025.1_ASM32902v1_genomic      Salmonella enterica subsp. enterica serovar Enteritidis str. CHS44      0.0     0.0     400/400 Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Salmonella; enterica; subsp. enterica; serovar Enteritidis; str. CHS44  Salmonella enterica subsp. enterica  Salmonella enterica     Salmonella      Enterobacteriaceae      Enterobacterales        Gammaproteobacteria     Proteobacteria  Bacteria    enterica Enteritidis             PRJNA185053     SAMN01041154    702979  NZ_ALFF ./rcn/refseq-NZ-702979-PRJNA185053-SAMN01041154-NZ_ALFF-.-Salmonella_enterica_subsp._enterica_serovar_Enteritidis_str._CHS44.fna
... 
[truncated output]
```

The top match is *Salmonella enterica* subsp. enterica serovar Enteritidis str. CHS44 with a distance of 0.0 and 400/400 sketches matching, which is what we expected. There's other taxonomic information available in the results table that may be useful. 


### `contains` - find what NCBI RefSeq Genomes are contained in your input sequences

If you have a metagenomic sample or maybe a sample with some contamination, you may be interested in seeing what's in your sample. You can do this with `refseq_masher contains <INPUT>`.

```
Usage: refseq_masher contains [OPTIONS] INPUT...

  Find the NCBI RefSeq genomes contained in your sequence files using Mash
  Screen

  Input is expected to be one or more FASTA/FASTQ files or one or more
  directories containing FASTA/FASTQ files. Files can be Gzipped.

Options:
  --mash-bin TEXT              Mash binary path (default="mash")
  -o, --output PATH            Output file path (default="-"/stdout)
  --output-type [tab|csv]      Output file type (tab|csv)
  -n, --top-n-results INTEGER  Output top N results sorted by identity in
							   ascending order (default=0/all)
  -i, --min-identity FLOAT     Mash screen min identity to report
							   (default=0.9)
  -v, --max-pvalue FLOAT       Mash screen max p-value to report
							   (default=0.01)
  -p, --parallelism INTEGER    Mash screen parallelism; number of threads to
							   spawn (default=1)
  -h, --help                   Show this message and exit.
```

#### Example - metagenomic a sample SAMEA1877339

For this example, we're going to see what RefSeq genomes are contained within sample [SAMEA1877340](https://www.ebi.ac.uk/ena/data/view/SAMEA1877340) from BioProject [PRJEB1775](https://www.ebi.ac.uk/ena/data/view/PRJEB1775).


Description from BioProject PRJEB1775:

> Design, Setting and Patients Forty-five samples were selected from a set of fecal specimens obtained from patients with diarrhea during the 2011 outbreak of STEC O104:H4 in Germany. Samples were chosen to represent STEC-positive patients with a range of clinical conditions and colony counts together with a small number of patients with other infections (Campylobacter jejnuni, Clostridium difficile and Salmonella enterica). Samples were subjected to high-throughput sequencing on the Illumina MiSeq and HiSeq 2500, followed by bioinformatics analysis.


We're going to download the FASTQ files for [ERR260489](https://www.ebi.ac.uk/ena/data/view/ERR260489&display=html):

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR260/ERR260489/ERR260489_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR260/ERR260489/ERR260489_2.fastq.gz
```

We're going to run `refseq_masher` against these FASTQ files:

```bash
refseq_masher -vv contains --top-n-results 50 -p 12 -o containment-ERR260489.tab ERR260489_1.fastq.gz ERR260489_2.fastq.gz
```



Log:

```
2018-01-29 10:59:25,849 INFO: Grouped 2 fastqs into 1 groups [in ...refseq_masher/refseq_masher/utils.py:174]
2018-01-29 10:59:25,849 INFO: Collected 0 FASTA inputs and 1 read sets [in ...refseq_masher/refseq_masher/utils.py:185]
2018-01-29 10:59:25,849 INFO: Running Mash Screen with NCBI RefSeq sketch database against sample "ERR260489" with inputs: ['../ERR260489_1.fastq.gz', '../ERR260489_2.fastq.gz'] [in ...refseq_masher/refseq_masher/mash/screen.py:44]
Loading ...refseq_masher/refseq_masher/data/RefSeqSketches.msh...
   4669418 distinct hashes.
Streaming from 2 inputs...
   Estimated distinct k-mers in pool: 206836855
Summing shared...
Computing coverage medians...
Writing output...
2018-01-29 11:00:19,665 INFO: Ran Mash Screen on all input. Merging NCBI taxonomic information into results output. [in ...refseq_masher/refseq_masher/cli.py:134]
2018-01-29 11:00:19,666 INFO: Fetching all taxonomy info for 23 unique NCBI Taxonomy UIDs [in ...refseq_masher/refseq_masher/taxonomy.py:35]
2018-01-29 11:00:19,669 INFO: Dropping columns with all NA values (ncol=32) [in ...refseq_masher/refseq_masher/taxonomy.py:38]
2018-01-29 11:00:19,671 INFO: Columns with all NA values dropped (ncol=12) [in ...refseq_masher/refseq_masher/taxonomy.py:40]
2018-01-29 11:00:19,671 INFO: Merging Mash results with relevant taxonomic information [in ...refseq_masher/refseq_masher/taxonomy.py:41]
2018-01-29 11:00:19,674 INFO: Merged Mash results with taxonomy info [in ...refseq_masher/refseq_masher/taxonomy.py:43]
2018-01-29 11:00:19,674 INFO: Merged taxonomic information into results output [in ...refseq_masher/refseq_masher/cli.py:136]
2018-01-29 11:00:19,674 INFO: Reordering output columns [in ...refseq_masher/refseq_masher/cli.py:137]
2018-01-29 11:00:19,677 INFO: Wrote output to "containment-ERR260489.tab" [in ...refseq_masher/refseq_masher/writers.py:20]
```


Output:

```
sample	top_taxonomy_name	identity	shared_hashes	median_multiplicity	pvalue	full_taxonomy	taxonomic_subspecies	taxonomic_species	taxonomic_genus	taxonomic_family	taxonomic_order	taxonomic_class	taxonomic_phylum	taxonomic_superkingdom	subspecies	serovar	plasmid	bioproject	biosample	taxid	assembly_accession	match_id	taxonomic_species group	match_comment
ERR260489	Bacteroides fragilis	1.0	400/400	786	0.0	Bacteria; FCB group; Bacteroidetes/Chlorobi group; Bacteroidetes; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; fragilis		Bacteroides fragilis	Bacteroides	Bacteroidaceae	Bacteroidales	Bacteroidia	Bacteroidetes	Bacteria			pLV22a			817		./rcn/refseq-NG-817-.-.-.-pLV22a-Bacteroides_fragilis.fna		
... [1 row] ...
ERR260489	Escherichia coli O104:H4 str. E92/11	1.0	400/400	48	0.0	Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia; coli; O104:H4; str. E92/11		Escherichia coli	Escherichia	Enterobacteriaceae	Enterobacterales	Gammaproteobacteria	Proteobacteria	Bacteria			pE9211p3			1090927	NZ_AHAU	./rcn/refseq-NZ-1090927-.-.-NZ_AHAU-pE9211p3-Escherichia_coli_O104_H4_str._E92_11.fna		
... [3 rows] ...
ERR260489	Kingella kingae KKC2005004457	1.0	400/400	5	0.0	Bacteria; Proteobacteria; Betaproteobacteria; Neisseriales; Neisseriaceae; Kingella; kingae; KKC2005004457		Kingella kingae	Kingella	Neisseriaceae	Neisseriales	Betaproteobacteria	Proteobacteria	Bacteria			unnamed			1229911		./rcn/refseq-NG-1229911-.-.-.-unnamed-Kingella_kingae_KKC2005004457.fna		
ERR260489	Bacteroides cellulosilyticus WH2	0.9998440000000001	399/400	772	0.0	Bacteria; FCB group; Bacteroidetes/Chlorobi group; Bacteroidetes; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; cellulosilyticus; WH2		Bacteroides cellulosilyticus	Bacteroides	Bacteroidaceae	Bacteroidales	Bacteroidia	Bacteroidetes	Bacteria			pBWH2B			1268240	NZ_ATFI	./rcn/refseq-NZ-1268240-.-.-NZ_ATFI-pBWH2B-Bacteroides_cellulosilyticus_WH2.fna		
... [1 row] ...
ERR260489	Klebsiella pneumoniae	0.9998440000000001	399/400	4	0.0	Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Klebsiella; pneumoniae		Klebsiella pneumoniae	Klebsiella	Enterobacteriaceae	Enterobacterales	Gammaproteobacteria	Proteobacteria	Bacteria			pMRC151			573		./rcn/refseq-NG-573-.-.-.-pMRC151-Klebsiella_pneumoniae.fna		
... [37 rows] ...
```


Some of the top genomes contained in this sample are sorted by identity and median multiplicity are:

- *Bacteroides fragilis* - fully contained (400/400) and high multiplicity (768)
- *Escherichia coli* O104:H4 - fully contained (400/400) and median multiplicity of 48
- *Kingella kingae* - fully contained (400/400) and median multiplicity of 5
- *Klebsiella pneumoniae* - 399/400 sketches contained with median multiplicity of 4

So with Mash we are able to find that the sample contained the expected genomic data (especially *E. coli* O104:H4). 



## Legal 

Copyright Government of Canada 2017

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

## Contact

**Gary van Domselaar**: gary.vandomselaar@phac-aspc.gc.ca

