#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Main CLI script for refseq_masher with commands for running Mash dist
(matches) and Mash screen (contains) against a bundled NCBI RefSeq genomes
sketch database of 54,925 genomes that were Mash sketched at k=16, s=400.
"""

import click
import logging
from typing import List

import pandas as pd

import refseq_masher.mash.dist as mash_dist
import refseq_masher.mash.screen as mash_screen
from .const import MASH_DIST_ORDERED_COLUMNS, MASH_SCREEN_ORDERED_COLUMNS
from .taxonomy import merge_ncbi_taxonomy_info
from .utils import collect_inputs, init_console_logger, order_output_columns
from .writers import write_dataframe, OUTPUT_TYPES
from .utils import exc_exists

SCRIPT_NAME = 'refseq_masher'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def validate_mash_binary_exists(ctx, param, value):
    try:
        assert exc_exists(value)
        return value
    except AssertionError:
        raise click.BadParameter('Mash does not exist at "{}". '
                                 'Please install Mash to your $PATH'.format(value))


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
@click.option('-v', '--verbose', count=True,
              help="Logging verbosity (-v for logging warnings; -vvv for logging debug info)")
def cli(verbose):
    """Find the closest matching NCBI RefSeq genomes or the genomes contained in your contigs or reads.
    """
    lvl = init_console_logger(verbose)
    logging.debug('Initialized logging with %s level', lvl)


@cli.command()
@click.option('--mash-bin', default='mash',
              callback=validate_mash_binary_exists,
              help='Mash binary path (default="mash")')
@click.option('-o', '--output', default='-',
              type=click.Path(exists=False, writable=True),
              help='Output file path (default="-"/stdout)')
@click.option('--output-type', default='tab',
              type=click.Choice(OUTPUT_TYPES.keys()),
              help='Output file type ({})'.format('|'.join(OUTPUT_TYPES.keys())))
@click.option('-n', '--top-n-results', default=5, type=int,
              help='Output top N results sorted by distance in ascending order (default=5)')
@click.option('-m', '--min-kmer-threshold', type=int, default=8,
              help='Mash sketch of reads: "Minimum copies of each k-mer '
                   'required to pass noise filter for reads" (default=8)')
@click.option('-T', '--tmp-dir',
              type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True),
              default='/tmp',
              help='Temporary analysis files path (where to save temp Mash sketch file) (default="/tmp")')
@click.argument('input', type=click.Path(exists=True), nargs=-1, required=True)
def matches(mash_bin, output, output_type, top_n_results, min_kmer_threshold, tmp_dir, input):
    """Find NCBI RefSeq genome matches for an input genome fasta file

    Input is expected to be one or more FASTA/FASTQ files or one or more
    directories containing FASTA/FASTQ files. Files can be Gzipped.
    """
    dfs = []  # type: List[pd.DataFrame]
    contigs, reads = collect_inputs(input)
    logging.debug('contigs: %s', contigs)
    logging.debug('reads: %s', reads)
    for fasta_path, sample_name in contigs:
        df = mash_dist.fasta_vs_refseq(fasta_path,
                                       mash_bin=mash_bin,
                                       sample_name=sample_name,
                                       tmp_dir=tmp_dir)
        if top_n_results > 0:
            df = df.head(top_n_results)
        dfs.append(df)
    for fastq_paths, sample_name in reads:
        df = mash_dist.fastq_vs_refseq(fastq_paths,
                                       mash_bin=mash_bin,
                                       sample_name=sample_name,
                                       m=min_kmer_threshold,
                                       tmp_dir=tmp_dir)
        if top_n_results > 0:
            df = df.head(top_n_results)
        dfs.append(df)
    logging.info('Ran Mash dist on all input. Merging NCBI taxonomic information into results output.')
    dfout = merge_ncbi_taxonomy_info(pd.concat(dfs))
    logging.info('Merged taxonomic info into results output')
    logging.info('Reordering output columns')
    dfout = order_output_columns(dfout, MASH_DIST_ORDERED_COLUMNS)
    write_dataframe(dfout, output, output_type)


@cli.command()
@click.option('--mash-bin', default='mash',
              callback=validate_mash_binary_exists,
              help='Mash binary path (default="mash")')
@click.option('-o', '--output', default='-',
              type=click.Path(exists=False, writable=True),
              help='Output file path (default="-"/stdout)')
@click.option('--output-type', default='tab',
              type=click.Choice(OUTPUT_TYPES.keys()),
              help='Output file type ({})'.format('|'.join(OUTPUT_TYPES.keys())))
@click.option('-n', '--top-n-results', default=0, type=int,
              help='Output top N results sorted by identity in ascending order (default=0/all)')
@click.option('-i', '--min-identity', default=0.9, type=float,
              help='Mash screen min identity to report (default=0.9)')
@click.option('-v', '--max-pvalue', default=0.01, type=float,
              help='Mash screen max p-value to report (default=0.01)')
@click.option('-p', '--parallelism', default=1, type=int,
              help='Mash screen parallelism; number of threads to spawn (default=1)')
@click.argument('input', type=click.Path(exists=True), nargs=-1, required=True)
def contains(mash_bin, output, output_type, top_n_results, min_identity, max_pvalue, parallelism, input):
    """Find the NCBI RefSeq genomes contained in your sequence files using Mash Screen

    Input is expected to be one or more FASTA/FASTQ files or one or more
    directories containing FASTA/FASTQ files. Files can be Gzipped.
    """
    dfs = []
    contigs, reads = collect_inputs(input)

    for input_paths, sample_name in (contigs + reads):
        df = mash_screen.vs_refseq(inputs=input_paths,
                                   mash_bin=mash_bin,
                                   sample_name=sample_name,
                                   max_pvalue=max_pvalue,
                                   min_identity=min_identity,
                                   parallelism=parallelism)
        if top_n_results > 0:
            df = df.head(top_n_results)
        dfs.append(df)
    logging.info('Ran Mash Screen on all input. Merging NCBI taxonomic information into results output.')
    dfout = merge_ncbi_taxonomy_info(pd.concat(dfs))
    logging.info('Merged taxonomic information into results output')
    logging.info('Reordering output columns')
    dfout = order_output_columns(dfout, MASH_SCREEN_ORDERED_COLUMNS)
    write_dataframe(dfout, output_path=output, output_type=output_type)
