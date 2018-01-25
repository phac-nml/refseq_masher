# -*- coding: utf-8 -*-

import logging

import click
import pandas as pd

OUTPUT_TYPES = {'tab': '\t',
                'csv': ','}

def write_dataframe(dfout: pd.DataFrame,
                    output_path: str,
                    output_type: str) -> None:
    df_write_args = dict(sep=OUTPUT_TYPES[output_type], index=None)
    if output_path == '-':
        logging.info('Writing output to stdout')
        click.echo(dfout.to_csv(**df_write_args))
    else:
        dfout.to_csv(output_path, **df_write_args)
        logging.info('Wrote output to "%s"', output_path)


