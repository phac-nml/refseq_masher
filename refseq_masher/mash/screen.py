import logging
from typing import Union, List

import pandas as pd

from .const import MASH_REFSEQ_MSH
from .parser import mash_screen_output_to_dataframe
from .utils import run_command


def mash_screen_against_refseq(inputs: Union[str, List[str]],
                               mash_bin: str = 'mash',
                               sample_name: str = None,
                               max_pvalue: float = 0.01,
                               min_identity: float = 0.9,
                               parallelism: int = 1) -> pd.DataFrame:
    cmd_list = [mash_bin, 'screen',
                '-v', str(max_pvalue),
                '-p', str(parallelism),
                '-i', str(min_identity),
                MASH_REFSEQ_MSH]
    if isinstance(inputs, list):
        cmd_list += inputs
    elif isinstance(inputs, str):
        cmd_list.append(inputs)
    else:
        raise TypeError('Unexpected type "{}" for "inputs": {}'.format(type(inputs), inputs))
    logging.info('Running Mash Screen with NCBI RefSeq sketch database '
                 'against sample "%s" with inputs: %s', sample_name, inputs)
    exit_code, stdout, stderr = run_command(cmd_list, stderr=None)
    df = mash_screen_output_to_dataframe(stdout)
    df['sample'] = sample_name
    return df