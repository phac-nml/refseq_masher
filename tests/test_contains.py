# -*- coding: utf-8 -*-

from io import StringIO

from click.testing import CliRunner
import pytest
import pandas as pd

import refseq_masher.cli as cli


@pytest.fixture(scope='module')
def runner():
    return CliRunner()


def test_fasta(runner):
    fasta = 'tests/data/Se-Enteritidis.fasta'
    result = runner.invoke(cli.contains, ['--top-n-results', '1',
                                          '-p', '4',
                                          fasta])

    assert result.exit_code == 0, 'Exit code should be 0'

    expected_top_tax_name = 'Salmonella enterica subsp. enterica serovar Enteritidis'
    df = pd.read_table(StringIO(result.output))
    assert df.top_taxonomy_name.str.contains(expected_top_tax_name).all(), \
        'Top Mash RefSeq result should have "{}" in the top_taxonomy_name field'.format(expected_top_tax_name)
