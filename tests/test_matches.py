# -*- coding: utf-8 -*-

from io import StringIO

from click.testing import CliRunner
import pytest
import pandas as pd

from refseq_masher import __version__
import refseq_masher.cli as cli


@pytest.fixture(scope='module')
def runner():
    return CliRunner()


def test_print_version(runner):
    result = runner.invoke(cli.cli, ['--version'])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_fasta(runner):
    fasta = 'tests/data/Se-Enteritidis.fasta'
    result = runner.invoke(cli.matches, [fasta])

    assert result.exit_code == 0, 'Exit code should be 0'
    expected_top_tax_name = 'Salmonella enterica subsp. enterica serovar Enteritidis'
    df = pd.read_table(StringIO(result.output))
    assert df.top_taxonomy_name.str.contains(expected_top_tax_name).all(), \
        'All top 5 Mash RefSeq results should have "{}" in the top_taxonomy_name field'.format(expected_top_tax_name)
