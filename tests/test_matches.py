# -*- coding: utf-8 -*-
import os
from click.testing import CliRunner
import pytest

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
    fasta = 'tests/data/contigs.fasta'
    result = runner.invoke(cli.matches, [fasta])
    print(result)
    assert result.exit_code == 0
    print(result.output)