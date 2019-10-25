import os
import tempfile
import subprocess
import pytest

from predictPAM import main

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
pamseq = 'ATCG'

def test_map_pam():
	fasta = os.path.join(TEST_DIR, "test_data", "genome.fasta")
	sobj = main.map_pam(threads=2)
