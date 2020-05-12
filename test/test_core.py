"""Pytest unit tests for the core module of guidefinder
"""
import os
import pytest
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from typing import List, Set, Dict, Tuple
import guidefinder


# Pam Class
from guidefinder import TargetList


def test_pam_pam():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    assert getattr(pamobj, "pam") == "NGG"


def test_pam_orientation():
    pamobj = guidefinder.core.Pam("GATN", "3prime")
    assert getattr(pamobj, "pam_orientation") == "3prime"


pamobj = guidefinder.core.Pam("NGG", "5prime")


def test_pam_rc():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    assert pamobj.reverse_complement() == "CCN"


def test_pam_extend_ambiguous_dna():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    exset = pamobj.extend_ambiguous_dna(pamobj.pam)
    assert exset == frozenset(["AGG", "CGG", "GGG", "TGG"])


def test_pam_find_targets_3f():
    pamobj = guidefinder.core.Pam("NGG", "3prime")
    testseq1 = SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTAGCAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")
    target = pamobj.find_targets(seqrecord_obj=testseq1, strand="forward", target_len=6)
    assert str(target[0].seq) == "ATGATC"


def test_pam_find_targets_5f():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    testseq1 = SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTAGCAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")
    target = pamobj.find_targets(seqrecord_obj=testseq1, strand="forward", target_len=6)
    assert str(target[0].seq) == "ATGCAC"


def test_pam_find_targets_3r():
    pamobj = guidefinder.core.Pam("NGG", "3prime")
    testseq1 = SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")
    target = pamobj.find_targets(seqrecord_obj=testseq1, strand="reverse", target_len=6)
    assert str(target[0].seq) == "GCAGCT"


def test_pam_find_targets_5r():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    testseq1 = SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")
    target = pamobj.find_targets(seqrecord_obj=testseq1, strand="reverse", target_len=6)
    assert str(target[0].seq) == "AGCAGT"


def test_pam_find_targets_5b():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    testseq1 = SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")
    target = pamobj.find_targets(seqrecord_obj=testseq1, strand="both", target_len=6)
    assert str(target[0].seq) == "ATGCAC"
    assert str(target[1].seq) == "AGCAGT"

# def test_pam_find_targets_fullgenome():
#     pamobj = guidefinder.core.Pam("NGG", "3prime")
#     gb = SeqIO.read("test/test_data/Pseudomonas_aeruginosa_PAO1_107.fasta", "fasta")
#     target = pamobj.find_targets(seqrecord_obj=gb, strand="forward", target_len=6)
#     assert str(target[0].seq) == "AGAGAC"


# Target Class
def test_target():
    tl = guidefinder.core.Target(seq="ATGCACATGCACTGCTCCA",
                                 exact_pam="NGG",
                                 strand="forward",
                                 pam_orientation="5prime",
                                 seqid="NC_002516",
                                 start=10,
                                 stop=30)
    assert tl.seq == "ATGCACATGCACTGCTCCA"


targets = [guidefinder.core.Target(seq="ATGCACATGCACTGCTGGAT",
                                   exact_pam="NGG", strand="forward",
                                   pam_orientation="5prime",
                                   seqid="NC_002516",
                                   start=410,
                                   stop=430),
          guidefinder.core.Target(seq="ATGCAAATTCTTGTGCTCCA",
                                  exact_pam="NGG",
                                  strand="forward",
                                  pam_orientation="5prime",
                                  seqid="NC_002516",
                                  start=1050,
                                  stop=1071)]


# TargetList Class

def test_find_unique_near_pam():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    assert len(tl.unique_targets) == 2


def test_create_index():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()


def test_get_neighbors():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    print(tl.neighbors)
    assert tl.neighbors["ATGCAAATTCTTGTGCTCCA"]["neighbors"]["dist"][1] == 6


def test_export_bed():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    df = tl.export_bed()


# Annotation class tests


tl: TargetList = guidefinder.core.TargetList(targets=targets, lcp=10, hammingdist=2, knum=2)
tl.find_unique_near_pam()
tl.create_index()
tl.get_neighbors()
tf_df = tl.export_bed()
anno = guidefinder.core.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"],
                                       target_bed_df=tf_df)

def test_get_genbank_features():
    anno._get_genbank_features()
    assert 7 == len(anno.feature_dict)
    assert 5584 == len(anno.genbank_bed_df)


def test_get_qualifiers():
    anno = guidefinder.core.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_qualifiers()
    assert anno.qualifiers.shape == (5584, 7)

def test_get_nearby_features(tmp_path):
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    gb = SeqIO.read("test/test_data/Pseudomonas_aeruginosa_PAO1_107.sample.fasta", "fasta")
    targets = pamobj.find_targets(seqrecord_obj=gb, strand="forward", target_len=20)
    tl = guidefinder.core.TargetList(targets=targets, lcp=10, hammingdist=2, knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    tf_df = tl.export_bed()
    anno = guidefinder.core.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_nearby_features()
    assert anno.nearby.shape == (10644, 14)

def test_filter_features():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    gb = SeqIO.read("test/test_data/Pseudomonas_aeruginosa_PAO1_107.sample.fasta", "fasta")
    targets = pamobj.find_targets(seqrecord_obj=gb, strand="both", target_len=20)
    tl = guidefinder.core.TargetList(targets=targets, lcp=10, hammingdist=2, knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    tf_df = tl.export_bed()
    anno = guidefinder.core.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_nearby_features()
    anno._filter_features()
    anno._get_qualifiers()
    prettydf = anno._format_guide_table(tl)
    return prettydf

def test_get_fastas(tmp_path):
    gbfiles = ["test/test_data/Burkholderia_thailandensis_E264__ATCC_700388_133.gbk.gz",
               "test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"]
    guidefinder.core.get_fastas(gbfiles, tmp_path)

