"""Pytest unit tests for the core module of guidefinder
"""
import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import guidefinder

def test_pam_pam():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    assert getattr(pamobj, "pam") == "NGG"

def test_pam_orientation():
    pamobj = guidefinder.core.Pam("GATN", "3prime")
    assert getattr(pamobj, "pam_orientation") == "3prime"

pamobj = guidefinder.core.Pam("NGG", "5prime")

def test_pam_rc():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    assert pamobj.reverse_complement().pam == "CCN"

def test_pam_extend_ambiguous_dna():
    pamobj = guidefinder.core.Pam("NGG", "5prime")
    exset = pamobj.extend_ambiguous_dna()
    assert exset == frozenset(["AGG", "CGG","GGG","TGG"])


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

def test_target():
    tl = guidefinder.core.Target(seq="ATGCACATGCACTGCTCCA",
                                  pam="NGG",
                                  strand="forward",
                                  pam_orientation="5prime",
                                  seqid="testseq1",
                                  start=10,
                                  stop=30)
    assert tl.seq == "ATGCACATGCACTGCTCCA"

targets = [guidefinder.core.Target(seq="ATGCACATGCACTGCTGGA",
                              pam="NGG",
                              strand="forward",
                              pam_orientation="5prime",
                              seqid="NC_002516",
                              start=10,
                              stop=30),
     guidefinder.core.Target(seq="ATGCAAATTCTTGTGCTCCA",
                                   pam="NGG",
                                   strand="forward",
                                   pam_orientation="5prime",
                                   seqid="NC_002516",
                                   start=50,
                                   stop=71)]


def test_find_unique_near_pam():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     levindist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    assert len(tl.unique_targets) == 2

def test_create_index():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     levindist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()

def test_get_neighbors():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     levindist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    assert tl.neighbors["ATGCAAATTCTTGTGCTCCA"]["neighbors"]["dist"][1] == 7

def test_export_bed():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     levindist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    df = tl.export_bed()
    print(df)

def test_get_genbank_features():
    gbfiles = ["test/test_data/Burkholderia_thailandensis_E264__ATCC_700388_133.gbk",
               "test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"]
    gb_df = guidefinder.core.get_genbank_features(gbfiles)
    print(gb_df)

def test_get_nearby_feature():
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     levindist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    df = tl.export_bed()
    gbfiles = ["test/test_data/Burkholderia_thailandensis_E264__ATCC_700388_133.gbk",
               "test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"]
    gb_df = guidefinder.core.get_genbank_features(gbfiles)
    down, up = guidefinder.core.get_nearby_feature(targets=df, features=gb_df)
    print(down)
    print(up)

def test_get_fastas(tmp_path):
    gbfiles = ["test/test_data/Burkholderia_thailandensis_E264__ATCC_700388_133.gbk",
               "test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk.gz"]
    guidefinder.core.get_fastas(gbfiles, tmp_path)

def test_merge_downstream_upstream(tmp_path):
    tl = guidefinder.core.TargetList(targets=targets,
                                     lcp=10,
                                     levindist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    df = tl.export_bed()
    gbfiles = ["test/test_data/Burkholderia_thailandensis_E264__ATCC_700388_133.gbk",
               "test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"]
    gb_df = guidefinder.core.get_genbank_features(gbfiles)
    down, up = guidefinder.core.get_nearby_feature(targets=df, features=gb_df)
    outputfilename = os.path.join(tmp_path, "out.txt")
    targetfile_columns = list(list(df.values())[0].keys())
    featurefile_columns = list(gb_df.columns)
    joined_columns = targetfile_columns + featurefile_columns
    joined_columns.append("distance")
    ofn = merge_downstream_upstream(downsfile=down,
                                    upsfile=up,
                                    columns_name=joined_columns,
                                    outputfilename=outputfilename)
