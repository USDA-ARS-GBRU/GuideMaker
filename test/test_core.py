"""Pytest unit tests for the core module of guidefinder
"""
import os
import pytest
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd
from typing import List, Set, Dict, Tuple
import guidemaker





# yaml loader
def test_load_parameters():
    yaml_dict_obj = guidemaker.core.load_parameters("test/test_data/parameters.yaml")
    expected_dict = {'NMSLIB': {'M': 8, 'efc': 64, 'post': 2, 'ef': 256},
                     'CONTROL': {'MINIMUM_HMDIST': 7,'CONTROL_SEARCH_MULTIPLE': [1, 2, 3, 4]},
                     'MINIMUM_PROPORTION': 0.5}
    assert(yaml_dict_obj.keys()==expected_dict.keys())


# Pam Class

def test_pam_pam():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime")
    assert getattr(pamobj, "pam") == "NGG"


def test_pam_orientation():
    pamobj = guidemaker.core.PamTarget("GATN", "3prime")
    assert getattr(pamobj, "pam_orientation") == "3prime"


pamobj = guidemaker.core.PamTarget("NGG", "5prime")


def test_pam_find_targets_5p():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime")
    testseq1 = [SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, target_len=6)
    assert target['target'][0] == "ATGCAC"
    assert target['target'][1] == "AGCAGT"

def test_pam_find_targets_3p():
    pamobj = guidemaker.core.PamTarget("NGG", "3prime")
    testseq1 = [SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, target_len=6)
    assert target['target'][0] == "ATGATC"
    assert target['target'][1] == "GCAGCT"

def test_pam_find_targets_fullgenome():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime")
    #gb = SeqIO.parse("forward.fasta", "fasta")
    gb = SeqIO.parse("test/test_data/Carsonella_ruddii.fasta", "fasta")
    target = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    assert target['target'][0] == "AAATGGTACGTTATGTGTTA"

tardict = {'target': ['ATGCACATGCACTGCTGGAT','ATGCAAATTCTTGTGCTCCA','CAAGCACTGCTGGATCACTG'],
        'exact_pam': ["AGG","TGG","CGG"],
        'start': [410, 1050, 1150],
        'stop': [430, 1070, 1170],
        'strand': [True, True, False],   # forward =True, reverse = Fasle
        'pam_orientation': [False,False, False], # 5prime =True, 3prime = Fasle
        'seqid': ['AP009180.1','AP009180.2','AP009180.1']}
    

targets = pd.DataFrame(tardict)
targets = targets.astype({"target":'str', "exact_pam": 'category', "start": 'uint32', "stop": 'uint32',"strand": 'bool', "pam_orientation": 'bool',"seqid": 'category'})



#guidemaker.getsize(targets)

# TargetProcessor Class


def test_check_restriction_enzymes():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NRAGCA'])
    assert tl.targets.shape == (2, 5)


def test_find_unique_near_pam():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    assert len(tl.unique_targets) == 2


def test_create_index():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index()


def test_get_neighbors():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    print(tl.neighbors)
    assert tl.neighbors["ATGCAAATTCTTGTGCTCCA"]["neighbors"]["dist"][1] == 12


def test_export_bed():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    df = tl.export_bed()
    assert df.shape == (2, 6)

def test_get_control_seqs():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime")
    gb = SeqIO.parse("test/test_data/Carsonella_ruddii.fasta", "fasta")
    targets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=targets, lu=10, hammingdist=2, knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index()
    gb = SeqIO.parse("test/test_data/Carsonella_ruddii.fasta", "fasta")
    data = tl.get_control_seqs(gb,length=20, n=100, num_threads=2)
    assert data[2].shape == (100, 3)


# Annotation class tests

tl = guidemaker.TargetProcessor(targets=targets, lu=10, hammingdist=2, knum=2)
tl.check_restriction_enzymes(['NRAGCA'])
tl.find_unique_near_pam()
tl.create_index()
tl.get_neighbors()
tf_df = tl.export_bed()
anno = guidemaker.Annotation(genbank_list=["test/test_data/Carsonella_ruddii.gbk"],
                                           target_bed_df=tf_df)

def test_get_genbank_features():
    anno._get_genbank_features()
    assert 7 == len(anno.feature_dict)
    assert 182 == len(anno.genbank_bed_df)


def test_get_qualifiers():
    anno = guidemaker.core.Annotation(genbank_list=["test/test_data/Carsonella_ruddii.gbk"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_qualifiers()
    assert anno.qualifiers.shape == (182, 7)

def test_get_nearby_features(tmp_path):
    pamobj = guidemaker.core.PamTarget("NGG", "5prime")
    gb = SeqIO.parse("test/test_data/Carsonella_ruddii.fasta", "fasta")
    targets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    for tar in targets:
        if len(tar.seq) < 20:
            print(str(tar.start) + ", " + str(tar.stop) + ", " + tar.seq)
    tl = guidemaker.core.TargetProcessor(targets=targets, lu=10, hammingdist=2, knum=2)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    tf_df = tl.export_bed()
    anno = guidemaker.core.Annotation(genbank_list=["test/test_data/Carsonella_ruddii.gbk"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_nearby_features()
    assert anno.nearby.shape == (6486, 14)


def test_filter_features():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime")
    gb = SeqIO.parse("test/test_data/Carsonella_ruddii.fasta", "fasta")
    pamtargets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=pamtargets, lu=10, hammingdist=2, knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    tf_df = tl.export_bed()
    anno = guidemaker.core.Annotation(genbank_list=["test/test_data/Carsonella_ruddii.gbk"],
                                        target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_nearby_features()
    anno._filter_features()
    anno._get_qualifiers()
    prettydf = anno._format_guide_table(tl)
    assert prettydf.shape == (783, 21)

# Function : get_fastas
def test_get_fastas(tmp_path):
    gbfiles = ["test/test_data/Burkholderia_thailandensis_E264__ATCC_700388_133.gbk.gz",
               "test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"]
    guidemaker.core.get_fastas(gbfiles, tmp_path)

# Function : extend_ambiguous_dna
def test_extend_ambiguous_dna():
    extend_seq = guidemaker.core.extend_ambiguous_dna('NGG')
    expected_seq = ['GGG', 'AGG', 'TGG', 'CGG']
    assert all([a == b for a, b in zip(extend_seq, expected_seq)])




