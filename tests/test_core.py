"""Pytest unit tests for the core module of GuideMaker
"""
import os
import pytest

import numpy as np
import pandas as pd
from Bio import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


from guidemaker import doench_predict
from guidemaker import cfd_score_calculator
import guidemaker
from guidemaker.definitions import ROOT_DIR

# from typing import List, Dict, Tuple, TypeVar, Generator
# import altair as alt
# from pathlib import Path

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
configpath = os.path.join(ROOT_DIR,"data","config_default.yaml")


#PamTarget Class

def test_pam_pam():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime", "hamming")
    assert getattr(pamobj, "pam") == "NGG"


def test_pam_orientation():
    pamobj = guidemaker.core.PamTarget("GATN", "3prime", "hamming")
    assert getattr(pamobj, "pam_orientation") == "3prime"


pamobj = guidemaker.core.PamTarget("NGG", "5prime","hamming")


def test_pam_find_targets_5p():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime", "hamming")
    testseq1 = [SeqRecord(Seq.Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAAGTACAAAGCACGTTATTAGATGGTAACAATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAAGTACAAAGCACGTTATTAGATGGTGGGAAC"), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, target_len=6)
    assert target['target'][0] == "ATGCAC"
    assert target['target'][1] == "TAACAA"


#  5'-[guide of 25 nt][exact pam, 3nt][next two]-3'


def test_pam_find_targets_3p():
    pamobj = guidemaker.core.PamTarget("NGG", "3prime", "hamming")
    testseq1 = [SeqRecord(Seq.Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAAGTACAAAGCACGTTATTAGATGGTAACAATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAAGTACAAAGCACGTTATTAGATGGTGGGAAC]"), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, target_len=6)
    assert target['target'][0] == "ATGATC"
    assert target['target'][1] == "ATTAGA"

def test_pam_find_targets_fullgenome():
    file =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.fasta")
    pamobj = guidemaker.core.PamTarget("NGG", "5prime","hamming")
    #gb = SeqIO.parse("forward.fasta", "fasta")
    gb = SeqIO.parse(file, "fasta")
    target = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    assert target['target'][0] == "AAATGGTACGTTATGTGTTA"

tardict = {'target': ['AAATGGTACGTTATGTGTTA','AAATGGTACGTTATGTGTTA','AACAGTAAAATGGTTTAATG'],
        'exact_pam': ["AGG","TGG","CGG"],
        'start': [35, 41, 158572],
        'stop': [55, 61, 158592],
        'strand': [True, True, False],   # forward =True, reverse = Fasle
        'pam_orientation': [False,False, False], # 5prime =True, 3prime = Fasle
        'target_seq30': ['TTAGGAAATGGTACGTTATGTGTTATAAGA', 'AATGGTACGTTATGTGTTATAAGAATTTCT', 'AACGGAACAGTAAAATGGTTTAATGATACA'],
        'seqid': ['AP009180.1','AP009180.2','AP009180.1'],
        'seedseq': [np.nan, np.nan, np.nan],
        'isseedduplicated': [np.nan, np.nan, np.nan],
        'hasrestrictionsite': [np.nan, np.nan, np.nan],
        'dtype': ['hamming','hamming','hamming']}
    

targets = pd.DataFrame(tardict)
targets = targets.astype({"target":'str', "exact_pam": 'category', "start": 'uint32', "stop": 'uint32',"strand": 'bool', "pam_orientation": 'bool',"seqid": 'category'})


# TargetProcessor Class
def test_check_restriction_enzymes():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lsr=10,
                                     editdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NGGTAB'])
    assert tl.targets['hasrestrictionsite'][0] == True


def test_find_unique_near_pam():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lsr=10,
                                     editdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NGGTAB'])
    tl.find_unique_near_pam()
    assert tl.targets[tl.targets['isseedduplicated'] == False].shape == (2,12)


def test_create_index():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lsr=10,
                                     editdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)



def test_get_neighbors():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lsr=10,
                                     editdist=2,
                                     knum=2)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    tl.get_neighbors(configpath=configpath)
    print(tl.neighbors)
    assert tl.neighbors["AAATGGTACGTTATGTGTTA"]["neighbors"]["dist"][1] == 12


def test_export_bed():
    tl = guidemaker.core.TargetProcessor(targets=targets,
                                     lsr=10,
                                     editdist=2,
                                     knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    tl.get_neighbors(configpath=configpath)
    df = tl.export_bed()
    assert df.shape == (2, 5)




def test_get_control_seqs():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime","hamming")
    file =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.fasta")
    gb = SeqIO.parse(file, "fasta")
    targets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=targets, lsr=10, editdist=2, knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    gb = SeqIO.parse(file, "fasta")
    data = tl.get_control_seqs(gb,length=20, n=100, num_threads=2, configpath=configpath)
    assert data[2].shape == (100, 3)


# Annotation class tests
filegbk =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.gbk")
tl = guidemaker.TargetProcessor(targets=targets, lsr=10, editdist=2, knum=2)
tl.check_restriction_enzymes(['NRAGCA'])
tl.find_unique_near_pam()
tl.create_index(configpath=configpath)
tl.get_neighbors(configpath=configpath)
tf_df = tl.export_bed()
anno = guidemaker.Annotation(annotation_list=[filegbk], annotation_type="genbank",
                                           target_bed_df=tf_df)

def test_get_genbank_features():
    anno.get_annotation_features()
    assert 7 == len(anno.feature_dict)
    assert 182 == len(anno.genbank_bed_df)


def test_get_qualifiers():
    filegbk =os.path.join(TEST_DIR, "test_data", "Carsonella_ruddii.gbk")
    anno = guidemaker.Annotation(annotation_list=[filegbk], annotation_type="genbank",
                                 target_bed_df=tf_df)
    anno.get_annotation_features()
    anno._get_qualifiers(configpath=configpath)
    assert anno.qualifiers.shape == (182, 7)

def test_get_nearby_features(tmp_path):
    pamobj = guidemaker.core.PamTarget("NGG", "5prime","hamming")
    filegbk =os.path.join(TEST_DIR,"test_data", "Carsonella_ruddii.gbk")
    file =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.fasta")
    gb = SeqIO.parse(file, "fasta")
    targets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=targets, lsr=10, editdist=2, knum=2)
    tl.check_restriction_enzymes(['NGGRT'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    tl.get_neighbors(configpath=configpath)
    tf_df = tl.export_bed()
    anno = guidemaker.Annotation(annotation_list=[filegbk], annotation_type="genbank",
                                    target_bed_df=tf_df)
    anno.get_annotation_features()
    anno._get_nearby_features()
    assert anno.nearby.shape == (7074, 12)


def test_filter_features():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime","hamming")
    filegbk =os.path.join(TEST_DIR,"test_data", "Carsonella_ruddii.gbk")
    file =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.fasta")
    gb = SeqIO.parse(file, "fasta")
    pamtargets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=pamtargets, lsr=10, editdist=2, knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    tl.get_neighbors(configpath=configpath)
    tf_df = tl.export_bed()
    anno = guidemaker.Annotation(annotation_list=[filegbk], annotation_type="genbank",
                                 target_bed_df=tf_df)
    anno.get_annotation_features()
    anno._get_nearby_features()
    anno._filter_features()
    anno._get_qualifiers(configpath=configpath)
    anno._format_guide_table(tl)
    assert anno.pretty_df.shape == (900, 23)


def test_filterlocus():
    pamobj = guidemaker.core.PamTarget("NGG", "5prime","hamming")
    filegbk =os.path.join(TEST_DIR,"test_data", "Carsonella_ruddii.gbk")
    file =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.fasta")
    gb = SeqIO.parse(file, "fasta")
    pamtargets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=pamtargets, lsr=10, editdist=2, knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    tl.get_neighbors(configpath=configpath)
    tf_df = tl.export_bed()
    anno = guidemaker.core.Annotation(annotation_list=[filegbk],annotation_type="genbank",
                                        target_bed_df=tf_df)
    anno.get_annotation_features()
    anno._get_nearby_features()
    anno._filter_features()
    anno._get_qualifiers(configpath=configpath)
    anno._format_guide_table(tl)
    filter_prettydf = anno._filterlocus(attribute='locus_tag', filter_by_locus=['CRP_001'])
    assert filter_prettydf.shape == (4, 23)

# Function : get_fastas
@pytest.fixture(scope='session')
def test_get_fastas(tmpdir_factory):
    filegb =os.path.join(TEST_DIR,"test_data", "Carsonella_ruddii.gbk")
    guidemaker.core.get_fastas(filegb, tmpdir_factory)


# Function : extend_ambiguous_dna
def test_extend_ambiguous_dna():
    extend_seq = guidemaker.core.extend_ambiguous_dna('NGG')
    expected_seq = ['GGG', 'AGG', 'TGG', 'CGG']
    assert all([a == b for a, b in zip(extend_seq, expected_seq)])

def test_predict_guides():
    seqs = np.array(['GTACAAAGCACGTTATTAGATGGTGGGAAC', 'TCTAATCACGACAGCATCACTATTAGGCCG', 'TGAAATGTCTCTTATCTCTGTGTAAGGCTC'])
    exp_scores = np.array([[0.59383124], [0.28157765], [0.5276569]], dtype='float32')
    scores = doench_predict.predict(seqs)
    assert (exp_scores == scores).all()

def test_cdf_calc():
    result = cfd_score_calculator.calc_cfd("GCATGCACAGCTAGCATGCATGCAGCT", "GCATGCACAGCTAGCATGCATGCAGCG")
    assert abs(result - 0.176470588) < 0.0001


# doench_predict
def test_get_doench_efficiency_score():
    pamobj = guidemaker.core.PamTarget("NGG", "3prime","hamming")
    filegbk =os.path.join(TEST_DIR,"test_data", "Carsonella_ruddii.gbk")
    file =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.fasta")
    gb = SeqIO.parse(file, "fasta")
    pamtargets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=pamtargets, lsr=10, editdist=2, knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    tl.get_neighbors(configpath=configpath)
    tf_df = tl.export_bed()
    anno = guidemaker.core.Annotation(annotation_list=[filegbk],annotation_type="genbank",
                                        target_bed_df=tf_df)
    anno.get_annotation_features()
    anno._get_nearby_features()
    anno._filter_features()
    anno._get_qualifiers(configpath=configpath)
    anno._format_guide_table(tl)
    filter_pretty_30mer_df = anno._filterlocus(attribute = "locus_tag")
    doench_df = guidemaker.core.get_doench_efficiency_score(df=filter_pretty_30mer_df, pam_orientation=pamobj.pam_orientation)
    assert abs(doench_df.Efficiency[213] - 0.3245381) < 0.0001

def test_cfd_score():
    pamobj = guidemaker.core.PamTarget("NGG", "3prime","hamming")
    filegbk =os.path.join(TEST_DIR,"test_data", "Carsonella_ruddii.gbk")
    file =os.path.join(TEST_DIR, "test_data","Carsonella_ruddii.fasta")
    gb = SeqIO.parse(file, "fasta")
    pamtargets = pamobj.find_targets(seq_record_iter=gb, target_len=20)
    tl = guidemaker.core.TargetProcessor(targets=pamtargets, lsr=10, editdist=2, knum=10)
    tl.check_restriction_enzymes(['NRAGCA'])
    tl.find_unique_near_pam()
    tl.create_index(configpath=configpath)
    tl.get_neighbors(configpath=configpath)
    tf_df = tl.export_bed()
    anno = guidemaker.core.Annotation(annotation_list=[filegbk], annotation_type="genbank",
                                        target_bed_df=tf_df)
    anno.get_annotation_features()
    anno._get_nearby_features()
    anno._filter_features()
    anno._get_qualifiers(configpath=configpath)
    anno._format_guide_table(tl)
    filter_pretty_30mer_df = anno._filterlocus(attribute='locus_tag')
    cfd_df = guidemaker.core.cfd_score(df=filter_pretty_30mer_df)
    assert abs(cfd_df['Max CFD'][0] > 0)



def test_levin_dist():
    """Test that levin and hamming distance are being claculated correctly.

        The test sequence has three guides embedded in it:
        1. CGTAGCTAG[T]CACTAGCTGACA_GCA|AGG
        2. CGTAGCTAG[A]CACTAGCTGACA_GCA|AGG
        3. CGTAGCTAG[T]CACTAGCTGACTAGCA|AGG
        guide 2 has 1 substitution and guide 3 has 1 insertion relative to guide 1.
         The levin distances for seq 1 vs  [2, 3] are [1, 2], while the hamming distances for seq 1 vd/ [2, 3] are [1,16] 
    """
    distseq = [SeqRecord(Seq.Seq(
        "CGTAGCTAGTCACTAGCTGACAGCAAGGTTTTTCGTAGCTAGACACTAGCTGACAGCAAGGTTTTTTCGTAGCTAGTCACTAGCTGACTAGCAAGG"),
        id="distseq")]
    pt_levin = guidemaker.core.PamTarget("NGG", "3prime","levin")
    pt_hamming = guidemaker.core.PamTarget("NGG", "3prime","hamming")
    pd_levin = pt_levin.find_targets(seq_record_iter=distseq, target_len=20)
    pd_hamming = pt_hamming.find_targets(seq_record_iter=distseq, target_len=20)
    tp_levin =  guidemaker.core.TargetProcessor(targets=pd_levin, lsr=0, editdist=1, knum=3 )
    tp_hamming = guidemaker.core.TargetProcessor(targets=pd_hamming, lsr=0, editdist=1, knum=3)
    tp_levin.find_unique_near_pam()
    tp_hamming.find_unique_near_pam()
    tp_hamming.check_restriction_enzymes()
    tp_levin.check_restriction_enzymes()
    tp_levin.create_index(configpath=configpath)
    tp_hamming.create_index(configpath=configpath)
    tp_levin.get_neighbors(configpath=configpath)
    tp_hamming.get_neighbors(configpath=configpath)
    assert(tp_levin.neighbors['CTAGTCACTAGCTGACAGCA']['neighbors']['dist'] == [0, 1, 2])
    assert(tp_hamming.neighbors['CTAGTCACTAGCTGACAGCA']['neighbors']['dist'] == [0, 1, 16])
    
