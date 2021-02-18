"""Pytest unit tests for the core module of guidefinder
"""
import os
import pytest
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from typing import List, Set, Dict, Tuple
import guidemaker


# Pam Class

def test_pam_pam():
    pamobj = guidemaker.core.Pam("NGG", "5prime")
    assert getattr(pamobj, "pam") == "NGG"


def test_pam_orientation():
    pamobj = guidemaker.core.Pam("GATN", "3prime")
    assert getattr(pamobj, "pam_orientation") == "3prime"


pamobj = guidemaker.core.Pam("NGG", "5prime")



def test_pam_find_targets_3f():
    pamobj = guidemaker.core.Pam("NGG", "3prime")
    testseq1 = [SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTAGCAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, strand="forward", target_len=6)
    assert str(target[0].seq) == "ATGATC"


def test_pam_find_targets_5f():
    pamobj = guidemaker.core.Pam("NGG", "5prime")
    testseq1 = [SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTAGCAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, strand="forward", target_len=6)
    assert str(target[0].seq) == "ATGCAC"


def test_pam_find_targets_3r():
    pamobj = guidemaker.core.Pam("NGG", "3prime")
    testseq1 = [SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, strand="reverse", target_len=6)
    assert str(target[0].seq) == "GCAGCT"


def test_pam_find_targets_5r():
    pamobj = guidemaker.core.Pam("NGG", "5prime")
    testseq1 = [SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, strand="reverse", target_len=6)
    assert str(target[0].seq) == "AGCAGT"


def test_pam_find_targets_5b():
    pamobj = guidemaker.core.Pam("NGG", "5prime")
    testseq1 = [SeqRecord(Seq("AATGATCTGGATGCACATGCACTGCTCCAAGCTGCATGAAAA",
                             alphabet=IUPAC.ambiguous_dna), id="testseq1")]
    target = pamobj.find_targets(seq_record_iter=testseq1, strand="both", target_len=6)
    assert str(target[0].seq) == "ATGCAC"
    assert str(target[1].seq) == "AGCAGT"

# def test_pam_find_targets_fullgenome():
#     pamobj = guidemaker.core.Pam("NGG", "3prime")
#     gb = SeqIO.parse("test/test_data/Pseudomonas_aeruginosa_PAO1_107.fasta", "fasta")
#     target = pamobj.find_targets(seq_record_iter=gb, strand="forward", target_len=6)
#     assert str(target[0].seq) == "AGAGAC"


# Target Class
def test_target():
    tl = guidemaker.core.Target(seq="ATGCACATGCACTGCTCCA",
                                 exact_pam="NGG",
                                 strand="forward",
                                 pam_orientation="5prime",
                                 seqid="NC_002516",
                                 start=10,
                                 stop=30)
    assert tl.seq == "ATGCACATGCACTGCTCCA"


targets = [guidemaker.Target(seq="ATGCACATGCACTGCTGGAT",
                                   exact_pam="NGG", strand="forward",
                                   pam_orientation="5prime",
                                   seqid="NC_002516",
                                   start=410,
                                   stop=430),
          guidemaker.Target(seq="ATGCAAATTCTTGTGCTCCA",
                                  exact_pam="NGG",
                                  strand="forward",
                                  pam_orientation="5prime",
                                  seqid="NC_002516",
                                  start=1050,
                                  stop=1071)]
guidemaker.getsize(targets)

targets = [guidemaker.Target(seq="ATGCACATGCACTGCTGGAT",
                                   exact_pam="NGG", strand=0,
                                   pam_orientation=True,
                                   seqid="NC_002516",
                                   start=410,
                                   stop=430),
          guidemaker.Target(seq="ATGCAAATTCTTGTGCTCCA",
                                  exact_pam="NGG",
                                  strand=0,
                                  pam_orientation=True,
                                  seqid="NC_002516",
                                  start=1050,
                                  stop=1071)]

guidemaker.getsize(targets)
# TargetList Class

def test_find_unique_near_pam():
    tl = guidemaker.core.TargetList(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    assert len(tl.unique_targets) == 2


def test_create_index():
    tl = guidemaker.core.TargetList(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()


def test_get_neighbors():
    tl = guidemaker.core.TargetList(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    print(tl.neighbors)
    assert tl.neighbors["ATGCAAATTCTTGTGCTCCA"]["neighbors"]["dist"][1] == 12


def test_export_bed():
    tl = guidemaker.core.TargetList(targets=targets,
                                     lu=10,
                                     hammingdist=2,
                                     knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    df = tl.export_bed()


# Annotation class tests


tl = guidemaker.TargetList(targets=targets, lu=10, hammingdist=2, knum=2)
tl.find_unique_near_pam()
tl.create_index()
tl.get_neighbors()
tf_df = tl.export_bed()
anno = guidemaker.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"],
                                           target_bed_df=tf_df)

def test_get_genbank_features():
    anno._get_genbank_features()
    assert 7 == len(anno.feature_dict)
    assert 5584 == len(anno.genbank_bed_df)


def test_get_qualifiers():
    anno = guidemaker.core.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_qualifiers()
    assert anno.qualifiers.shape == (5584, 7)

def test_get_nearby_features(tmp_path):
    pamobj = guidemaker.core.Pam("NGG", "5prime")
    gb = SeqIO.parse("test/test_data/Pseudomonas_aeruginosa_PAO1_107.sample.fasta", "fasta")
    targets = pamobj.find_targets(seq_record_iter=gb, strand="forward", target_len=20)
    for tar in targets:
        if len(tar.seq) < 20:
            print(str(tar.start) + ", " + str(tar.stop) + ", " + tar.seq)
    tl = guidemaker.core.TargetList(targets=targets, lu=10, hammingdist=2, knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    tf_df = tl.export_bed()
    anno = guidemaker.core.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_nearby_features()
    assert anno.nearby.shape == (10792, 14)

def test_get_control_seqs():
    pamobj = guidemaker.core.Pam("NGG", "5prime")
    gb = SeqIO.parse("test/test_data/Pseudomonas_aeruginosa_PAO1_107.sample.fasta", "fasta")
    targets = pamobj.find_targets(seq_record_iter=gb, strand="forward", target_len=20)
    tl = guidemaker.core.TargetList(targets=targets, lu=10, hammingdist=2, knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    gb = SeqIO.parse("test/test_data/Pseudomonas_aeruginosa_PAO1_107.sample.fasta", "fasta")
    data = tl.get_control_seqs(gb,length=20, n=100, num_threads=2)

def test_filter_features():
    pamobj = guidemaker.core.Pam("NGG", "5prime")
    gb = SeqIO.parse("test/test_data/Pseudomonas_aeruginosa_PAO1_107.sample.fasta", "fasta")
    pamtargets = pamobj.find_targets(seq_record_iter=gb, strand="both", target_len=20)
    tl = guidemaker.core.TargetList(targets=pamtargets, lu=10, hammingdist=2, knum=2)
    tl.find_unique_near_pam()
    tl.create_index()
    tl.get_neighbors()
    tf_df = tl.export_bed()
    anno = guidemaker.core.Annotation(genbank_list=["test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk.gz"],
                                       target_bed_df=tf_df)
    anno._get_genbank_features()
    anno._get_nearby_features()
    anno._filter_features()
    anno._get_qualifiers()
    prettydf = anno._format_guide_table(tl)
    assert prettydf.shape == (2141, 21)

def test_get_fastas(tmp_path):
    gbfiles = ["test/test_data/Burkholderia_thailandensis_E264__ATCC_700388_133.gbk.gz",
               "test/test_data/Pseudomonas_aeruginosa_PAO1_107.gbk"]
    guidemaker.core.get_fastas(gbfiles, tmp_path)


def findall(pam, seq):
            """ Find occurrence of substring(PAM) in a string(sequence)
            """
            i = 0
            try:
                while True:
                    i = seq.index(pam, i)
                    yield i
                    i += 1
            except ValueError:
                pass


from multiprocessing.dummy import Pool as ThreadPool
hay='ATATATATTATTATATTATTATATTATTATATTATTTTTTTTTTTTTTTATATTAT'
needle="AT"
aa = findall(needle, hay)
for i in aa:
    print(i)

pool = ThreadPool(10)
results = pool.map(findall(needle, hay))


with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    executor.map(findall(needle, hay), range(3))

def multi_threads(threads=3):
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        aa = executor.map(findall(needle, hay), range(threads))
        return aa

futures = multi_threads(threads=3)

found=[]
for future in futures:
            #print('type(future)=',type(future))
            for f in future:
                if f:
                    try:
                        found.append(f)
                    except:
                        print_exc()

foundsize = len(found)

############
import concurrent.futures as cf
import itertools
from time import time
from traceback import print_exc
from itertools import chain

def _findmatch(nmin, nmax, number):
    '''Function to find the occurrence of number in range nmin to nmax and return
       the found occurrences in a list.'''
    print('\n def _findmatch', nmin, nmax, number)
    start = time()
    match=[]
    for n in range(nmin, nmax):
        if number in str(n):
            match.append(n)
    end = time() - start
    print("found {0} in {1:.4f}sec".format(len(match),end))
    return match



def _concurrent_map(nmax, number, workers):
    '''Function that utilises concurrent.futures.ProcessPoolExecutor.map to
       find the occurrences of a given number in a number range in a parallelised
       manner.'''
    # 1. Local variables
    chunk = nmax // workers
    futures = []
    found =[]
    #2. Parallelization
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        # 2.1. Discretise workload and submit to worker pool
        for i in range(workers):
            cstart = chunk * i
            cstop = chunk * (i + 1) if i != workers - 1 else nmax
            futures.append(executor.submit(_findmatch, cstart, cstop, number))
    return chain.from_iterable(f.result() for f in cf.as_completed(futures))


nmax = int(1E8) # Number range maximum.
number = str(5) # Number to be found in number range.
workers = 4     # Pool of workers

start = time()
a = _concurrent_map(nmax, number, workers)
end = time() - start
print('\n main')
print('workers = ', workers)
print("found {0} in {1:.4f}sec".format(sum(1 for x in a),end))



guidemaker.