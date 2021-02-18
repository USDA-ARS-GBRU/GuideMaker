"""Core classes and functions for GuideMaker

"""
import os
import yaml
from typing import List, Dict, Tuple
import logging
from itertools import product
import gzip
import hashlib
import statistics
from collections import Counter
import numpy as np
#from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
import nmslib
from pybedtools import BedTool
import pandas as pd
import re
import gc
from Bio import Seq
from collections import deque


logger = logging.getLogger('guidemaker.core')



def load_parameters(yaml_file):
    """load yaml file with global variables and hyper paramteres
    """
    with open(yaml_file, "r") as f:
        y = yaml.load(stream=f, Loader=yaml.FullLoader)
        return y


# load Global variable from paramters.yaml
yaml_dict = load_parameters("parameters.yaml")

def is_gzip(filename):
    try:
        with open(filename, "rb") as f:
            logging.info("check if %s is gzipped" % filename)
            return f.read(2) == b'\x1f\x8b'
    except IOError as e:
        logging.error("Could not open the file %s to determine if it was gzipped" % filename)
        raise e
class Pam:
    """A Class representing a Protospacer Adjacent Motif (PAM)

    """

    def __init__(self, pam: str, pam_orientation: str) -> None:
        """Pam __init__

        Args:
            pam (str): A DNA string in ambiguous IUPAC format
            pam_orientation (str): [5prime | 3prime ]
                5prime means the order is 5'-[pam][target]-3'
                3prime means the order is 5'-[target][pam]-3'
        """
        for letter in pam.upper():
            assert letter in ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'X', 'N']
        assert pam_orientation in ["3prime", "5prime"]
        self.pam: str = pam.upper()
        self.pam_orientation: str = pam_orientation

    def __str__(self) -> str:
        return "A PAM object: {self.pam}".format(self=self)
    

    def find_targets(self, seq_record_iter: object, strand: str, target_len: int) -> List[object]:
        """Find all targets on a sequence that match for the PAM on the requested strand(s)

        Args:
            seq_record_iter (object): A Biopython SeqRecord iterator from SeqIO.parse
            strand (str): The strand to search choices: ["forward", "reverse", "both"]
            target_len (int): The length of the target sequence
        Returns:
            list: A list of Target class instances

        """
        target_list = []

        def reverse_complement(seq: str) -> object:
            """reverse complement of the PAM sequence
            """
            pamseq = Seq.Seq(seq)
            return str(pamseq.reverse_complement())

        def check_target(seq, target_len):
            if len(seq) == target_len and all(letters in ['A','T','C','G'] for letters in seq): # if not ATCG in the target then ignore those targets
                return True
            return False


        #https://stackoverflow.com/questions/18933711/find-all-occurrences-of-a-substring-including-overlap
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
            

        # run_for_5p
        def run_for_5p(expandpam, seq):
            for eachpam in expandpam:
                exact_pam = eachpam
                matches_g = findall(exact_pam, seq)
                for i in matches_g:
                    start = i + len(self.pam)
                    stop = i + len(self.pam) + target_len
                    target_seq = seq[start:stop]
                    if check_target(target_seq, target_len):
                        target_list.append(Target(seq=target_seq,
                                                        exact_pam=exact_pam,
                                                        # forward = 0, reverse =1
                                                        strand=0,
                                                        # 5prime =True, 3prime = False
                                                        pam_orientation=True,
                                                        seqid=id,
                                                        start = start,
                                                        stop=stop))
            #print("f5p: " + str(len(target_list)),getsize(target_list)/1e+9 , "GB")


        def run_for_3p(expandpam, seq):
            for eachpam in expandpam:
                exact_pam = eachpam
                matches_g = findall(exact_pam, seq)
                for i in matches_g:
                    start = i - target_len
                    stop = i 
                    target_seq = seq[start:stop]
                    if check_target(target_seq, target_len):
                        target_list.append(Target(seq=target_seq,
                                                        exact_pam=exact_pam,
                                                        # forward = 0, reverse =1
                                                        strand=0,
                                                        # 5prime =True, 3prime = False
                                                        pam_orientation=False,
                                                        seqid=id,
                                                        start = start,
                                                        stop=stop))
            #print("f3p: " + str(len(target_list)),getsize(target_list)/1e+9 , "GB")


        def run_rev_5p(expandpam, seq):
            for eachpam in expandpam:
                exact_pam = eachpam
                matches_g = findall(exact_pam, seq)
                for i in matches_g:
                    start = i - target_len
                    stop = i
                    target_seq = reverse_complement(seq[start:stop])
                    if check_target(target_seq, target_len):
                        target_list.append(Target(seq=target_seq,
                                                        exact_pam=reverse_complement(exact_pam),
                                                        # forward = 0, reverse =1
                                                        strand=1,
                                                        # 5prime =True, 3prime = False
                                                        pam_orientation=True,
                                                        seqid=id,
                                                        start = start,
                                                        stop=stop))
            #print("r5p: " + str(len(target_list)),getsize(target_list)/1e+9 , "GB")

        def run_rev_3p(expandpam, seq):
            for eachpam in expandpam:
                exact_pam = eachpam
                matches_g = findall(exact_pam, seq)
                for i in matches_g:
                    start = i + len(self.pam)
                    stop = i + len(self.pam) + target_len
                    target_seq = reverse_complement(seq[start : stop])
                    if check_target(target_seq, target_len):
                        target_list.append(Target(seq=target_seq,
                                                        exact_pam=reverse_complement(exact_pam),
                                                        # forward = 0, reverse =1
                                                        strand=1,
                                                        # 5prime =True, 3prime = False
                                                        pam_orientation=False,
                                                        seqid=id,
                                                        start = start,
                                                        stop=stop))
            #print("r3p: " + str(len(target_list)), getsize(target_list)/1e+9 ,"GB")


        for record in seq_record_iter:
                id = record.id
                seq =str(record.seq)
                if self.pam_orientation == "5prime":
                    run_for_5p(extend_ambiguous_dna(self.pam), seq)
                    run_rev_5p(extend_ambiguous_dna(reverse_complement(self.pam)), seq)
                elif self.pam_orientation == "3prime":
                    run_for_3p(extend_ambiguous_dna(self.pam), seq)
                    run_rev_3p(extend_ambiguous_dna(reverse_complement(self.pam)), seq)
        gc.collect() # clear memory after each chromosome
        return target_list


class Target:
    """A class representing a candidate target sequence for a PAM

    This is an object for holding data on possible target sequences
    adjacent to PAM sites.
    """
    __slots__ = ['seq', 'exact_pam','strand','pam_orientation','seqid','start','stop','md5']
    def __init__(self, seq: str, exact_pam: str, strand: str, pam_orientation: str,
                 seqid: str, start: int, stop: int) -> object:
        self.seq: str = seq
        self.exact_pam: str = exact_pam
        self.strand: str = strand
        self.pam_orientation: str = pam_orientation
        self.seqid: str = seqid
        self.start: int = start
        self.stop: int = stop
        self.md5: str = hashlib.md5(seq.encode()).hexdigest()
    def __str__(self) -> str:
        return "A Target object: {self.seq} on sequence {self.seqid} \
                position {self.start}".format(self=self)
    def __eq__(self, other):
        return other == self.seq
    def __ne__(self, other):
        return not self.__eq__(other)
    def __len__(self):
        return len(self.seq)


class TargetList:
    """A Class representing a set of guide RNA targets

    The class includes all targets in a set, a dict filtered to be unique in
    the region near the PAM, and a dict with edit distances for sequences.

    """
    def __init__(self, targets: List, lu: int, hammingdist: int=2, knum: int=2) -> None:
        self.targets: List = targets
        self.lu: int = lu # lenght of unique zone
        self.hammingdist: int = hammingdist
        self.knum: int = knum
        self.unique_targets: dict = {}
        self.nmslib_index: object = None
        self.neighbors: dict = {}
        self.ncontrolsearched: int = None
        self.gc_percent: float = None
        self.genomesize: float = None

    def __str__(self):
        info = "TargetList: contains a set of {} potential PAM targets".format(len(self.targets))
        return info

    def __len__(self):
        return len(self.targets)
    
    
    def check_restriction_enzymes(self, restriction_enzyme_list: list=[]):
        """Check for restriction enzymes
        
        Returns:
            None (but restriction enzyme checked tagets to self and targets are updated)
        """
        element_to_exclude = []
        for record in set(restriction_enzyme_list):
            element_to_exclude.append(extend_ambiguous_dna(record))
        element_to_exclude = sum(element_to_exclude, []) # flatout list of list to list
        # retrive targets if doesnot contain the strings specify in the restriction enzyme list
        #self.targets = [x for x in  self.targets if not any(restenzyme in x.seq for restenzyme in element_to_exclude)]
        to_delete = []
        for tobj in self.targets:
            for restenzyme in element_to_exclude:
                if restenzyme in tobj.seq:
                    to_delete.append(tobj)
        # update self.targets- such that target object if occures in to_delete list, is not selected
        self.targets = [tobj for tobj in self.targets if tobj not in to_delete]

    def _one_hot_encode(self, seq_list: List[object])-> List[str]:
        """One hot encode Target DNA as a binary string representation for LMSLIB

        """
        charmap = {'A': '1 0 0 0', 'C': '0 1 0 0', 'G': '0 0 1 0', 'T': '0 0 0 1'}
        def seq_to_bin(target_obj):
            seq = target_obj.seq
            charlist = [charmap[letter] for letter in seq]
            return " ".join(charlist)
        return list(map(seq_to_bin, seq_list))


    def find_unique_near_pam(self) -> None:
        """identify unique sequences in the target list

        The function filters a list of Target objects for targets that
        are unique in the region closest to the PAM. The region length is defined
        by the lu (length of unique zone).

        Args:
            lu (int): Length of conserved sequence close to PAM
        """
        def _get_prox(target):
            if target.pam_orientation == True: # 5prime = Truem 3prime=False
            	if self.lu == 0:
            		return target.seq
            	else:
            		return target.seq[0:self.lu]
            elif target.pam_orientation == False:# 5prime = Truem 3prime=False
            	if self.lu == 0:
            		return target.seq
            	else:
                	return target.seq[(len(target) - self.lu):]
        lu_dict ={}
        for target in self.targets:
            proximal = _get_prox(target)
            if proximal in lu_dict.keys():
                lu_dict[proximal].append(target)
            else:
                lu_dict[proximal] = [target]
        filteredlist = deque()
        for lkey, lval in lu_dict.items():
            if len(lval) == 1:
                filteredlist.append(lval[0])
        self.unique_targets = list(filteredlist)

    def _make_full_unique_targets(self):
        full_targerts= []
        for unq_target in self.unique_targets:
            full_targerts.append(str(unq_target.seq))
        full_unique = set(full_targerts)
        return full_unique
        

    def create_index(self, M: int=yaml_dict['NMSLIB']['M'], num_threads=2, efC: int=yaml_dict['NMSLIB']['efc'], post: int=yaml_dict['NMSLIB']['post']) -> None:
        """Create nmslib index

        Converts converts self.targets to binary one hot encoding and returns. NMSLIB index in

        Args:
            num_threads (int): cpu threads
            M (int): Controls the number of bi-directional links created for each element
                during index construction. Higher values lead to better results at the expense
                of memory consumption. Typical values are 2 -100, but for most datasets a
                range of 12 -48 is suitable. Can’t be smaller than 2.
            efC (int): Size of the dynamic list used during construction. A larger value means
                   a better quality index, but increases build time. Should be an integer value
                   between 1 and the size of the dataset.

        Returns:
            None (but writes NMSLIB index to self)
        """


        logging.info("unique targets for index: %s" % len(self.targets))
        bintargets = self._one_hot_encode(self.targets)
        index_params = {'M': M, 'indexThreadQty': num_threads,'efConstruction': efC, 'post': post}
        index = nmslib.init(space='bit_hamming',
                            dtype=nmslib.DistType.INT,
                            data_type=nmslib.DataType.OBJECT_AS_STRING,
                            method='hnsw')
        index.addDataPointBatch(bintargets)
        index.createIndex(index_params, print_progress=True)
        self.nmslib_index = index

    def get_neighbors(self, ef: int=yaml_dict['NMSLIB']['ef'], num_threads=2) -> None:
        """Get nearest neighbors for sequences removing sequences that
         have neighbors less than the Hamming distance threshold

        For the list of all targets calculate the (knum) nearest neighbors.
        filter out targets with close neighbors and
        Writes a dictionary to self.neighbors:
        self.neighbors[seq]{target: seq_obj, neighbors: {seqs:[s1, s1, ...], dist:[d1, d1,...]}}

        Args: None
        Returns: None
        """
        if not self.unique_targets:
            self.find_unique_near_pam()
        unique_bintargets = self._one_hot_encode(self.unique_targets)
        self.nmslib_index.setQueryTimeParams({'efSearch': ef})
        results_list = self.nmslib_index.knnQueryBatch(unique_bintargets,
                                               k=self.knum, num_threads = num_threads)
        neighbor_dict = {}
        for i, entry in enumerate(results_list):

            queryseq = self.unique_targets[i].seq
            hitseqidx = list(entry[0])
            hammingdist = list(entry[1])
            if hammingdist[1] >= 2 * self.hammingdist: # multiply by 4 b/c each base is one hot encoded in 4 bits
                neighbors = {"seqs": [self.targets[x].seq for x in hitseqidx], # reverse this?
                             "dist": [int(x/2) for x in hammingdist]}
                neighbor_dict[queryseq] = {"target": self.unique_targets[i],
                                           "neighbors": neighbors}
        self.neighbors = neighbor_dict

    def export_bed(self) -> object:
        """export the targets in self.neighbors to a bed format file

        Args:
            file (str): the name and location of file to export

        Returns:
            (obj): A Pandas Dataframe in Bed format
        """
        bdict = dict(chrom = [], chromstart = [], chromend = [], name = [], score = [], strand = [])
        for rec in self.neighbors.values():
            bdict['chrom'].append(rec["target"].seqid)
            bdict['chromstart'].append(rec["target"].start)
            bdict['chromend'].append(rec["target"].stop)
            bdict['name'].append(rec["target"].seq)
            bdict['score'].append(0)
            if rec["target"].strand == 0:
                bdict['strand'].append("+")
            elif rec["target"].strand == 1:
                bdict['strand'].append("-")
        df = pd.DataFrame.from_dict(bdict)
        df.sort_values(by=['chrom', 'chromstart'], inplace=True)
        return df

    def get_control_seqs(self, seq_record_iter: object, length: int=20, n: int=1000,
                         num_threads: int=2) -> Tuple[float, float, object]:
        """Create random sequences with a specified GC probability and find seqs with the greatest
         distance to any sequence flanking a PAM site

        Args:
            seq_record_iter (Bio.SeqIO): an iterator of Fastas
            length (int): length of the sequence, must match the index
            n = number of sequences to  return
            num_threads (int) nuer of processor threads
        """
        # get GC percent
        totlen = 0
        gccnt = 0
        for record in seq_record_iter:
            gccnt += GC(record.seq) * len(record)
            totlen += len(record)
        gc = gccnt/(totlen*100)
        #print("Percentage of GC content in the input genome: "+"{:.2f}".format(gc * 100))
        self.gc_percent = gc * 100
        self.genomesize = totlen / (1024 * 1024)
        
        minimum_hmdist=0
        sm_count = 0
        search_mult = 0
        
        #  search_mult (int): search this times n sequences
        search_multiple = yaml_dict['CONTROL_SEARCH_MULTIPLE']
       
        try:
            while  minimum_hmdist < 7 or search_mult ==  10000:
                # generate random sequences
                seqs = []
                search_mult = search_multiple[sm_count]
                for i in range(n  * search_mult):
                    seqs.append("".join(np.random.choice(a=["G", "C", "A", "T"], size=length,
                                                         replace=True, p=[gc/2, gc/2, (1 - gc)/2, (1 - gc)/2])))
                # one hot encode sequences
                binseq = []
                charmap = {'A': '1 0 0 0', 'C': '0 1 0 0', 'G': '0 0 1 0', 'T': '0 0 0 1'}
                for seq in seqs:
                    charlist = [charmap[letter] for letter in seq]
                    binseq.append(" ".join(charlist))
                rand_seqs = self.nmslib_index.knnQueryBatch(binseq, k=2, num_threads=num_threads)
                distlist = []
                for i in rand_seqs:
                    distlist.append(i[1][0])
                zipped = list(zip(seqs, distlist))
                dist_seqs = sorted(zipped, reverse=True, key=lambda x: x[1])
                sort_seq = [item[0] for item in dist_seqs][0:n]
                sort_dist = [item[1]/2 for item in dist_seqs][0:n]
                minimum_hmdist = int(min(sort_dist))
                sm_count += 1
        except IndexError as e:
           # print("Number of random control searched: ", search_mult * n)
            pass
        
        total_ncontrolsearched = search_mult * n
        self.ncontrolsearched = total_ncontrolsearched
        randomdf = pd.DataFrame(data={"Sequences":sort_seq, "Hamming distance":sort_dist})
        def create_name(seq):
            return "Cont-" + hashlib.md5(seq.encode()).hexdigest()
        randomdf['name'] = randomdf["Sequences"].apply(create_name)
        randomdf = randomdf[["name", "Sequences", "Hamming distance"]]
        return (min(sort_dist),
                statistics.median(sort_dist),
                randomdf)

class Annotation:
    def __init__(self, genbank_list: List[str], target_bed_df: object) -> None:
        """Annotation class for data and methods on targets and gene annotations

        Args:
            genbank_list (List[str]): A list of genbank files from a single genome
            target_bed_df (object): A pandas dataframe in Bed format with the
                locations of targets in the genome
        Returns:
            None
        """
        self.genbank_list: List[str] = genbank_list
        self.target_bed_df: object = target_bed_df
        self.genbank_bed_df: object = None
        self.feature_dict: Dict = None
        self.nearby: object = None
        self.filtered_df: object = None
        self.qualifiers: object = None
        
    def _get_genbank_features(self, feature_types: List[str] = None) -> None:
        """Parse genbank records into pandas DF/Bed format and dict format saving to self

        Args:
            feature_types (List[str]): a list of Genbank feature types to use

        Returns:
            None
        """
        if feature_types is None:
            feature_types = ["CDS"]
        feature_dict = {}
        pddict = dict(chrom=[], chromStart=[], chromEnd=[], name=[], score=[], strand=[])
        for gbfile in self.genbank_list:
            try:
                if is_gzip(gbfile):
                    f = gzip.open(gbfile, mode='rt')
                else:
                    f = open(gbfile, mode='r')
            except IOError as e:
                logging.error("The genbank file %s could not be opened" % gbfile)
                raise e
            genbank_file = SeqIO.parse(f, "genbank")
            for entry in genbank_file:
                for record in entry.features:
                    if record.type in feature_types:
                        if record.strand in [1, -1, "+","-"]:
                            pddict["strand"].append("-" if record.strand < 0 else "+")
                            featid = hashlib.md5(str(record).encode()).hexdigest()
                            pddict['chrom'].append(entry.id)
                            pddict["chromStart"].append(record.location.start.position)
                            pddict["chromEnd"].append(record.location.end.position)
                            pddict["name"].append(featid)
                            pddict["score"].append(0)
                            for qualifier_key, qualifier_val in record.qualifiers.items():
                                if not qualifier_key in feature_dict:
                                    feature_dict[qualifier_key] = {}
                                feature_dict[qualifier_key][featid] = qualifier_val
            genbankbed = pd.DataFrame.from_dict(pddict)
            self.genbank_bed_df = genbankbed
            self.feature_dict = feature_dict
        f.close()

    def _get_qualifiers(self, min_prop: float = yaml_dict['MINIMUM_PROPORTION'], excluded: List[str] = None) -> object:
        """Create a dataframe with features and their qualifier values

        Create a dataframe with features and their qualifier values for
        all qualifiers over the minimum threshold (except 'translation'). Add
        to self.qualifiers

        Args:
            min_prop (float): A float between 0-1 representing the fraction of
            features the qualifier must be present in to be included in the dataframe
            excluded (List(str)): A list of genbank qualifiers to exclude, Default ["translation"]
        Returns:
            None

        """
        if excluded is None:
            excluded = ["translation"]
        final_quals = []
        qual_df = pd.DataFrame(data ={"Feature id":[]})
        for quals in self.feature_dict:
            if len(quals)/len(self.feature_dict) > min_prop:
                final_quals.append(quals)
        for qualifier in final_quals:
            if qualifier not in excluded:
                featlist = []
                quallist = []
                for feat, qual in self.feature_dict[qualifier].items():
                    featlist.append(feat)
                    quallist.append(";".join([str(i) for i in qual]))
                tempdf = pd.DataFrame({'Feature id': featlist, qualifier: quallist})
                qual_df = qual_df.merge(tempdf, how="outer", on="Feature id")
        self.qualifiers = qual_df

    def _get_nearby_features(self) -> None:
        """Adds downstream information to the given target sequences and mapping information

        Args:
            None

        Returns:
            None

        Note:
            writes a dataframe of nearby features to self.nearby
        """
        # Import Features and sort by chromosome and then by start position in ascending order
        featurebed = BedTool.from_dataframe(self.genbank_bed_df)
        featurebed = featurebed.sort()
        # import guide files and sort by chromosome and then by start position in ascending order
        mapbed = BedTool.from_dataframe(self.target_bed_df)
        mapbed = mapbed.sort()
        # get feature downstream of target sequence
        downstream = mapbed.closest(featurebed, d=True, fd=True, D="a", t="first")
        # get feature upstream of target sequence
        upstream = mapbed.closest(featurebed, d=True, id=True, D="a", t="first")
        headers = {0: "Accession", 1: "Guide start", 2: "Guide end", 3:"Guide sequence",
                   4: "Score", 5:"Guide strand", 6: "Feature Accession",
                   7: "Feature start", 8:"Feature end", 9:"Feature id",
                   10:"Feature score", 11:"Feature strand", 12: "Feature distance"}
        downstream: pd.DataFrame = downstream.to_dataframe(disable_auto_names=True, header=None)
        downstream['direction'] = 'downstream'
        upstream = upstream.to_dataframe(disable_auto_names=True, header=None)
        upstream['direction'] = 'upstream'
        upstream = upstream.append(downstream)
        self.nearby = upstream.rename(columns=headers)

    def _filter_features(self, before_feat: int=100, after_feat: int=200) -> None:
        """merge targets with Feature list and filter for guides close enough to interact.

        Args:
            before_feat (int): The maximum distance before the start of a feature measured from closest point to guide
            after_feat (int): The maximum distance after the start codon (into the gene)

        Returns:
            None
        """
        # for guides in the same orientation as the targets ( +/+ or -/-) select guides that are within
        #  before_feat of the gene start
        filtered_df = self.nearby.query('`Guide strand` == `Feature strand` and 0 < `Feature distance` < @before_feat')
        # for guides in the +/+ orientation select guides where the end is within [before_feat] of the gene start
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "+" and `Feature strand` == "+" \
                                             and `Feature distance` == 0 and \
                                             `Guide end` - `Feature start` < @after_feat'))
        # for guides in the -/- orientation select guides where the end is within [before_feat] of the gene start
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "-" and `Feature strand` == "-" \
                                                     and `Feature distance` == 0 \
                                                     and `Feature end` - `Guide start` < @after_feat'))
        # Select guides where target is + and guide is - and the guide is infront of the gene
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "-" and `Feature strand` == "+" and \
                                                     0 <`Feature start` - `Guide end` < @before_feat' ))
        # Select guides where target is - and guide is + and the guide is infront of the gene
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "+" and `Feature strand` == "-" and \
                                                     0 <`Guide start` - `Feature end` < @before_feat' ))
        # Select guides where target is + and guide is - and the guide is is within [before_feat] of the gene start
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "-" and `Feature strand` == "+" and \
                                                             0 <`Guide end` -`Feature start`  < @after_feat'))
        # Select guides where target is - and guide is + and the guide is is within [before_feat] of the gene start
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "+" and `Feature strand` == "-" and \
                                                             0 <`Feature end` - `Guide start` < @after_feat'))

        self.filtered_df = filtered_df

    def _format_guide_table(self, targetlist: List[object]) -> object :
        """Create guide table for output

        """
        def gc(seq):
            cnt = 0
            for letter in seq:
                if letter in ["G", "C"]:
                    cnt += 1
            return cnt/len(seq)
        def get_exact_pam(seq):
            return targetlist.neighbors[seq]["target"].exact_pam
        def get_guide_hash(seq):
            return targetlist.neighbors[seq]["target"].md5
        def get_off_target_score(seq):
            dlist = targetlist.neighbors[seq]["neighbors"]["dist"]
            s = [str(i) for i in dlist]
            return ";".join(s)
        def get_off_target_seqs(seq):
            slist = targetlist.neighbors[seq]["neighbors"]["seqs"]
            return ";".join(slist)
        pretty_df = self.filtered_df.copy()
        pretty_df['GC'] = pretty_df['Guide sequence'].apply(gc)
        pretty_df['PAM'] = pretty_df['Guide sequence'].apply(get_exact_pam)
        pretty_df['Guide name'] = pretty_df['Guide sequence'].apply(get_guide_hash)
        pretty_df['Target strand'] = np.where(pretty_df['Guide strand'] == pretty_df['Feature strand'], 'coding', 'non-coding')
        pretty_df['Similar guide distances'] = pretty_df['Guide sequence'].apply(get_off_target_score)
        pretty_df['Similar guides'] = pretty_df['Guide sequence'].apply(get_off_target_seqs)
        pretty_df = pretty_df[['Guide name', "Guide sequence", 'GC', "Accession","Guide start", "Guide end",
                    "Guide strand", 'PAM',  "Feature id",
                    "Feature start", "Feature end", "Feature strand",
                    "Feature distance", 'Similar guides', 'Similar guide distances']]
        pretty_df: object = pretty_df.merge(self.qualifiers, how="left", on="Feature id")
        pretty_df = pretty_df.sort_values(by=['Accession', 'Feature start'])
        pretty_df['Guide start'] = pretty_df['Guide start'] + 1
        pretty_df['Feature start'] = pretty_df['Feature start'] + 1
        return pretty_df
    
    def locuslen(self):
        locus_count = len(self.feature_dict['locus_tag' or 'locus'].keys())
        return(locus_count)

        
def get_fastas(filelist, tempdir=None):
    """Saves a Fasta and from 1 or more Genbank files (may be gzipped)
    Args:
        filelist (str): Genbank file to process

    Returns:
        None
    """
    try:
        fastpath = os.path.join(tempdir, "forward.fasta")
        with open(fastpath, "w") as f1:
            for file in filelist:
                if is_gzip(file):
                    with gzip.open(file, 'rt') as f:
                        records = SeqIO.parse(f, "genbank")
                        SeqIO.write(records, f1, "fasta")
                else:
                    with open(file, 'r') as f:
                        records = (SeqIO.parse(f, "genbank"))
                        SeqIO.write(records, f1, "fasta")
        return fastpath
    except Exception as e:
        print("An error occurred in input genbank file %s" % file)
        raise e

def extend_ambiguous_dna(seq):
            """return list of all possible sequences given an ambiguous DNA input
            """
            dna_dict = Seq.IUPAC.IUPACData.ambiguous_dna_values
            extend_list = []
            for i in product(*[dna_dict[j] for j in seq]):
                extend_list.append("".join(i))
            return extend_list