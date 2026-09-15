"""Core classes and functions for GuideMaker."""
import os
import re
import yaml
import logging
import gzip
import hashlib
import statistics
import nmslib
import regex
import gc
from typing import List, Dict, TypeVar, Generator
from itertools import product
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from pybedtools import BedTool
from Bio import Seq
from copy import deepcopy
import pandas as pd
import numpy as np
import altair as alt
from guidemaker import doench_predict
from guidemaker import cfd_score_calculator

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None

def is_gzip(filename: str):
    try:
        with open(filename, "rb") as f:
            logger.info("check if %s is gzipped" % filename)
            return f.read(2) == b'\x1f\x8b'
    except IOError as e:
        logger.error("Could not open the file %s to determine if it was gzipped" % filename)
        raise e


class PamTarget:

    """
    A Class representing a Protospacer Adjacent Motif (PAM) and targets.

    The classincludes all targets for given PAM as a dataframe,PAM and target attributes,
    and methods to find target and control sequences.

    """

    def __init__(self, pam: str, pam_orientation: str, dtype: str) -> None:
        """
        Pam __init__

        Args:
            pam (str): A DNA string in ambiguous IUPAC format
            pam_orientation (str): [5prime | 3prime ]
                5prime means the order is 5'-[pam][target]-3'
                3prime means the order is 5'-[target][pam]-3'
            dtype (str): hamming or leven

        Returns:
            None
        """
        for letter in pam.upper():
            assert letter in ['A', 'C', 'G', 'T', 'M', 'R', 'W',
                'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'X', 'N']
        assert pam_orientation in ["3prime", "5prime"]
        self.pam: str = pam.upper()
        self.pam_orientation: str = pam_orientation
        self.dtype: str = dtype

    def __str__(self) -> str:
        """
        str __init__

        Args:
            self

        Returns:
            self(str)
        """
        return "A PAM object: {self.pam}".format(self=self)

    def find_targets(self, seq_record_iter: object, target_len: int) -> pd.DataFrame:
        """
        Find all targets on a sequence that match for the PAM on both strand(s)

        Args:
            seq_record_iter (object): A Biopython SeqRecord iterator from SeqIO.parse
            target_len (int): The length of the target sequence

        Returns:
            PandasDataFrame: A pandas dataframe with of matching targets
        """

        def reverse_complement(seq: str) -> str:
            """
            Reverse complement of the PAM sequence

            Args:
                seq (str): A DNA string

            Returns:
                str: A reverse complement of DNA string
            """
            bpseq = Seq.Seq(seq)
            return str(bpseq.reverse_complement())

        def pam2re(pam: str) -> str:
            """
            Convert an IUPAC ambiguous PAM to a Regex expression

            Args:
                pam (str): A DNA string

            Returns:
                str: A Regex expression
            """
            dnaval = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
                      'M': '[A|C]', 'R': '[A|G]', 'W': '[A|T]', 'S': '[C|G]',
                      'Y': '[C|T]', 'K': '[G|T]', 'V': '[A|C|G]', 'H': '[A|C|T]',
                      'D': '[A|G|T]', 'B': '[C|G|T]', 'X': '[G|A|T|C]', 'N': '[G|A|T|C]'}
            return "".join([dnaval[base] for base in pam])

        #                5prime means the order is 5'-[pam][target]-3'
        #                3prime means the order is 5'-[target][pam]-3'

        def check_target(seq: str, target_len: int) -> bool:
            """
            Check targets for guidelength and DNA bases

            Args:
                seq (str): A DNA string
                target_len(int): Guide length

            Returns:
                bool: True or False
            """
            if len(seq) == target_len and all(letters in ['A', 'T', 'C', 'G'] for letters in seq):  # if not ATCG in the target then ignore those targets
                return True
            return False

        def run_for_5p(pam_pattern: str, dnaseq: str, target_len: int) -> Generator:
            """
            Search for guides with 5prime pam orientation in the forward strand

            Args:
                pam_pattern (str): A DNA string representing PAM
                dnaseq (str): A DNA string representing genome
                target_len (int): Guide length

            Returns:
                (Generator): A generator with target_seq, exact_pam, start, stop, strand, and pam_orientation
            """
            for match_obj in regex.finditer(pattern=pam_pattern, string=dnaseq, overlapped=True):
                target_seq = dnaseq[match_obj.end(): match_obj.end() + target_len]
                target_seq30 = dnaseq[match_obj.start()-3: match_obj.start()+27]
                ## 5'-[guide of 25 nt][exact pam, 3nt][next two]-3'
                if check_target(target_seq, target_len):
                    exact_pam = match_obj.group(0)
                    start = match_obj.end()
                    stop = match_obj.end() + target_len
                    # 5prime =True, 3prime = False
                    pam_orientation = True
                    # forward =True, reverse = False
                    strand = True
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation, target_seq30



        def run_for_3p(pam_pattern, dnaseq, target_len) -> Generator:
            """
            Search for guides with 3prime pam orientation in the reverse strand

            Args:
                pam_pattern (str): A DNA string representing PAM
                dnaseq (str): A DNA string representing genome
                target_len (int): Guide length

            Returns:
                (Generator): A generator with target_seq, exact_pam, start, stop, strand, and pam_orientation
            """
            for match_obj in regex.finditer(pattern=pam_pattern, string=dnaseq, overlapped=True):
                target_seq = dnaseq[match_obj.start() - target_len: match_obj.start()]
                target_seq30 = dnaseq[match_obj.end()-27 :match_obj.end()+3]
                if check_target(target_seq, target_len):
                    exact_pam = match_obj.group(0)
                    start = match_obj.start() - target_len
                    stop = match_obj.start()
                    # 5prime =True, 3prime = False
                    pam_orientation = False
                    # forward =True, reverse = False
                    strand = True
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation, target_seq30

        def run_rev_5p(pam_pattern, dnaseq, target_len) -> Generator:
            """
            Search for guides with 5prime pam orientation in the reverse strand

            Args:
                pam_pattern (str): A DNA string representing PAM
                dnaseq (str): A DNA string representing genome
                target_len (int): Guide length

            Returns:
                (Generator): A generator with target_seq, exact_pam, start, stop, strand, and pam_orientation
            """
            for match_obj in regex.finditer(pattern=pam_pattern, string=dnaseq, overlapped=True):
                target_seq = reverse_complement(
                    dnaseq[match_obj.start() - target_len: match_obj.start()])
                target_seq30 = reverse_complement(
                    dnaseq[match_obj.end()-27:match_obj.end()+3])
                if check_target(target_seq, target_len):
                    exact_pam = reverse_complement(match_obj.group(0))
                    start = match_obj.start() - target_len
                    stop = match_obj.start()
                    # 5prime =True, 3prime = False
                    pam_orientation = True
                    # forward =True, reverse = False
                    strand = False
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation, target_seq30

        def run_rev_3p(pam_pattern, dnaseq, target_len) -> Generator:
            """
            Search for guides with 3prime pam orientation in the reverse strand

            Args:
                pam_pattern (str): A DNA string representing PAM
                dnaseq (str): A DNA string representing genome
                target_len (int): Guide length

            Returns:
                (Generator): A generator with target_seq, exact_pam, start, stop, strand, and pam_orientation
            """
            for match_obj in regex.finditer(pattern=pam_pattern, string=dnaseq, overlapped=True):
                target_seq = reverse_complement(
                    dnaseq[match_obj.end(): match_obj.end() + target_len])
                target_seq30 = reverse_complement(dnaseq[match_obj.start()-3:match_obj.start()+27])
                if check_target(target_seq, target_len):
                    exact_pam = reverse_complement(match_obj.group(0))
                    start = match_obj.end()
                    stop = match_obj.end() + target_len
                    # 5prime =True, 3prime = False
                    pam_orientation = False
                    # forward =True, reverse = False
                    strand = False
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation, target_seq30

        target_list = []
        for record in seq_record_iter:
            record_id = record.id
            seq = str(record.seq)
            if self.pam_orientation == "5prime":
                # forward
                for5p = pd.DataFrame(run_for_5p(pam2re(self.pam), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation", "target_seq30"])
                for5p["seqid"] = record_id
                # string to boolean conversion is not straight - as all string were set to Trues- so change the encoding in functions above.
                # https://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python/715455#715455
                for5p = for5p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(for5p)
                # reverse
                rev5p = pd.DataFrame(run_rev_5p(pam2re(reverse_complement(self.pam)), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation","target_seq30"])
                rev5p["seqid"] = record_id
                rev5p = rev5p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(rev5p)
                # Question? Append directly vs. concat then append? https://ravinpoudel.github.io/AppendVsConcat/
            elif self.pam_orientation == "3prime":
                # forward
                for3p = pd.DataFrame(run_for_3p(pam2re(self.pam), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation","target_seq30"])
                for3p["seqid"] = record_id
                for3p = for3p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(for3p)
                # reverse
                rev3p = pd.DataFrame(run_rev_3p(pam2re(reverse_complement(self.pam)), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation","target_seq30"])
                rev3p["seqid"] = record_id
                rev3p = rev3p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(rev3p)
            gc.collect()  # clear memory after each chromosome
        target_list = [item for item in target_list if not item.empty]
        df_targets = pd.concat(target_list, ignore_index=True)
        df_targets = df_targets.assign(seedseq=np.nan, hasrestrictionsite=np.nan, isseedduplicated=np.nan)
        df_targets = df_targets.astype({"seedseq": 'str', "isseedduplicated": 'bool'})
        df_targets = df_targets.assign(dtype=self.dtype)
        df_targets = df_targets.astype({"dtype": 'category'})
        return df_targets


class TargetProcessor:

    """
    A Class representing a set of guide RNA targets.

    The class includes all targets in a dataframe, methods to process target and a dict with edit distances for sequences.

    """

    def __init__(self, targets: pd.DataFrame, lsr: int, editdist: int = 2, knum: int = 2) -> None:
        """
        TargetProcessor __init__

        Args:
            targets (PandasDataFrame): Dataframe with output from class PamTarget
            lsr (int): Length of seed region
            editdist (int): Edit distance
            knum (int): Number of negative controls

        Returns:
            None
        """
        self.targets = targets  # pandas dataframe
        self.lsr: int = lsr  # length of seed region
        self.editdist: int = editdist
        self.knum: int = knum
        self.nmslib_index: object = None
        self.neighbors: dict = {}
        self.closest_neighbor_df: pd.DataFrame = None
        self.ncontrolsearched: int = None
        self.gc_percent: float = None
        self.genomesize: float = None
        self.pam_orientation: bool = targets['pam_orientation'].iat[0]

    def __str__(self) -> None:
        """
        str __init__

        Args:
            self

        Return:
            None
        """
        info = "TargetList: contains a set of {} potential PAM targets".format(len(self.targets))
        return info

    def __len__(self) -> int:
        """
        len __init__ to display length of self.targets

        Args:
            self.targets

        Return:
            (int): Length of the self.targets
        """
        return len(self.targets)

    def check_restriction_enzymes(self, restriction_enzyme_list: list = []) -> None:
        """
        Check for restriction enzymes and its reverse complement within gRNA sequence

        Args:
            restriction_enzyme_list (list): A list with sequence for restriction enzymes

        Returns:
            None
        """
        element_to_exclude = []
        for record in set(restriction_enzyme_list):
            for letter in record.upper():
                assert letter in ['A', 'C', 'G', 'T', 'M', 'R', 'W',
                    'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'X', 'N']
            record_seq = Seq.Seq(record.upper())
            element_to_exclude.append(extend_ambiguous_dna(str(record_seq)))
            element_to_exclude.append(extend_ambiguous_dna(
                str(record_seq.reverse_complement())))  # reverse complement
        element_to_exclude = sum(element_to_exclude, [])  # flatout list of list to list with restriction enzyme sites
        if len(element_to_exclude) > 0:
            self.targets['hasrestrictionsite'] = self.targets['target'].str.contains('|'.join(element_to_exclude))
        else:
            self.targets['hasrestrictionsite'] = False

    def _one_hot_encode(self, seq_list: List[object]) -> List[str]:
        """One hot encode Target DNA as a binary string representation for NMSLIB."""
        charmap = {'A': '1 0 0 0', 'C': '0 1 0 0', 'G': '0 0 1 0', 'T': '0 0 0 1'}

        def seq_to_bin(seq):
            charlist = [charmap[letter] for letter in seq]
            return " ".join(charlist)
        return list(map(seq_to_bin, seq_list))

    def find_unique_near_pam(self) -> None:
        """
        Identify unique sequences in the target list

        The function filters a list of Target objects for targets that
        are unique in the region closest to the PAM. The region length is defined
        by the lsr (length of seed region that need to be unique).

        Args:
            lsr (int): Length of seed region that is close to PAM

        Returns:
            None
        """
        def _get_prox(tseq):  # get target sequence as input
            if self.pam_orientation == True:  # 5prime = True 3prime=False
                if self.lsr == 0:
                    return tseq
                else:
                    return tseq[0:self.lsr]
            elif self.pam_orientation == False:  # 5prime = True 3prime=False
                if self.lsr == 0:
                    return tseq
                else:
                    return tseq[(len(tseq) - self.lsr):]
        # https://stackoverflow.com/questions/12555323/adding-new-column-to-existing-dataframe-in-python-pandas
        self.targets = deepcopy(self.targets)
        self.targets.loc[:, 'seedseq'] = self.targets.loc[:, 'target'].apply(_get_prox)
        self.targets.loc[:, 'isseedduplicated'] = self.targets.loc[:, 'seedseq'].duplicated()

    def create_index(self, configpath: str, num_threads=2):
        """
        Create nmslib index

        Converts self.targets to binary one hot encoding and returns NMSLIB index

        Args:
            num_threads (int): cpu threads
            configpath (str): Path to config file which contains hyper parameters for NMSLIB

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
        with open(configpath) as cf:
            config = yaml.safe_load(cf)

        M, efC, post = config['NMSLIB']['M'], config['NMSLIB']['efc'], config['NMSLIB']['post']

        # index everything but not duplicates
        notduplicated_targets = list(set(self.targets['target'].tolist()))
        #mod_logger.info("unique targets for index: %s" % len(notduplicated_targets))
        if self.targets['dtype'].iat[0] == "hamming":
            bintargets = self._one_hot_encode(notduplicated_targets)
            index_params = {'M': M, 'indexThreadQty': num_threads, 'efConstruction': efC, 'post': post}
            index = nmslib.init(space='bit_hamming',
                            dtype=nmslib.DistType.INT,
                            data_type=nmslib.DataType.OBJECT_AS_STRING,
                            method='hnsw')
            index.addDataPointBatch(bintargets) # notduplicated_targets
            index.createIndex(index_params, print_progress=True)
            self.nmslib_index = index
        else:
            bintargets = notduplicated_targets
            index_params = {'M': M, 'indexThreadQty': num_threads, 'efConstruction': efC, 'post': post}
            index = nmslib.init(space='leven',
                            dtype=nmslib.DistType.INT,
                            data_type=nmslib.DataType.OBJECT_AS_STRING,
                            method='hnsw')
            index.addDataPointBatch(bintargets) # notduplicated_targets
            index.createIndex(index_params, print_progress=True)
            self.nmslib_index = index



    def get_neighbors(self, configpath, num_threads=2) -> None:
        """
        Get nearest neighbors for sequences removing sequences that
        have neighbors less than the Hamming distance threshold.
        For the list of all targets calculate the (knum) nearest neighbors.
        filter out targets with close neighbors and
        Writes a dictionary to self.neighbors:
        self.neighbors[seq]{target: seq_obj, neighbors: {seqs:[s1, s1, ...], dist:[d1, d1,...]}}

        Args:
            configpath (str): Path to a parameter config file
            num_threads (int): Number of threads

        Returns:
            None
        """
        with open(configpath) as cf:
            config = yaml.safe_load(cf)

        ef = config['NMSLIB']['ef']

        # unique_targets = self.targets.loc[self.targets['isseedduplicated']
        #     == False]['target'].tolist()
        # For indexing we need to use all targets -- for checking off-targets. For searching neighbours remove seed duplicated and one wiht restriction site.
        unique_targets = self.targets.loc[(self.targets['isseedduplicated']==False) | (self.targets['hasrestrictionsite']==False)]['target'].tolist()
        if self.targets['dtype'].iat[0] == "hamming":
            unique_bintargets = self._one_hot_encode(unique_targets)  # search unique seed one
        else:
            unique_bintargets = unique_targets

        self.nmslib_index.setQueryTimeParams({'efSearch': ef})
        results_list = self.nmslib_index.knnQueryBatch(unique_bintargets,
                                               k=self.knum, num_threads=num_threads)
        neighbor_dict = {}
        for i, entry in enumerate(results_list):
            queryseq = unique_targets[i]
            hitseqidx = entry[0].tolist()
            editdist = entry[1].tolist()
            if self.targets['dtype'].iat[0] == "hamming":
                # check that the closest sequence meets the min. dist. requirment. We multiply by 2 b/c each
                # base is one hot encoded. e.g. 1000 vs 0100 = 2 differences
                if editdist[1] >= 2 * self.editdist:
                    neighbors = {"seqs": [self.targets['target'].values[x] for x in hitseqidx],  # reverse this?
                                "dist": [int(x / 2) for x in editdist]}
                    neighbor_dict[queryseq] = {"target": unique_targets[i],
                                            "neighbors": neighbors}
            else:
               if editdist[1] >= self.editdist:
                    neighbors = {"seqs": [self.targets['target'].values[x] for x in hitseqidx],  # reverse this?
                                "dist": [int(x) for x in editdist]}
                    neighbor_dict[queryseq] = {"target": unique_targets[i],
                                            "neighbors": neighbors}
        self.neighbors = neighbor_dict

    def export_bed(self) -> object:
        """
        Export the targets in self.neighbors to a bed format file

        Args:
            file (str): the name and location of file to export

        Returns:
            (obj): A Pandas Dataframe in Bed format
        """
        # df = self.targets.copy()
        # why deepcopy - https://stackoverflow.com/questions/55745948/why-doesnt-deepcopy-of-a-pandas-dataframe-affect-memory-usage
        # select only guides that are not duplecated in the seedseq
        df = deepcopy(self.targets.loc[self.targets['isseedduplicated'] == False])
        df = df[["seqid", "start", "stop", "target", "strand"]]
        df.loc[:, 'strand'] = df.loc[:, 'strand'].apply(lambda x: '+' if x == True else '-')
        df.columns = ["chrom", "chromstart", "chromend", "name", "strand"]
        df.sort_values(by=['chrom', 'chromstart'], inplace=True)
        return df

    def get_control_seqs(self, seq_record_iter: object, configpath, length: int = 20, n: int = 10,
                         num_threads: int = 2) -> pd.DataFrame:
        """
        Create random sequences with a specified GC probability and find seqs with the greatest
        distance to any sequence flanking a PAM site

        Args:
            seq_record_iter (Bio.SeqIO): An iterator of fastas
            length (int): Length of the sequence, must match the index
            n (int): Number of sequences to  return
            num_threads (int): Number of processor threads

        Returns:
            (PandasDataFrame): A pandas dataframe with control sequence
        """

        with open(configpath) as cf:
            config = yaml.safe_load(cf)

        MINIMUM_HMDIST = config['CONTROL']['MINIMUM_HMDIST']

        MAX_CONTROL_SEARCH_MULTIPLE = max(config['CONTROL']['CONTROL_SEARCH_MULTIPLE'])

        #  search_mult (int): search this times n sequences
        CONTROL_SEARCH_MULTIPLE = config['CONTROL']['CONTROL_SEARCH_MULTIPLE']

        # get GC percent
        totlen = 0
        gccnt = 0
        for record in seq_record_iter:
            gccnt += gc_fraction(record.seq) * len(record)
            totlen += len(record)
        gc = gccnt / (totlen)
        self.gc_percent = gc * 100
        self.genomesize = totlen / (1024 * 1024)

        minimum_hmdist = 0
        sm_count = 0
        search_mult = 0

        try:
            while minimum_hmdist < MINIMUM_HMDIST or search_mult == MAX_CONTROL_SEARCH_MULTIPLE:
                # generate random sequences
                seqs = []
                search_mult = CONTROL_SEARCH_MULTIPLE[sm_count]
                for i in range(n * search_mult):
                    seqs.append("".join(np.random.choice(a=["G", "C", "A", "T"], size=length,
                                                         replace=True, p=[gc / 2, gc / 2, (1 - gc) / 2, (1 - gc) / 2])))
                # one hot encode sequences
                binseq = []
                charmap = {'A': '1 0 0 0', 'C': '0 1 0 0', 'G': '0 0 1 0', 'T': '0 0 0 1'}
                for seq in seqs:
                    if self.targets['dtype'].iat[0] == "hamming":
                        charlist = [charmap[letter] for letter in seq]
                        binseq.append(" ".join(charlist))
                    else: # leven
                        binseq.append(seq)

                rand_seqs = self.nmslib_index.knnQueryBatch(binseq, k=2, num_threads=num_threads)
                distlist = []
                for i in rand_seqs:
                    distlist.append(i[1][0])
                zipped = list(zip(seqs, distlist))
                dist_seqs = sorted(zipped, reverse=True, key=lambda x: x[1])
                sort_seq = [item[0] for item in dist_seqs][0:n]

                #sort_dist
                if self.targets['dtype'].iat[0] == "hamming":
                    sort_dist = [item[1] / 2 for item in dist_seqs][0:n]  ### ? does divide by 2 holds for leven???
                else:
                    sort_dist = [item[1] for item in dist_seqs][0:n]  ### ? does divide by 2 holds for leven???

                minimum_hmdist = int(min(sort_dist))
                sm_count += 1
        except IndexError as e:
            raise e

        total_ncontrolsearched = search_mult * n
        self.ncontrolsearched = total_ncontrolsearched
        randomdf = pd.DataFrame(data={"Sequences": sort_seq, "Hamming distance": sort_dist})

        def create_name(seq):
            return "Cont-" + hashlib.md5(seq.encode()).hexdigest()
        randomdf['name'] = randomdf["Sequences"].apply(create_name)
        randomdf = randomdf[["name", "Sequences", "Hamming distance"]]
        randomdf.head()
        return (min(sort_dist),
                statistics.median(sort_dist),
                randomdf)


class Annotation:

    """
    Annotation class for data and methods on targets and gene annotations.

    """

    def __init__(self, annotation_list: List[str], annotation_type: str, target_bed_df: object) -> None:
        """
        Annotation class for data and methods on targets and gene annotations

        Args:
            annotation_list (List[str]): A list of genbank files from a single genome
            annotation_type (str): "genbank" | "gff"
            target_bed_df (object): A pandas dataframe in Bed format with the
                locations of targets in the genome

        Returns:
            None
        """
        self.annotation_list: List[str] = annotation_list
        self.annotation_type = annotation_type
        self.target_bed_df: object = target_bed_df
        self.genbank_bed_df: object = None
        self.feature_dict: Dict = None
        self.nearby: object = None
        self.filtered_df: object = None
        self.qualifiers: object = None

    def check_annotation_type(self):
        """open GTF/GFF and determine if the file provided by the GFF argument is a GFF or GTF file

            Args: None

            Returns (str): ["gff" | "gtf"]
        """
        def search(f):
            line1 = f.readline()
            gffmatch = re.search("gff-version", line1)
            if gffmatch is not None:
                return "gff"
            gtfmatch = re.search("gtf-version", line1)
            if gtfmatch is not None:
                return "gtf"
            else:
                logger.error("Could not verify the GFF/GTF file type. Please make sure your GFF/GTF file starts with '#gtf-version' or '##gff-version'")
                raise ValueError
        testfile = self.annotation_list[0]
        if is_gzip(testfile):
            with gzip.open(testfile, 'rt') as f:
                return search(f)
        else:
            with open(testfile, 'r') as f:
                return search(f)

    def get_annotation_features(self, feature_types: List[str] = None) -> None:
        """
        Parse annotation records into pandas DF/Bed format and dict format saving to self

        Args:
            feature_types (List[str]): a list of Genbank feature types to use

        Returns:
            None
        """
        if feature_types is None:
            feature_types = ["CDS"]
        feature_dict = {}
        pddict = dict(chrom=[], chromStart=[], chromEnd=[], name=[], strand=[])
        if self.annotation_type == "genbank":
            for gbfile in self.annotation_list:
                try:
                    if is_gzip(gbfile):
                        f = gzip.open(gbfile, mode='rt')
                    else:
                        f = open(gbfile, mode='r')
                except IOError as e:
                    logger.error("The genbank file %s could not be opened" % gbfile)
                    raise e
                genbank_file = SeqIO.parse(f, "genbank")
                for entry in genbank_file:
                    for record in entry.features:
                        if record.type in feature_types:
                            if record.location.strand in [1, -1, "+", "-"]:
                                pddict["strand"].append("-" if str(record.location.strand) in ['-1', '-' ] else "+")
                            featid = hashlib.md5(str(record).encode()).hexdigest()
                            pddict['chrom'].append(entry.id)
                            pddict["chromStart"].append(int(record.location.start))
                            pddict["chromEnd"].append(int(record.location.end))
                            pddict["name"].append(featid)
                            for qualifier_key, qualifier_val in record.qualifiers.items():
                                if not qualifier_key in feature_dict:
                                    feature_dict[qualifier_key] = {}
                                feature_dict[qualifier_key][featid] = qualifier_val
            genbankbed = pd.DataFrame.from_dict(pddict)
            self.genbank_bed_df = genbankbed
            self.feature_dict = feature_dict
            f.close()
        elif self.annotation_type == "gff":
            anno_format = self.check_annotation_type()
            for gff in self.annotation_list:
                bedfile = BedTool(gff)
                for rec in bedfile:
                    if rec[2] in feature_types:
                        pddict["chrom"].append(rec[0])
                        pddict["chromStart"].append(rec[3])
                        pddict["chromEnd"].append(rec[4])
                        pddict["strand"].append(rec[6])
                        featid = hashlib.md5(str(rec).encode()).hexdigest()
                        pddict["name"].append(featid)
                        featlist = rec[8].split(';')
                        for feat in featlist:
                            try:
                                if feat.isspace(): # this handles whitespace strings
                                    continue
                                if not feat: # this handles empty strings
                                    continue
                                if anno_format == 'gtf':
                                    fl = re.search('^[^"]*', feat)
                                    fv = re.search('"([^"]*)"', feat)
                                    feat_key = fl.group(0).strip()
                                    feat_val = fv.group(0).strip('"')
                                elif anno_format =='gff':
                                    fl = feat.split('=')
                                    feat_key = fl[0]
                                    feat_val = fl[1]
                                if not feat_key in feature_dict:
                                    feature_dict[feat_key] = {}
                                feature_dict[feat_key][featid] = feat_val
                            except:
                                logger.warning("There appears to be an error in the formatting of an attribute in the "
                                               "record below. Please check your input GFF or GTF file. The record is: {rec} "
                                               "and the attribute is: {att}. Skipping this feature.".format(rec=featlist, att=feat))
                                continue
            genbankbed = pd.DataFrame.from_dict(pddict)
            self.genbank_bed_df = genbankbed
            self.feature_dict = feature_dict


    def _get_qualifiers(self, configpath, excluded: List[str] = None) -> object:
        """
        Create a dataframe with features and their qualifier values

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
        with open(configpath) as cf:
            config = yaml.safe_load(cf)

        min_prop = config['MINIMUM_PROPORTION']

        if excluded is None:
            excluded = ["translation"]
        final_quals = []
        qual_df = pd.DataFrame(data={"Feature id": []})
        for featkey, quals in self.feature_dict.items():
            if len(quals) / len(self.feature_dict[featkey]) > min_prop:
                final_quals.append(featkey)
        for qualifier in final_quals:
            if qualifier not in excluded:
                featlist = []
                quallist = []
                for feat, qual in self.feature_dict[qualifier].items():
                    featlist.append(feat)
                    if isinstance(qual, list):
                        quallist.append(";".join([str(i) for i in qual]))
                    else:
                        quallist.append(qual)
                tempdf = pd.DataFrame({'Feature id': featlist, qualifier: quallist})
                qual_df = qual_df.merge(tempdf, how="outer", on="Feature id")
        self.qualifiers = qual_df

    def _get_nearby_features(self) -> None:
        """
        Adds downstream information to the given target sequences and mapping information

        Args:
            None

        Returns:
            None

        Note:
            Writes a dataframe of nearby features to self.nearby
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
        headers = {0: "Accession", 1: "Guide start", 2: "Guide end", 3: "Guide sequence",
                   4: "Guide strand", 5: "Feature Accession", 6: "Feature start", 7: "Feature end", 8: "Feature id",
                   9: "Feature strand", 10: "Feature distance"}
        downstream: pd.DataFrame = downstream.to_dataframe(disable_auto_names=True, header=None)
        downstream['direction'] = 'downstream'
        upstream = upstream.to_dataframe(disable_auto_names=True, header=None)
        upstream['direction'] = 'upstream'
        upstream = pd.concat([downstream, upstream], axis = 0)
        self.nearby = upstream.rename(columns=headers)


    def _filter_features(self, before_feat: int = 100, after_feat: int = 200 ) -> None:
        """
        Merge targets with Feature list and filter for guides close enough to interact.

        Args:
            before_feat (int): The maximum distance before the start of a feature measured from closest point to guide
            after_feat (int): The maximum distance after the start codon (into the gene)

        Returns:
            None
        """
        # for guides in the same orientation as the targets ( +/+ or -/-) select guides that are within
        #  before_feat of the gene start
        filtered_df = self.nearby.query(
            '`Guide strand` == `Feature strand` and 0 < `Feature distance` < @before_feat')

        # for guides in the +/+ orientation select guides where the end is within [before_feat] of the gene start
        p1 = (self.nearby.query('`Guide strand` == "+" and `Feature strand` == "+" \
                                 and `Feature distance` == 0 \
                                 and `Guide end` - `Feature start` < @after_feat \
                                 and `Guide end` <= `Feature end`'))

        # for guides in the -/- orientation select guides where the end is within [before_feat] of the gene start
        p2 = (self.nearby.query('`Guide strand` == "-" and `Feature strand` == "-" \
                                 and `Feature distance` == 0 \
                                 and `Feature end` - `Guide start` < @after_feat \
                                 and `Guide end` <= `Feature end`'))

        # Select guides where target is + and guide is - and the guide is infront of the gene
        p3 = (self.nearby.query('`Guide strand` == "-" and `Feature strand` == "+" and \
                                 0 <`Feature start` - `Guide end` < @before_feat'))

        # Select guides where target is - and guide is + and the guide is infront of the gene
        p4 = (self.nearby.query('`Guide strand` == "+" and `Feature strand` == "-" and \
                                 0 <`Guide start` - `Feature end` < @before_feat'))

        # Select guides where target is + and guide is - and the guide is is within [before_feat] of the gene start
        p5 = (self.nearby.query('`Guide strand` == "-" and `Feature strand` == "+" and \
                                 0 <`Guide end` -`Feature start` < @after_feat \
                                 and `Guide end` <= `Feature end`'))

        # Select guides where target is - and guide is + and the guide is is within [before_feat] of the gene start
        p6 = (self.nearby.query('`Guide strand` == "+" and `Feature strand` == "-" and \
                                 0 <`Feature end` - `Guide start` < @after_feat \
                                 and `Guide end` <= `Feature end`'))
        self.filtered_df = pd.concat([filtered_df, p1, p2, p3, p4, p5, p6], axis=0)

    def _format_guide_table(self, targetprocessor_object) -> pd.DataFrame:
        """
        Create guide table for output

        Args:
            target- a dataframe with targets from targetclass

        Returns:
            (PandasDataFrame): A formated pandas dataframe
        """
        def gc(seq):
            cnt = 0
            for letter in seq:
                if letter in ["G", "C"]:
                    cnt += 1
            return cnt / len(seq)

        def get_guide_hash(seq):
            return hashlib.md5(seq.encode()).hexdigest()

        def checklen30(seq):
            if len(seq) == 30:
                return True
            return False

        def get_off_target_score(seq):
            dlist = targetprocessor_object.neighbors[seq]["neighbors"]["dist"]
            s = [str(i) for i in dlist]
            return ";".join(s)

        def get_off_target_seqs(seq):
            slist = targetprocessor_object.neighbors[seq]["neighbors"]["seqs"]
            return ";".join(slist)
        pretty_df = deepcopy(self.filtered_df)  # anno class object
        # retrive only guides that are in neighbors keys.
        pretty_df = pretty_df[pretty_df["Guide sequence"].isin(
            list(targetprocessor_object.neighbors.keys()))]
        pretty_df['GC'] = pretty_df['Guide sequence'].apply(gc)
        pretty_df['Guide name'] = pretty_df['Guide sequence'].apply(get_guide_hash)
        pretty_df['Target strand'] = np.where(
            pretty_df['Guide strand'] == pretty_df['Feature strand'], 'coding', 'non-coding')
        pretty_df['Similar guide distances'] = pretty_df['Guide sequence'].apply(
            get_off_target_score)
        pretty_df['Similar guides'] = pretty_df['Guide sequence'].apply(get_off_target_seqs)
        pretty_df = pd.merge(pretty_df, targetprocessor_object.targets, how="left",
         left_on=['Guide sequence', 'Guide start', 'Guide end', 'Accession'],
            right_on=['target', 'start', 'stop', 'seqid'])
        # rename exact_pam to PAM
        pretty_df = pretty_df.rename(columns={"exact_pam": "PAM"})

        pretty_df = pretty_df[['Guide name', 'Guide sequence', 'GC', 'dtype', 'Accession', 'Guide start', 'Guide end',
                    'Guide strand', 'PAM', 'Feature id',
                    'Feature start', 'Feature end', 'Feature strand',
                    'Feature distance', 'Similar guides', 'Similar guide distances','target_seq30']]
        pretty_df = pretty_df.merge(self.qualifiers, how="left", on="Feature id")
        pretty_df = pretty_df.sort_values(by=['Accession', 'Feature start'])
        # to match with the numbering with other tools- offset
        pretty_df['Guide start'] = pretty_df['Guide start'] + 1
        pretty_df['Feature start'] = pretty_df['Feature start'] + 1
        pretty_df=pretty_df.loc[pretty_df['target_seq30'].apply(checklen30)==True]
        self.pretty_df = pretty_df

    def _filterlocus(self, attribute:str , filter_by_locus:list = []) -> pd.DataFrame:
        """
        Create guide table for output for a selected attribute type

        Args:
            attribute: The key in the attributes column (column 9) of the GFF/GTF file to filter on
            filter_by_locus: A list of Identifiers to filter the full data frame by

        Returns:
            (PandasDataFrame): A formated pandas dataframe
        """

        df = deepcopy(self.pretty_df)  # anno class object
        if len (filter_by_locus) > 0:
            df = df[df[attribute].isin(filter_by_locus)]
        return df

    def locuslen(self) -> int:
        """
        Count the number of locus tag in the genebank file

        Args:
            None

        Returns:
            (int): Number of locus tag
        """
        da_keys = self.feature_dict.keys()
        firsttag = (list(da_keys)[0])
        if firsttag:
            locus_count = len(self.feature_dict[firsttag].keys())
            return firsttag, locus_count
        else:
            logger.warning("A locus key could not be found.")
            return "notag", 0



class GuideMakerPlot:

    """
    A class with functions to plot guides over genome cooridinates.

    """

    def __init__(self, prettydf: pd.DataFrame, outdir: str) -> None:
        """
        GuideMakerPlot class for visualizing distrubution of gRNA, features, and locus.

        Args:
            prettydf (PandasDataFrame): Final output from GuideMaker
            outdir (str): Output Directory

        Returns:
            None
        """
        self.prettydf = prettydf
        self.accession = list(set(self.prettydf['Accession']))

        def _singleplot(df):
            """
            Returns guidemaker plot describing PAM targets

            Args:
                df(PandasDataFrame): Final output from GuideMaker for a single accession

            Return:
                None
            """
            source = df
            brush = alt.selection(type='interval', encodings=['x'])
            binNum = int(round(source['Feature end'].max() / 200, 0))
            display_info = source.columns.tolist()

            # Feature density
            densityF = alt.Chart(source).transform_density(
            'Feature start',
            as_=['Feature start', 'Feature Density'],
            extent=[1, source['Feature end'].max()],
            bandwidth=binNum,
            ).mark_area(color='black', opacity=0.6).encode(
            x=alt.X('Feature start', axis=alt.Axis(title='Genome Coordinates (bp)', tickCount=5)),
            y='Feature Density:Q',
            ).properties(height=50, width=500)

            # Guide density
            densityG = alt.Chart(source).transform_density(
            'Guide start',
            as_=['Guide start', 'Guide Density'],
            extent=[1, source['Feature end'].max()],
            bandwidth=binNum,
            ).mark_area(color='pink', opacity=0.6).encode(
            x=alt.X('Guide start', axis=alt.Axis(title='Genome Coordinates (bp)', tickCount=5)),
            y='Guide Density:Q',
            ).properties(height=50, width=500).add_selection(brush)

            # locus bar
            locus = alt.Chart(source).mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3).encode(
            x='count(locus_tag):Q',
            y=alt.Y('locus_tag', axis=alt.Axis(title='Locus')),
            color='PAM:N',
            tooltip=display_info
            ).transform_filter(
            brush
            ).interactive().properties(height=500, width=500)
            guidemakerChart = (densityF & densityG & locus)
            return(guidemakerChart)

        for accession in self.accession:
            df = self.prettydf[self.prettydf['Accession'] == accession]
            accession_plot = _singleplot(df)
            plot_file_name = f"{outdir}/{accession}.html"
            accession_plot.save(plot_file_name)


def get_fastas(filelist, input_format="genbank", tempdir=None):
    """
    Saves a Fasta and from 1 or more Genbank files (may be gzipped)

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
                        records = SeqIO.parse(f, input_format)
                        SeqIO.write(map(lambda obj: obj.upper(), records), f1, "fasta")
                else:
                    with open(file, 'r') as f:
                        records = (SeqIO.parse(f, input_format))
                        SeqIO.write(map(lambda obj: obj.upper(), records), f1, "fasta")
        return fastpath
    except Exception as e:
        logger.exception("An error occurred in the input file %s" % file)
        raise e


def extend_ambiguous_dna(seq: str) -> List[str]:
    """
    Return list of all possible sequences given an ambiguous DNA input

    Args:
        seq(str): A DNA string

    Return:
        List[str]: A list of DNA string with expanded ambiguous DNA values
    """
    ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }
    extend_list = []
    for i in product(*[ambiguous_dna_values[j] for j in seq]):
        extend_list.append("".join(i))
    return extend_list


# add CDF and Doench Azimuth scores

def cfd_score(df):
    def cfd_calculator(knnstrlist, guide, mm_scores):
            knnlist =  knnstrlist.split(';')
            cfd_list=[]
            for item in knnlist:
                result=cfd_score_calculator.calc_cfd(guide, item, mm_scores=mm_scores)
                cfd_list.append(result)
            s = [str(i) for i in cfd_list]
            return s

    def get_max_cfd(cfdlist):
        newlist = [float(x) for x in cfdlist]
        newlist.sort()
        maxcfd = newlist[-1]
        return(maxcfd)
    mm_scores, _ = cfd_score_calculator.get_mm_pam_scores()
    df['CFD Similar Guides'] = df.apply(lambda x: cfd_calculator(x['Similar guides'], x['Guide sequence'], mm_scores=mm_scores), axis=1)
    # Add a column with max CFD score
    df['Max CFD'] = df['CFD Similar Guides'].apply(get_max_cfd)
    return df



def get_doench_efficiency_score(df, pam_orientation, num_threads=1):
    checkset={'AGG','CGG','TGG','GGG'}
    # filter out lines with N'safter the PAM, these cannot be scored
    df2 = df[-df.target_seq30.str.contains("WSKMYRVHDBN")]
    if len(df) != len(df2):
        n_removed = len(df) - len(df2)
        logger.warning("{} guides were removed from consideration becasue there were ambiguous nucleotides in the region flanking the PAM site. These cannot be scored.".format(n_removed) )
    if pam_orientation == "3prime" and set(df2.PAM)==checkset:

        doenchscore = doench_predict.predict(np.array([x.upper() for x in df2.target_seq30]), num_threads=num_threads)
        df2["Efficiency"] = doenchscore
    else:
        logger.warning("NOTE: doench_efficiency_score based on Doench et al. 2016 - can only  be used for NGG PAM).Check PAM sequence and PAM orientation")
        df2["Efficiency"] = "Not Available"
    return df2.drop('target_seq30', axis=1)
