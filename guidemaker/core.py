"""Core classes and functions for GuideMaker."""
import os
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
from Bio.SeqUtils import GC
from pybedtools import BedTool
from Bio import Seq
from copy import deepcopy
import pandas as pd
import numpy as np
import altair as alt


logger = logging.getLogger('guidemaker.core')
PandasDataFrame = TypeVar('pandas.core.frame.DataFrame')


def is_gzip(filename: str):
    try:
        with open(filename, "rb") as f:
            logging.info("check if %s is gzipped" % filename)
            return f.read(2) == b'\x1f\x8b'
    except IOError as e:
        logging.error("Could not open the file %s to determine if it was gzipped" % filename)
        raise e


class PamTarget:
    """A Class representing a Protospacer Adjacent Motif (PAM) and targets. The class includes all targets for given PAM as a dataframe, PAM and target attributes, and methods to find target and control sequences."""

    def __init__(self, pam: str, pam_orientation: str) -> None:
        """Pam __init__

        Args:
            pam (str): A DNA string in ambiguous IUPAC format
            pam_orientation (str): [5prime | 3prime ]
                5prime means the order is 5'-[pam][target]-3'
                3prime means the order is 5'-[target][pam]-3'

        Returns:
            None
        """
        for letter in pam.upper():
            assert letter in ['A', 'C', 'G', 'T', 'M', 'R', 'W',
                'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'X', 'N']
        assert pam_orientation in ["3prime", "5prime"]
        self.pam: str = pam.upper()
        self.pam_orientation: str = pam_orientation

    def __str__(self) -> str:
        return "A PAM object: {self.pam}".format(self=self)

    def find_targets(self, seq_record_iter: object, target_len: int) -> PandasDataFrame:
        """Find all targets on a sequence that match for the PAM on both strand(s)

        Args:
            seq_record_iter (object): A Biopython SeqRecord iterator from SeqIO.parse
            target_len (int): The length of the target sequence

        Returns:
            PandasDataFrame: A pandas dataframe with of matching targets
        """

        def reverse_complement(seq: str) -> str:
            """Reverse complement of the PAM sequence

            Args:
                seq (str): A DNA string

            Returns:
                str: A reverse complement of DNA string
            """
            bpseq = Seq.Seq(seq)
            return str(bpseq.reverse_complement())

        def pam2re(pam: str) -> str:
            """Convert an IUPAC ambiguous PAM to a Regex expression

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
            """Check targets for guidelength and DNA bases

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
            """Search for guides with 5prime pam orientation in the forward strand

            Args:
                pam_pattern (str): A DNA string representing PAM
                dnaseq (str): A DNA string representing genome
                target_len (int): Guide length

            Returns:
                (Generator): A generator with target_seq, exact_pam, start, stop, strand, and pam_orientation
            """
            for match_obj in regex.finditer(pattern=pam_pattern, string=dnaseq, overlapped=True):
                target_seq = dnaseq[match_obj.end(): match_obj.end() + target_len]
                if check_target(target_seq, target_len):
                    exact_pam = match_obj.group(0)
                    start = match_obj.end()
                    stop = match_obj.end() + target_len
                    # 5prime =True, 3prime = False
                    pam_orientation = True
                    # forward =True, reverse = False
                    strand = True
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation

        def run_for_3p(pam_pattern, dnaseq, target_len) -> Generator:
            """Search for guides with 3prime pam orientation in the reverse strand

            Args:
                pam_pattern (str): A DNA string representing PAM
                dnaseq (str): A DNA string representing genome
                target_len (int): Guide length

            Returns:
                (Generator): A generator with target_seq, exact_pam, start, stop, strand, and pam_orientation
            """
            for match_obj in regex.finditer(pattern=pam_pattern, string=dnaseq, overlapped=True):
                target_seq = seq[match_obj.start() - target_len: match_obj.start()]
                if check_target(target_seq, target_len):
                    exact_pam = match_obj.group(0)
                    start = match_obj.start() - target_len
                    stop = match_obj.start()
                    # 5prime =True, 3prime = False
                    pam_orientation = False
                    # forward =True, reverse = False
                    strand = True
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation

        def run_rev_5p(pam_pattern, dnaseq, target_len) -> Generator:
            """Search for guides with 5prime pam orientation in the reverse strand

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
                if check_target(target_seq, target_len):
                    exact_pam = reverse_complement(match_obj.group(0))
                    start = match_obj.start() - target_len
                    stop = match_obj.start()
                    # 5prime =True, 3prime = False
                    pam_orientation = True
                    # forward =True, reverse = False
                    strand = False
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation

        def run_rev_3p(pam_pattern, dnaseq, target_len) -> Generator:
            """Search for guides with 3prime pam orientation in the reverse strand

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
                if check_target(target_seq, target_len):
                    exact_pam = reverse_complement(match_obj.group(0))
                    start = match_obj.end()
                    stop = match_obj.end() + target_len
                    # 5prime =True, 3prime = False
                    pam_orientation = False
                    # forward =True, reverse = False
                    strand = False
                    yield target_seq, exact_pam, start, stop, strand, pam_orientation

        target_list = []
        for record in seq_record_iter:
            record_id = record.id
            seq = str(record.seq)
            if self.pam_orientation == "5prime":
                # forward
                for5p = pd.DataFrame(run_for_5p(pam2re(self.pam), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation"])
                for5p["seqid"] = record_id
                # string to boolean conversion is not straight - as all string were set to Trues- so change the encoding in functions above.
                # https://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python/715455#715455
                for5p = for5p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(for5p)
                # reverse
                rev5p = pd.DataFrame(run_rev_5p(pam2re(reverse_complement(self.pam)), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation"])
                rev5p["seqid"] = record_id
                rev5p = rev5p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(rev5p)
                # Question? Append directly vs. concat then append? https://ravinpoudel.github.io/AppendVsConcat/
            elif self.pam_orientation == "3prime":
                # forward
                for3p = pd.DataFrame(run_for_3p(pam2re(self.pam), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation"])
                for3p["seqid"] = record_id
                for3p = for3p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(for3p)
                # reverse
                rev3p = pd.DataFrame(run_rev_3p(pam2re(reverse_complement(self.pam)), seq, target_len), columns=[
                                     "target", "exact_pam", "start", "stop", "strand", "pam_orientation"])
                rev3p["seqid"] = record_id
                rev3p = rev3p.astype({"target": 'str', "exact_pam": 'category', "start": 'uint32',
                                     "stop": 'uint32', "strand": 'bool', "pam_orientation": 'bool', "seqid": 'category'})
                target_list.append(rev3p)
            gc.collect()  # clear memory after each chromosome
        df_targets = pd.concat(target_list, ignore_index=True)
        df_targets = df_targets.assign(seedseq=np.nan, isseedduplicated=np.nan)
        return df_targets


class TargetProcessor:
    """A Class representing a set of guide RNA targets

    The class includes all targets in a dataframe, methods to process target and a dict with edit distances for sequences.

    """

    def __init__(self, targets, lsr: int, hammingdist: int = 2, knum: int = 2) -> None:
        self.targets = targets  # pandas dataframe
        self.lsr: int = lsr  # length of seed region
        self.hammingdist: int = hammingdist
        self.knum: int = knum
        self.nmslib_index: object = None
        self.neighbors: dict = {}
        self.ncontrolsearched: int = None
        self.gc_percent: float = None
        self.genomesize: float = None
        self.pam_orientation: bool = targets['pam_orientation'].iat[0]

    def __str__(self):
        info = "TargetList: contains a set of {} potential PAM targets".format(len(self.targets))
        return info

    def __len__(self):
        return len(self.targets)

    def check_restriction_enzymes(self, restriction_enzyme_list: list = []) -> None:
        """Check for restriction enzymes and its reverse complement within gRNA sequence

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
        element_to_exclude = sum(element_to_exclude, [])  # flatout list of list to list
        if len(element_to_exclude) > 0:
            self.targets = self.targets.loc[self.targets['target'].str.contains(
                '|'.join(element_to_exclude)) == False]
        else:
            self.targets

    def _one_hot_encode(self, seq_list: List[object]) -> List[str]:
        """One hot encode Target DNA as a binary string representation for NMSLIB

        """
        charmap = {'A': '1 0 0 0', 'C': '0 1 0 0', 'G': '0 0 1 0', 'T': '0 0 0 1'}

        def seq_to_bin(seq):
            charlist = [charmap[letter] for letter in seq]
            return " ".join(charlist)
        return list(map(seq_to_bin, seq_list))

    def find_unique_near_pam(self) -> None:
        """identify unique sequences in the target list

        The function filters a list of Target objects for targets that
        are unique in the region closest to the PAM. The region length is defined
        by the lsr (length of seed region that need to be unique).

        Args:
            lsr (int): Length of seed region that is close to PAM

        Returns:
            Nome
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
        """Create nmslib index

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
        #logging.info("unique targets for index: %s" % len(notduplicated_targets))
        bintargets = self._one_hot_encode(notduplicated_targets)
        index_params = {'M': M, 'indexThreadQty': num_threads, 'efConstruction': efC, 'post': post}
        index = nmslib.init(space='bit_hamming',
                            dtype=nmslib.DistType.INT,
                            data_type=nmslib.DataType.OBJECT_AS_STRING,
                            method='hnsw')
        index.addDataPointBatch(bintargets)
        index.createIndex(index_params, print_progress=True)
        self.nmslib_index = index

    def get_neighbors(self, configpath, num_threads=2) -> None:
        """Get nearest neighbors for sequences removing sequences that
         have neighbors less than the Hamming distance threshold

        For the list of all targets calculate the (knum) nearest neighbors.
        filter out targets with close neighbors and
        Writes a dictionary to self.neighbors:
        self.neighbors[seq]{target: seq_obj, neighbors: {seqs:[s1, s1, ...], dist:[d1, d1,...]}}

        Args: 
            None

        Returns:
            None
        """
        with open(configpath) as cf:
            config = yaml.safe_load(cf)

        ef = config['NMSLIB']['ef']

        unique_targets = self.targets.loc[self.targets['isseedduplicated']
            == False]['target'].tolist()
        unique_bintargets = self._one_hot_encode(unique_targets)  # search unique seed one
        self.nmslib_index.setQueryTimeParams({'efSearch': ef})
        results_list = self.nmslib_index.knnQueryBatch(unique_bintargets,
                                               k=self.knum, num_threads=num_threads)
        neighbor_dict = {}
        for i, entry in enumerate(results_list):
            queryseq = unique_targets[i]
            hitseqidx = entry[0].tolist()
            hammingdist = entry[1].tolist()
            # here we just check if the first element of hammingist list is >= 2 * self.hammingdist, as list is sorted- if first fails whole fails
            # to close guides.
            # this should be 0 or 1?
            # this should be 1 == b/c each guides will have exact match with itself at 0 position.
            if hammingdist[1] >= 2 * self.hammingdist:  # multiply by 4 b/c each base is one hot encoded in 4 bits
                neighbors = {"seqs": [self.targets['target'].values[x] for x in hitseqidx],  # reverse this?
                             "dist": [int(x / 2) for x in hammingdist]}
                neighbor_dict[queryseq] = {"target": unique_targets[i],
                                           "neighbors": neighbors}
        self.neighbors = neighbor_dict

    def export_bed(self) -> object:
        """Export the targets in self.neighbors to a bed format file

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
                         num_threads: int = 2) -> PandasDataFrame:
        """Create random sequences with a specified GC probability and find seqs with the greatest
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
            gccnt += GC(record.seq) * len(record)
            totlen += len(record)
        gc = gccnt / (totlen * 100)
        #print("Percentage of GC content in the input genome: "+"{:.2f}".format(gc * 100))
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
                    charlist = [charmap[letter] for letter in seq]
                    binseq.append(" ".join(charlist))
                rand_seqs = self.nmslib_index.knnQueryBatch(binseq, k=2, num_threads=num_threads)
                distlist = []
                for i in rand_seqs:
                    distlist.append(i[1][0])
                zipped = list(zip(seqs, distlist))
                dist_seqs = sorted(zipped, reverse=True, key=lambda x: x[1])
                sort_seq = [item[0] for item in dist_seqs][0:n]
                sort_dist = [item[1] / 2 for item in dist_seqs][0:n]
                minimum_hmdist = int(min(sort_dist))
                sm_count += 1
        except IndexError as e:
           # print("Number of random control searched: ", search_mult * n)
            pass

        total_ncontrolsearched = search_mult * n
        self.ncontrolsearched = total_ncontrolsearched
        randomdf = pd.DataFrame(data={"Sequences": sort_seq, "Hamming distance": sort_dist})

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
        pddict = dict(chrom=[], chromStart=[], chromEnd=[], name=[], strand=[])
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
                        if record.strand in [1, -1, "+", "-"]:
                            pddict["strand"].append("-" if record.strand < 0 else "+")
                            featid = hashlib.md5(str(record).encode()).hexdigest()
                            pddict['chrom'].append(entry.id)
                            pddict["chromStart"].append(record.location.start.position)
                            pddict["chromEnd"].append(record.location.end.position)
                            pddict["name"].append(featid)
                            for qualifier_key, qualifier_val in record.qualifiers.items():
                                if not qualifier_key in feature_dict:
                                    feature_dict[qualifier_key] = {}
                                feature_dict[qualifier_key][featid] = qualifier_val
            genbankbed = pd.DataFrame.from_dict(pddict)
            self.genbank_bed_df = genbankbed
            self.feature_dict = feature_dict
        f.close()

    def _get_qualifiers(self, configpath, excluded: List[str] = None) -> object:
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
        with open(configpath) as cf:
            config = yaml.safe_load(cf)

        min_prop = config['MINIMUM_PROPORTION']

        if excluded is None:
            excluded = ["translation"]
        final_quals = []
        qual_df = pd.DataFrame(data={"Feature id": []})
        for quals in self.feature_dict:
            if len(quals) / len(self.feature_dict) > min_prop:
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
        upstream = upstream.append(downstream)
        self.nearby = upstream.rename(columns=headers)

    def _filter_features(self, before_feat: int = 100, after_feat: int = 200) -> None:
        """Merge targets with Feature list and filter for guides close enough to interact.

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
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "+" and `Feature strand` == "+" \
                                             and `Feature distance` == 0 and \
                                             `Guide end` - `Feature start` < @after_feat'))
        # for guides in the -/- orientation select guides where the end is within [before_feat] of the gene start
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "-" and `Feature strand` == "-" \
                                                     and `Feature distance` == 0 \
                                                     and `Feature end` - `Guide start` < @after_feat'))
        # Select guides where target is + and guide is - and the guide is infront of the gene
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "-" and `Feature strand` == "+" and \
                                                     0 <`Feature start` - `Guide end` < @before_feat'))
        # Select guides where target is - and guide is + and the guide is infront of the gene
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "+" and `Feature strand` == "-" and \
                                                     0 <`Guide start` - `Feature end` < @before_feat'))
        # Select guides where target is + and guide is - and the guide is is within [before_feat] of the gene start
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "-" and `Feature strand` == "+" and \
                                                             0 <`Guide end` -`Feature start`  < @after_feat'))
        # Select guides where target is - and guide is + and the guide is is within [before_feat] of the gene start
        filtered_df = filtered_df.append(self.nearby.query('`Guide strand` == "+" and `Feature strand` == "-" and \
                                                             0 <`Feature end` - `Guide start` < @after_feat'))

        self.filtered_df = filtered_df

    def _format_guide_table(self, targetprocessor_object) -> PandasDataFrame:
        """Create guide table for output

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

        pretty_df = pretty_df[['Guide name', "Guide sequence", 'GC', "Accession", "Guide start", "Guide end",
                    "Guide strand", 'PAM', "Feature id",
                    "Feature start", "Feature end", "Feature strand",
                    "Feature distance", 'Similar guides', 'Similar guide distances']]
        pretty_df = pretty_df.merge(self.qualifiers, how="left", on="Feature id")
        pretty_df = pretty_df.sort_values(by=['Accession', 'Feature start'])
        # to match with the numbering with other tools- offset
        pretty_df['Guide start'] = pretty_df['Guide start'] + 1
        pretty_df['Feature start'] = pretty_df['Feature start'] + 1
        return pretty_df

    def locuslen(self) -> int:
        """ Count the number of locus tag in the genebank file

        Args:
            None

        Returns:
            (int): Number of locus tag
        """

        locus_count = len(self.feature_dict['locus_tag' or 'locus'].keys())
        return(locus_count)


class GuideMakerPlot:
    def __init__(self, prettydf: PandasDataFrame, outdir: str) -> None:
        """GuideMakerPlot class for visualizing distrubution of gRNA, features, and locus.

        Args:
            prettydf (PandasDataFrame): Final output from GuideMaker
            outdir (str): Output Directory

        Returns:
            None
        """
        self.prettydf = prettydf
        self.accession = list(set(self.prettydf['Accession']))

        def _singleplot(df):
            """Returns guidemaker plot describing PAM targets

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


def extend_ambiguous_dna(seq: str) -> List[str]:
    """Return list of all possible sequences given an ambiguous DNA input

    Args:
        seq(str): A DNA string

    Return:
        List[str]: A list of DNA string with expanded ambiguous DNA values
    """
    dna_dict = Seq.IUPAC.IUPACData.ambiguous_dna_values
    extend_list = []
    for i in product(*[dna_dict[j] for j in seq]):
        extend_list.append("".join(i))
    return extend_list
