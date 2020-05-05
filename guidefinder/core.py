"""Core classes and functions for Guidefinder

"""
import os
from typing import List, Set, Dict, Tuple
from itertools import product, tee, chain
import gzip
import hashlib
from Bio.Seq import Seq
from Bio import SeqIO
import nmslib
from pybedtools import BedTool
import pandas as pd
import uuid



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
            assert letter in ['G', 'A', 'T', 'C', 'N' ]
        assert pam_orientation in ["3prime", "5prime"]
        self.pam: str = pam.upper()
        self.pam_orientation: str = pam_orientation

    def __str__(self) -> str:
        return "A PAM object: {self.pam}".format(self=self)

    def reverse_complement(self) -> object:
        """reverse compliment of the PAM sequence
        Args:
            None
        Returns:
            str: a new Pam object Reverse complemented

        """
        pamseq = Seq(self.pam)
        return Pam(pam=str(pamseq.reverse_complement()),
                   pam_orientation=self.pam_orientation)

    def extend_ambiguous_dna(self) -> Set[str]:
        """convert ambiguous DNA input to return frozen set of all sequences

        Returns:
            (frozenset): all possible sequences given an ambiguous DNA input

        """
        dnaval = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
                  'M': 'AC', 'R': 'AG', 'W': 'AT', 'S': 'CG',
                  'Y': 'CT', 'K': 'GT', 'V': 'ACG', 'H': 'ACT',
                  'D': 'AGT', 'B': 'CGT', 'X': 'GATC', 'N': 'GATC'}
        dnalist = []
        for i in product(*[dnaval[j] for j in self.pam]):
            dnalist.append("".join(i))
        return frozenset(dnalist)

    def find_targets(self, seqrecord_obj: object, strand: str, target_len: int) -> List[object]:
        """Find all targets on a sequence that match for the PAM on the requested strand(s)

        Args:
            seqrecord_obj (object): A Biopython SeqRecord instance
            strand (str): The strand to search choices: ["forward", "reverse", "both"]
            target_len (int): The length of the target sequence
        Returns:
            list: A list of Target class instances

        """
        fset = self.extend_ambiguous_dna()
        rset = self.reverse_complement().extend_ambiguous_dna()

        def window(iterable, size):
            iters = tee(iterable, size)
            for i in range(1, size):
                for each in iters[i:]:
                    next(each, None)
            return zip(*iters)
#                5prime means the order is 5'-[pam][target]-3'
#                3prime means the order is 5'-[target][pam]-3'
        def compute_seq_coords(self, hitset, strand, target_len, seqrecord_obj, i):
            if hitset == "fset" and strand in ("forward", "both") and self.pam_orientation == "3prime":
                start = i - target_len
                stop = i
                seq = str(seqrecord_obj[start:stop].seq)
                exact_pam = str(seqrecord_obj[stop:(stop + len(self.pam))].seq)
            elif hitset == "fset" and strand in ("forward", "both") and self.pam_orientation == "5prime":
                start = i + len(self.pam)
                stop = i + len(self.pam)  + target_len
                seq = str(seqrecord_obj[start:stop].seq)
                exact_pam = str(seqrecord_obj[start - len(self.pam):(start)].seq)
            elif hitset == "rset" and strand in ("reverse", "both") and self.pam_orientation == "3prime":
                start = i + len(self.pam)
                stop = i + len(self.pam) + target_len
                seq = str(seqrecord_obj[start:stop].seq.reverse_complement())
                exact_pam == str(seqrecord_obj[stop:(stop + len(self.pam))].seq.reverse_complement())
            elif hitset == "rset" and strand in ("reverse", "both") and self.pam_orientation == "5prime":
                start = i - target_len
                stop = i
                seq = str(seqrecord_obj[start:stop].seq.reverse_complement())
                exact_pam =str(seqrecord_obj[start - len(self.pam):(start)].seq.reverse_complement())
            else:
                return None
            if 0 <= start <= len(seqrecord_obj) and  0 <= stop <= len(seqrecord_obj):
                return Target(seq=seq,
                              exact_pam=exact_pam,
                              strand=strand,
                              pam_orientation=self.pam_orientation,
                              seqid=seqrecord_obj.id,
                              start=start,
                              stop=stop)

        # Create iterable of PAM length windows across the sequence
        target_list = []
        kmer_iter = window(str(seqrecord_obj.seq), len(self.pam))
        for i, kmer in enumerate(kmer_iter):
            hitset = None
            if ''.join(kmer) in fset:
                hitset = "fset"
            elif ''.join(kmer) in rset:
                hitset = "rset"
            else:
                continue
            tar = compute_seq_coords(self, hitset=hitset,
                                     strand=strand,
                                     target_len=target_len,
                                     seqrecord_obj=seqrecord_obj,
                                     i=i)
            if tar:
                target_list.append(tar)
        return target_list


class Target:
    """A class representing a candidate target sequence for a PAM

    This is an object for holding data on possible target sequences
    adjacent to PAM sites.
    """
    def __init__(self, seq: str, exact_pam: str, strand: str, pam_orientation: str,
                 seqid: str, start: int, stop: int) -> None:
        self.seq: str= seq
        self.exact_pam: str = exact_pam
        self.strand: str = strand
        self.pam_orientation: str = pam_orientation
        self.seqid: str = seqid
        self.start: int = start
        self.stop: int = stop
        self.md5: str  = hashlib.md5(seq.encode()).hexdigest()

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
    def __init__(self, targets: List, lcp: int, levindist: int=2, knum: int=2) -> None:
        self.targets: List = targets
        self.lcp: int = lcp
        self.levindist: int = levindist
        self.knum: int = knum
        self.unique_targets: dict = {}
        self.nmslib_index: object = None
        self.neighbors: dict = {}

    def __str__(self):
        info = "TargetList: contains a set of {} potential PAM targets, {} \
                targets unique near the pam site and {} targets with an edit \
                distance of >= {}".format(len(self.targets),
                                              len(self.unique_targets),
                                              len(self.neighbors),
                                              self.levindist)
        return info

    def __len__(self):
        return len(self.targets)

    def find_unique_near_pam(self) -> None:
        """Filter a list of Target objects filtering for sequences

        The function filters a list of Target objects for targets that
        are unique in the region closest to the PAM. The region length is defined
        by the lcp.

        Args:
            lcp (int): Length of conserved sequence close to PAM
        """
        unique_lcp_list = []
        for target in self.targets:
            if getattr(target, "pam_orientation") == "5prime":
                proximal = target.seq[(len(target) - self.lcp):]
            elif getattr(target, "pam_orientation") == "3prime":
                proximal = target.seq[0:self.lcp]
            unique_lcp_list.append(proximal)
        unique_lcp_set = set(unique_lcp_list)
        final_dict = {}
        for target in self.targets:
            if getattr(target, "pam_orientation") == "5prime":
                proximal2 = target.seq[(len(target) - self.lcp):]
            elif getattr(target, "pam_orientation") == "3prime":
                proximal2 = target.seq[0:self.lcp]
            if proximal2 in  unique_lcp_set:
                final_dict[target.seq] = target
        self.unique_targets = final_dict

    def create_index(self):
        """Initializes and returns a NMSLIB index

        Args:
            strings (str): Strings to calculate distance

        Returns:
            None (but writes NMSLIB index to self)
        """
        seq_list = []
        index = nmslib.init(space='leven',
                            dtype=nmslib.DistType.INT,
                            data_type=nmslib.DataType.OBJECT_AS_STRING,
                            method='small_world_rand')
        index.addDataPointBatch(list(self.unique_targets.keys()))
        index.createIndex()
        self.nmslib_index = index

    def get_neighbors(self) -> None:
        """Get nearest neighbors for sequences.

        For the unique_sequences calculate the (knum) nearest neighbors.
        Writes a dictionary to self.neighbors:
        self.neighbors[seq]{target: seq_obj, neighbors: {seqs:[s1, s1, ...], dist:[d1, d1,...]}}

        Args: None
        Returns: None
        """
        seqs = list(self.unique_targets.keys())
        results_list = self.nmslib_index.knnQueryBatch(seqs,
                                               k=self.knum)
        neighbor_dict = {}
        for i, entry in enumerate(results_list):
            queryseq = seqs[ i - 1]
            hitseqidx = list(entry[0])
            editdist = list(entry[1])
            if editdist[1] >= self.levindist:
                neighbors = {"seqs" : [seqs[x] for x in hitseqidx],
                             "dist": editdist }
                neighbor_dict[queryseq] = {"target": self.unique_targets[queryseq],
                                           "neighbors": neighbors}
        self.neighbors =  neighbor_dict

    def export_bed(self) -> object:
        """export the targets in self.neighbors to a bed format file

        Args:
            file (str): the name and location of file to export

        Returns:
            (obj): A Pandas Dataframe in Bed format
        """
        chrom = []
        chromstart = []
        chromend = []
        name = []
        score = []
        strand = []
        for rec in self.neighbors.values():
            chrom.append(rec["target"].seqid)
            chromstart.append(rec["target"].start)
            chromend.append(rec["target"].stop)
            name.append(rec["target"].seq)
            score.append(0)
            if rec["target"].strand == "forward":
                strand.append("+")
            elif rec["target"].strand == "reverse":
                strand.append("-")
        df = pd.DataFrame.from_dict({"chrom": chrom,
                                     "chromstart":chromstart,
                                     "chromend": chromend,
                                     "name": name,
                                     "score": 0,
                                     "strand": strand})
        df.sort_values(by=['chrom','chromstart'], inplace=True)
        # df = df.replace(np.nan, "NA")
        return df
        #df.to_csv(path=outfile, sep="\t", index=False, header=False)


class Annotation:
    def __init__(self, genbank_list: List[str], target_bed_df: object ) -> None:
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

    def _get_genbank_features(self, feature_types: List[str]=["CDS"]) -> object:
        """Parse genbank records into pandas DF/Bed format and dict format saving to self

        Args:
            feature_types (List[str]): a list of Genbank feature types to use

        Returns:
            None
        """
        feature_dict = {}
        pddict = {'chrom':[], 'chromStart':[], 'chromEnd':[], 'name':[],'score':[],'strand':[]}
        for gbfile in self.genbank_list:
            if gbfile.lower().endswith(".gz"):  # uses file extension
                f = gzip.open(gbfile, mode='rt')
            else:
                f = open(gbfile, mode='r')
            genebank_file = SeqIO.parse(f, "genbank")
            for entry in genebank_file:
                for record in entry.features:
                    if record.type in feature_types:
                        featid = hashlib.md5(str(record).encode()).hexdigest()
                        pddict['chrom'].append(entry.id)
                        pddict["chromStart"].append(record.location.start.position)
                        pddict["chromEnd"].append(record.location.end.position)
                        pddict["name"].append(featid)
                        pddict["score"].append(0)
                        pddict["strand"].append("-" if record.strand < 0 else "+")
                        for qualifier_key, qualifier_val in record.qualifiers.items():
                            if not qualifier_key in feature_dict:
                                feature_dict[qualifier_key] = {}
                            feature_dict[qualifier_key][featid] = qualifier_val
            genbankbed = pd.DataFrame.from_dict(pddict)
            self.genbank_bed_df = genbankbed
            self.feature_dict = feature_dict
        f.close()

    def _get_qualifiers(self, min_prop: float=0.5, excluded: List[str]= ["translation"]) -> Object:
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
        final_quals = []
        qual_df = pd.DataFrame()
        for quals in self.feature_dict:
            if quals/len(self.feature_dict) > min_prop
            final_quals.append(quals)
        for qualifier in final_quals:
            if qualifier not in excluded:
                featlist = []
                quallist = []
                for feat, qual inself.feature_dict[qualifier].items()
                    featlist.append(feat)
                    quallist.append(qual)
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
        # format to featurefile
        featurebed = BedTool.from_dataframe(self.genbank_bed_df)
        # format to mapping file
        mapbed = BedTool.from_dataframe(self.target_bed_df)
        # get feature downstream of target sequence
        downstream = mapbed.closest(featurebed, d=True, fd=True, D="a", t="first")
        # get feature upstream of target sequence
        upstream = mapbed.closest(featurebed, d=True, id=True, D="a", t="first")
        headers = {0: "Accession", 1: "Guide start", 2: "Guide end", 3:"Guide sequence",
                   4: "Score", 5:"Guide strand", 6: "Feature Accession",
                   7: "Feature start", 8:"Feature end",9:"Feature id",
                   10:"Feature score", 11:"Feature strand", 12: "Feature distance"}
        downstream = downstream.to_dataframe(disable_auto_names=True, header=None)
        downstream = downstream.rename(columns=headers)
        downstream = downstream['direction'] = 'downstream'
        upstream = upstream.to_dataframe(disable_auto_names=True, header=None)
        upstream = upstream.rename(columns=headers)
        upstream = upstream['direction'] = 'upstream'
        self.nearby = upstream.append(downstream)

def _filter_features(before_feat=100, after_feat=200):
    """merge targets with Feature list and filter for guides close enough to interact.

    Args:
        before_feat (int): The maximum distance before the start of a feature
                           measured from closest point to guide
        after_feat (int): The maximum distance after the start codon (into the gene)


    Returns:
        None
    """
    # Select guide forward, gene forward targets
    filtered_df= self.nearby[(self.nearby["direction"] == "downstream") &
                             (self.nearby["Guide strand"] == "+") &
                             (self.nearby["Feature strand"] == "+") & (
                             (0 < self.nearby["Feature distance"] < before_feat) |
                             (0 == self.nearby["Feature distance"] &
                             self.nearby["Guide start"] - self.nearby["Feature start"] < after_feat)
                             )]
    filtered_df.append(self.nearby[(self.nearby["Guide strand"] == "-") &
                      (self.nearby["Feature strand"] == "-") & (
                      (0 < self.nearby["Feature distance"] < before_feat ) |
                      (0 == self.nearby["Feature distance"] &
                      (self.nearby["Feature start"] - self.nearby["Guide start"] < after_feat))
                      ])
    filtered_df.append(self.nearby[(self.nearby["Guide strand"] == "-") &
                      (self.nearby["Feature strand"] == "+") & (
                      ( 0 < self.nearby["Feature distance"] < before_feat ) |
                      (0 == self.nearby["Feature distance"] &
                      self.nearby["Feature start"] - self.nearby["Guide start"] < after_feat))
                      ])
    filtered_df.append(self.nearby[(self.nearby["Guide strand"] == "+") &
                      (self.nearby["Feature strand"] == "-") & (
                      (0 < self.nearby["Feature distance"] < before_feat ) |
                      (0 == self.nearby["Feature distance"] &
                      self.nearby["Feature start"] - self.nearby["Guide start"] < after_feat))
                      ])
    self.filtered_df = filtered_df

def _format_guide_table(targetlist str:) -> None:
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
        return targetlist.neighbors[seq]["neighbors"].dist
    def get_off_target_seqs(seq):
        return targetlist.neighbors[seq]["neighbors"].seqs
    pretty_df = self.filtered_df.copy()
    pretty_df['GC'] = pretty_df['Guide sequence'].apply(gc)
    pretty_df['PAM'] = pretty_df['Guide sequence'].apply(get_exact_pam)
    pretty_df['Guide name'] = pretty_df['Guide sequence'].apply(get_guide_hash)
    pretty_df['Target strand'] = np.where(pretty_df['Guide strand']==pretty_df['Feature strand'], 'coding', 'non-coding')
    pretty_df['Off target scores'] = pretty_df['Guide sequence'].apply(get_off_target_score)
    pretty_df['Off target guides'] = pretty_df['Guide sequence'].apply(get_off_target_seqs)
    pretty_df = pretty_df[['Guide name',"Guide sequence", 'GC', "Accession","Guide start","Guide end",,
                "Guide strand", 'PAM', "Feature id","Feature start",
                "Feature end", "Feature strand",
                "Feature distance", 'Off target guides', 'Off target scores']]
    pretty_df = pretty_df.merge(self.qualifiers, how="left", on="Feature id")


def get_fastas(filelist, tempdir):
    """Saves a Fasta and from 1 or more Genbank files (may be gzipped)

    Args:
        filelist (str): Genbank file to process

    Returns:
        None

    """
    try:
        with open(os.path.join(tempdir, "forward.fasta"), "w") as f1:
            gblist = []
            for file in filelist:
                if file.lower().endswith(".gz"):  # uses file extension
                    f = gzip.open(file, mode='rt')
                else:
                    f = open(file, mode='r')
                    records = (SeqIO.parse(f, "genbank"))
                    SeqIO.write(records, f1, "fasta")
    except Exception as e:
        print("An error occurred in input genbank file %" % file)
        raise e
    finally:
        f.close()
        f1.close()
