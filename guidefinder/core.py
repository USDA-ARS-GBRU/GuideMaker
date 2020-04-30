"""Core classes and functions for Guidefinder

"""
import os
from typing import List, Set, Dict, Tuple
from itertools import product, tee, chain
import gzip
from Bio.Seq import Seq
from Bio import SeqIO
import nmslib
from pybedtools import BedTool
import pandas as pd




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
            seq_obj (object): A Biopython SeqRecord instance
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
            elif hitset == "fset" and strand in ("forward", "both") and self.pam_orientation == "5prime":
                start = i + len(self.pam)
                stop = i + len(self.pam)  + target_len
                seq = str(seqrecord_obj[start:stop].seq)
            elif hitset == "rset" and strand in ("reverse", "both") and self.pam_orientation == "3prime":
                start = i + len(self.pam)
                stop = i + len(self.pam) + target_len
                seq = str(seqrecord_obj[start:stop].seq.reverse_complement())
            elif hitset == "rset" and strand in ("reverse", "both") and self.pam_orientation == "5prime":
                start = i - target_len
                stop = i
                seq = str(seqrecord_obj[start:stop].seq.reverse_complement())
            else:
                return None
            if 0 <= start <= len(seqrecord_obj) and  0 <= stop <= len(seqrecord_obj):
                return Target(seq=seq,
                              pam=self.pam,
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
    def __init__(self, seq: str, pam: str, strand: str, pam_orientation: str,
                 seqid: str, start: int, stop: int) -> None:
        self.seq: str= seq
        self.pam: str = pam
        self.strand: str = strand
        self.pam_orientation: str = pam_orientation
        self.seqid: str = seqid
        self.start: int = start
        self.stop: int = stop

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
                                     "strand": strand})
        df.sort_values(by=['chrom','chromstart'], inplace=True)
        # df = df.replace(np.nan, "NA")
        return df
        #df.to_csv(path=outfile, sep="\t", index=False, header=False)



def get_genbank_features(genbank_list: List[str]) -> object:
    """Return a list of features for a genbank file

    Args:
        filelist(genebank): Genbank file to process


    Returns:
        (obj): A Pandas Dataframe in Bed format
    """
    feature_list = []
    for gbfile in genbank_list:
    for file in filelist:
        if file.lower().endswith(".gz"):  # uses file extension
            f = gzip.open(file, mode='rt')
        else:
            f = open(file, mode='r')
        genebank_file = SeqIO.parse(f, "genbank")
        for entry in genebank_file:
            for record in entry.features:
                feature_dict = {}
                if record.type in ['CDS', 'gene']:
                    feature_dict["accession"] = entry.id
                    feature_dict["start"] = record.location.start.position
                    feature_dict["stop"] = record.location.end.position
                    feature_dict["type"] = record.type
                    feature_dict["strand_for_feature"] = 'reverse' if record.strand < 0 else 'forward'
                    for qualifier_key, qualifier_val in record.qualifiers.items():
                        feature_dict[qualifier_key] = qualifier_val
                    feature_list.append(feature_dict)
        genebankfeatures_df = pd.DataFrame(feature_list)
        #genebankfeatures_df = genebankfeatures_df.replace(np.nan, "NA")
        return genebankfeatures_df
        #genebankfeatures_df.to_csv(path=outfile, index=False, sep='\t',header=False)
    f.close()

def get_nearby_feature(targets: object, features: object) -> Tuple[object, object]:
    """Adds downstream information to the given target sequences and mapping information

    Args:
        targets (obj) : A Pandas dataframe in Bed format with targets
        features (lobj) : A Pandas dataframe in Bed format with GenBank information

    Returns:
        Tuple(object, object): Return a tuple with two BedTool objects
        containing target sequences, mapping information, and adjacent gene information
    """
    # format to featurefile
    featurebed = BedTool.from_dataframe(features)
    # format to mapping file
    mapbed = BedTool.from_dataframe(targets)
    # get feature downstream of target sequence
    downstream = mapbed.closest(featurebed, d=True, fd=True, D="a", t="first")
    # get feature upstream of target sequence
    upstream = mapbed.closest(featurebed, d=True, id=True, D="a", t="first")
    return downstream, upstream


def merge_downstream_upstream(downsfile, upsfile, columns_name, outputfilename):
    """Return a merged file

    Args:
        (tsv?? check type of pybed output): A file with target sequences, mapping information, and upstream information
        (tsv?? check type of pybed output): A file with target sequences, mapping information, and downstream information

    Returns:
        dataframe: A DataFrame with merged upstream information and downstream information for a target sequence
    """
    downstream_df = downsfile.to_dataframe(names=columns_name,low_memory=False)
    upstream_df = upsfile.to_dataframe(names=columns_name,low_memory=False)
    all_df = pd.merge(downstream_df, upstream_df,
                      right_on=columns_name[:8],
                      left_on=columns_name[:8],
                      how='outer')
    all_df.to_csv(outputfilename, index=False, sep='\t', header=True)
    return outputfilename

def get_fastas(filelist, tempdir):
    """Returns Fasta and from a Genbank file

    Args:
        filelist (str): Genbank file to process

    Returns:
        forward.fasta(str): Fasta file containing DNA sequences in the genome

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
        print("An error occurred in input genbank file")
        raise e
    finally:
        f.close()
        f1.close()
