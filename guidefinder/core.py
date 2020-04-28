"""Core classes and functions for Guidefinder

"""
from typing import List, Set, Dict, Tuple
from itertools import product, tee, chain

from Bio import Seq
import nmslib
from pybedtools import BedTool
import pandas as pd



class Pam:
    """A Class representing a Protospacer Adjacent Motif (PAM)
    """
    def __init__(self, pam: str, pam_orientation: str) -> None:
        self.pam: str = str(Seq.Seq(pam).seq)
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
        return Pam(pam=str(Seq.Seq(self.pam).reverse_complement().seq),
                   pam_orientation=self.pam_orientation)

    def extend_ambiguous_dna(self) -> Set[str]:
        """return list of all possible sequences given an ambiguous DNA input"""
        dnaval = Seq.IUPAC.IUPACData.ambiguous_dna_values
        dnalist = []
        for i in product(*[dnaval[j] for j in self.pam]):
            dnalist.append("".join(i))
        return frozenset(dnalist)

    def find_targets(self, seq_obj: object, strand: str, target_len: int) -> List[object]:
        """Find all targets on a sequence that match for the PAM on the requested strand(s)

        Args:
            seq_obj (object): A Biopython Seq instance
            strand (str): The strand to search choices: ["forward", "reverse", "both"]
            target_len (int): The length of the target sequence
        Returns:
            list: A list of Target class instances

        """
        if strand == "forward" or "both":
            fset = self.extend_ambiguous_dna()
        if strand == "reverse" or "both":
            rset = self.reverse_complement().extend_ambiguous_dna()
        target_list = []

        def window(iterable, size):
            iters = tee(iterable, size)
            for i in range(1, size):
                for each in iters[i:]:
                    next(each, None)
            return zip(*iters)

        def compute_seq_coords(self, hitset, strand, target_len, seq_obj, i):
            def tlappend(self, start, stop, seq):
                if 0 < start < len(seq_obj) and  0 < stop < len(seq_obj):
                    target_list.append(Target(seq=seq,
                                              pam=self.pam,
                                              strand=strand,
                                              pam_orientation=self.pam_orientation,
                                              seqid=seq_obj.id,
                                              start=start,
                                              stop=stop))
            if hitset == "fset" and strand in ("forward", "both") \
            and self.pam_orientation == "3prime":
                start = i - target_len - 1
                stop = i - 1
                seq = str(seq_obj[start:stop].seq)
                tlappend(self, start, stop, seq)
            elif hitset == "fset" and strand in ("forward", "both") \
            and self.pam_orientation == "5prime":
                start = i + len(self.pam) - 1
                stop = i + len(self.pam) + target_len
                seq = str(seq_obj[start:stop].seq)
                tlappend(self, start, stop, seq)
            elif hitset == "rset" and strand in ("reverse", "both") \
            and self.pam_orientation == "3prime":
                start = i + len(self.pam) - 1
                stop = i + len(self.pam) + target_len
                seq = str(seq_obj[start:stop].reverse_complement().seq)
                tlappend(self, start, stop, seq)
            elif hitset == "rset" and strand in ("reverse", "both") \
            and self.pam_orientation == "5prime":
                start = i - target_len - 1
                stop = i - 1
                seq = str(seq_obj[start:stop].reverse_complement().seq)
                tlappend(self, start, stop, seq)
        # Create iterable of PAM length windows across the sequence
        kmer_iter = window(seq_obj.seq, len(self.pam))
        for i, kmer in enumerate(kmer_iter):
            if ''.join(kmer) in fset:
                compute_seq_coords(self, hitset="fset",
                                   strand=strand,
                                   target_len=target_len,
                                   seq_obj=seq_obj,
                                   i=i)
            elif ''.join(kmer) in rset:
                compute_seq_coords(self, hitset="rset",
                                   strand=strand,
                                   target_len=target_len,
                                   seq_obj=seq_obj,
                                   i=i)
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
                distance of \u_2265 {}".format(len(self.targets),
                                              len(self.unique_targets),
                                              len(self.neighbors),
                                              self.levindist)
        return info

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
                final_dict[target.seq]=target
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
        index.addDataPointBatch(self.unique_targets.keys())
        index.createIndex(print_progress=True)
        self.nmslib_index = index

    def get_neighbors(self) -> None:
        """Get nearest neighbors for sequences.

        For the unique_sequences calculate the (knum) nearest neighbors.
        Writes a dictionary to self.neighbors:
        self.neighbors[seq]{target: seq_obj, neighbors: {seqs:[s1, s1, ...], dist:[d1, d1,...]}}

        Args: None
        Returns: None
        """
        seqs = self.unique_targets.keys()
        results_list = ref_index.knnQueryBatch(seqs,
                                               k=self.kmum)
        neighbor_dict = {}
        for entry, i in enumerate(results_list):
            if entry[1][1] <= self.levindist:
                seq = seqs[ i - 1]
                neighbors = {"seqs" : [seqs[x] for x in entry[0]], "dist": list(entry[1])}
                neighbor_dict[seqs[ i - 1]] = {"neighbors":neighbors}
            for target in self.unique_targets:
                if seq == target:
                    neighbor_dict[seqs[ i - 1]]["target"] = target

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
        #df = df.replace(np.nan, "NA")
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
    for file in filelist:
        genebank_file = SeqIO.parse(file,"genbank")
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
    mapbed = BedTool.from_dataframe(target_mapping)
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
        forward.fasta(str): Fasta file in forward orientation (5'-3')
        reverse_complement.fasta(str): Reverse Complement of Fasta file
    """
    try:
        with open(os.path.join(tempdir,"forward.fasta"), "w") as f1:
            gblist = []
            for file in filelist:
                gblist.append(SeqIO.parse(file, "genbank"))
            SeqIO.write(chain(gblist), f1, "fasta")
    except Exception as e:
        print("An error occurred in input genbank file")
        raise e
