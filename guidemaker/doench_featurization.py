"""doench_featurization.py. Simplified feature extraction to run the model 'V3_model_nopos' from Doench et al. 2016 for
   on-target scoring.

For use in Guidemaker https://guidemaker.org.
Adam Rivers, Unites States Department of Agriculture, Agricultural Research Service

Core code https://github.com/MicrosoftResearch/Azimuth is in Python2 and does not run well given changes to packages.
Miles Smith worked on porting to Python3 in this repo: https://github.com/milescsmith/Azimuth. including a new branch
that used Poetry to build. The work is not complete.

This work is derivitive of that BSD 3-clause licensed work. The key changes are:
 1. uch of the code needed for tasks other thant prediction of the V3_model_nopos was removed.
 2. The Calculation of NGGX features was re-written. A bug that prevented scaling to thousands guides efficiently.
 3. the Pickle model and scikit-learn were replaced with an Onnx model and onnxruntime for better persistance,
    security, and performance.

Reference:

John G. Doench*, Nicolo Fusi*, Meagan Sullender*, Mudra Hegde*, Emma W. Vaimberg*, Katherine F. Donovan, Ian Smith,
Zuzana Tothova, Craig Wilen , Robert Orchard , Herbert W. Virgin, Jennifer Listgarten*, David E. Root.
Optimized sgRNA design to maximize activity and minimize off-target effects for genetic screens with CRISPR-Cas9.
Nature Biotechnology Jan 2016, doi:10.1038/nbt.3437.
"""
from itertools import product, islice
from time import time
from typing import List
import logging
from functools import partial
from multiprocessing import Pool
from collections import deque

import numpy as np
import pandas as pd
from Bio.SeqUtils import MeltingTemp as Tm

logger = logging.getLogger(__name__)

def featurize_data(data, learn_options, pam_audit=True, length_audit=True) -> dict:
    """Creates a dictionary of feature data

    Args:
        data pd.DataFrame: of 30-mer sequences in column 1 and strand in column 2
        learn_options (dict): dict of model training parameters
        pam_audit (bool): should a check of GG  at position 25:27 be performed?
        length_audit (bool): should sequence length be checked?

    Returns:
        (dict): Returns a dict containing pandas dataframs of features
    """

    logging.info("Creating features for Doench et al. 2016 score prediction")
    if np.any(data["30mer"].str.len() != 30):
        raise AssertionError(f"Sequences should be 30 nt long")


    feature_sets = {}

    if learn_options["nuc_features"]:
        # spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        logging.info("Creating nucleotide features")
        feature_sets["_nuc_pd_Order1"], feature_sets["_nuc_pi_Order1"], feature_sets["_nuc_pd_Order2"], feature_sets["_nuc_pi_Order2"] = get_nuc_features(data)

    logging.info("Verifing nucleotide features")
    check_feature_set(feature_sets)

    if learn_options["gc_features"]:
        logging.info("Creating GC features")
        gc_above_10, gc_below_10, gc_count = gc_features(data, length_audit)
        feature_sets["gc_above_10"] = pd.DataFrame(gc_above_10)
        feature_sets["gc_below_10"] = pd.DataFrame(gc_below_10)
        feature_sets["gc_count"] = pd.DataFrame(gc_count)
        logging.info("gc features complete")

    if learn_options["include_NGGX_interaction"]:
        logging.info("Creating ggx features")
        feature_sets["NGGX"] = nggx_interaction_feature(data, pam_audit)

    if learn_options["include_Tm"]:
        logging.info("Creattng Tm features")
        feature_sets["Tm"] = Tm_feature(data, pam_audit, learn_options=None)


    check_feature_set(feature_sets)
    logging.info("final feature check complte")

    return feature_sets

def parallel_featurize_data(data, learn_options, pam_audit=True, length_audit=True, num_threads=1) -> dict:
    """ Use multprocessing to divide up the creation of ML features for Doench scoring
        Creates a dictionary of feature data

        Args:
            data pd.DataFrame:      of 30-mer sequences in column 1 and strand in column 2
            learn_options (dict):   dict of model training parameters
            pam_audit (bool):       should a check of GG  at position 25:27 be performed?
            length_audit (bool):    should sequence length be checked?

        Returns:
            (dict): Returns a dict containing pandas dataframs of features

    """
    if num_threads > 1:
        dflist = np.array_split(data, num_threads)
        partial_fd = partial(featurize_data,learn_options=learn_options, pam_audit=pam_audit, length_audit=length_audit )
        with Pool(processes=num_threads) as pool:
            result = pool.map(partial_fd, dflist)
        featdict = dict.fromkeys(result[0].keys())
        for featkey in featdict.keys():
            tempdflist = []
            for d1 in result:
                tempdflist.append(d1[featkey])
            featdict[featkey] = pd.concat(tempdflist)
        return featdict
    else:
        return featurize_data(data=data, learn_options=learn_options, pam_audit=pam_audit, length_audit=length_audit)


def get_nuc_features(data):
    """ Create first and second order nucleotide features

        Args:
            data pd.DataFrame:      of 30-mer sequences in column 1 and strand in column 2

        Returns:
            (tuple): Returns a tuple pwith 4 Pandas dataframes

    """
    seqlen = 30
    # create first header
    nuc_pi_Order1_header = [x[0] for x  in product('ATCG', repeat=1)]

    # create second header
    nuc_pd_Order1_header = []
    for i in range(seqlen):
        nuc_pd_Order1_header.extend([x + "_" + str(i) for x in  nuc_pi_Order1_header ])

    # create third header
    nuc_pi_Order2_header = [ "".join(x) for x in product('ATCG', repeat=2)]

    # create forth header
    nuc_pd_Order2_header = []
    for i in range(seqlen - 1):
        nuc_pd_Order2_header.extend([x + "_" + str(i) for x in  nuc_pi_Order2_header ])

    # Create lists for holding features
    nuc_pd_Order1_list = []
    nuc_pi_Order1_list = []
    nuc_pd_Order2_list = []
    nuc_pi_Order2_list = []

    # Create identity matricies and lookups for one hot encoding
    o1_id = np.eye(4)
    o2_id = np.eye(16)
    o1_lookup = {x : i for i, x in enumerate(nuc_pi_Order1_header)}
    o2_lookup = {x : i for i, x in enumerate(nuc_pi_Order2_header)}

    # not feasable ot use pandas.get_dummies for this
    def one_hot(seq, idmat, lookup):
        """One hot endode  an iterable matching the items in a lookup

        Args:
            seq (str): a 30mer sequence string
            idmat (np.array): a Numpy identity matrix
            lookup (dict): a dictionary matchin the substring to the row in idmat

        Returns:
            pd.Series: one hot encoding
        """
        featurevect = []
        for let in seq:
            pos  = lookup[let]
            featurevect.extend(list(idmat[pos,:]))
        return pd.Series(featurevect)

    def sliding_window(iterable, n):
        """Create a generator of substrings

            Args:
                iterable (iterable): an itterable object like a string or list
                n (int): the size of the window or kmer

            Returns:
                a genrator of substrings
        """
        # sliding_window('ABCDEFG', 4) -> ABCD BCDE CDEF DEFG
        it = iter(iterable)
        window = deque(islice(it, n), maxlen=n)
        if len(window) == n:
            yield tuple(window)
        for x in it:
            window.append(x)
            yield tuple(window)

    for seq in data["30mer"]:
        # add order 1 frequency features
        pi1dict = dict.fromkeys(nuc_pi_Order1_header, 0)
        for let in seq:
            pi1dict[let] += 1
        nuc_pi_Order1_list.append(pi1dict)
        # add order 1 positon features
        nuc_pd_Order1_list.append(one_hot(seq, o1_id, o1_lookup))
        # create list of 2mers
        seq_2mers =  [ "".join(x) for x in list(sliding_window(seq, 2))]
        # add order two frequency features
        pi2dict = dict.fromkeys(nuc_pi_Order2_header, 0)
        for let in seq_2mers:
            pi2dict[let] += 1
        nuc_pi_Order2_list.append(pi2dict)
        # add order 2 positon features
        nuc_pd_Order2_list.append(one_hot(seq_2mers, o2_id, o2_lookup))

    # Create DataFrames
    nuc_pd_Order1 = pd.DataFrame(data=nuc_pd_Order1_list)
    nuc_pd_Order1.columns = nuc_pd_Order1_header
    nuc_pi_Order1 = pd.DataFrame(data=nuc_pi_Order1_list, columns=nuc_pi_Order1_header)
    nuc_pd_Order2 = pd.DataFrame(data=nuc_pd_Order2_list)
    nuc_pd_Order2.columns = nuc_pd_Order2_header
    nuc_pi_Order2 = pd.DataFrame(data=nuc_pi_Order2_list, columns=nuc_pi_Order2_header)

    return nuc_pd_Order1, nuc_pi_Order1, nuc_pd_Order2, nuc_pi_Order2



def check_feature_set(feature_sets):
    """Ensure the number of features is the same in each feature set

    Args:
            feature_sets (dict): the feature set dictionary
    Returns:
        None
    """
    if feature_sets == {}:
        raise AssertionError("no feature sets present")

    num = None
    for ft in feature_sets:
        num2 = feature_sets[ft].shape[0]
        if num is None:
            num = num2
        else:
            if num < 1:
                raise AssertionError("should be at least one individual")
            if num != num2:
                raise AssertionError(
                    "Number of individuals do not match up across feature sets"
                )

    for item in feature_sets:
        if np.any(np.isnan(feature_sets[item])):
            raise Exception(f"found Nan in set {item}")


def nggx_interaction_feature(data, pam_audit=True):
    """ One hot encode the sequence of NX aroung pam site NGGX

    Args:
        data (pandas.DataFrame): A dataframe of 30-mer and strand (filled with NA)
        pam_audit (bool): should check of GG  at position 25:27 be performed?

    Returns:
        (pandas.DataFrame):  A dataframe with 16 columns containing  one hoe encoding of NX data
    """
    # function completly replaced from Doench et al. (2016). The old function had time complexity of O^2
    sequence = data["30mer"].values
    nxcombos = []
    for i in product('ACGT', repeat=2):
        nxcombos.append("".join(list(i)))
    nxlist = []
    # check that GG is where we think
    for seq in sequence:
        if pam_audit and seq[25:27] != "GG":
            raise Exception(f"expected GG but found {seq[25 :27]}")
        nxlist.append(seq[24] + seq[27])
    nxs = pd.Series(nxlist)
    feat_nx = pd.get_dummies(nxs, columns=nxcombos)
    # Add any missing columns
    for i, header in enumerate(nxcombos):
        if header not in feat_nx.columns:
            feat_nx.insert(i, header, 0)
    feat_nx.columns = ["NGGX" + i + "_0" for i in nxcombos]
    feat_nx = feat_nx.astype('float64')
    return feat_nx




def countGC(s, length_audit=True):
    """
    GC content for only the 20mer, as per the Doench paper/code
    """
    if length_audit:
        if len(s) != 30:
            raise AssertionError("seems to assume 30mer")
    return len(s[4:24].replace("A", "").replace("T", ""))


def organism_feature(data):
    """
    Human vs. mouse
    """
    organism = np.array(data["Organism"].values)
    feat = pd.DataFrame(pd.DataFrame(organism))
    import pdb

    pdb.set_trace()
    return feat


def gc_cont(seq):
    return (seq.count("G") + seq.count("C")) / float(len(seq))


def Tm_feature(data, pam_audit=True, learn_options=None):
    """
    assuming '30-mer'is a key
    get melting temperature features from:
        0-the 30-mer ("global Tm")
        1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA,
        i.e. the 5nts immediately proximal of the NGG PAM
        2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
        3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
    """

    if learn_options is None or "Tm segments" not in learn_options:
        segments = [(19, 24), (11, 19), (6, 11)]
    else:
        segments = learn_options["Tm segments"]

    sequence = data["30mer"].values
    featarray = np.ones((sequence.shape[0], 4))

    for i, seq in enumerate(sequence):
        if pam_audit and seq[25:27] != "GG":
            raise Exception(f"expected GG but found {seq[25:27]}")
        rna = False
        featarray[i, 0] = Tm.Tm_NN(seq, nn_table=Tm.RNA_NN2)  # 30mer Tm
        featarray[i, 1] = Tm.Tm_NN(
            seq[segments[0][0]: segments[0][1]], nn_table=Tm.RNA_NN2
        )  # 5nts immediately proximal of the NGG PAM
        featarray[i, 2] = Tm.Tm_NN(
            seq[segments[1][0]: segments[1][1]], nn_table=Tm.RNA_NN2
        )  # 8-mer
        featarray[i, 3] = Tm.Tm_NN(
            seq[segments[2][0]: segments[2][1]], nn_table=Tm.RNA_NN2
        )  # 5-mer

    feat = pd.DataFrame(
        featarray,
        index=data.index,
        columns=[
            f"Tm global_{rna}",
            f"5mer_end_{rna}",
            f"8mer_middle_{rna}",
            f"5mer_start_{rna}",
        ],
    )

    return feat


def gc_features(data, audit=True):
    gc_count = data["30mer"].apply(lambda seq: countGC(seq, audit))
    gc_count.name = "GC count"
    gc_above_10 = (gc_count > 10) * 1
    gc_above_10.name = "GC > 10"
    gc_below_10 = (gc_count < 10) * 1
    gc_below_10.name = "GC < 10"
    return gc_above_10, gc_below_10, gc_count


def normalize_features(data, axis):
    """
    input: pd.DataFrame of dtype=np.float64 array, of dimensions
    mean-center, and unit variance each feature
    """
    data -= data.mean(axis)
    data /= data.std(axis)
    # remove rows with NaNs
    data = data.dropna(1)
    if np.any(np.isnan(data.values)):
        raise Exception("found NaN in normalized features")
    return data
