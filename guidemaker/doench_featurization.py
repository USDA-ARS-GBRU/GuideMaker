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
from itertools import product
from time import time
from typing import List
import logging
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio.SeqUtils import MeltingTemp as Tm


def featurize_data(data, learn_options, pam_audit=True, length_audit=True, quiet=True) -> dict:
    """Creates a dictionary of feature data

    Args:
        data (np.array): of 30-mer sequences
        learn_options (dict): ict of model training parameters
        pam_audit (bool): should check of GG  at position 25:27 be performed?
        length_audit (bool): should sequence length be checked?
        quiet (bool):  Should progress be printed?

    Returns:
        (dict): Returns a dict containing feature types as keys and arrays as values
    """

    print("start features")
    if np.any(data["30mer"].str.len() != 30):
        raise AssertionError(f"should only have sequences 30 nt long")

    if not quiet:
        print("Constructing features...")
    time_0 = time()

    feature_sets = {}

    if learn_options["nuc_features"]:
        # spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        get_all_order_nuc_features(
            data["30mer"],
            feature_sets,
            learn_options,
            learn_options["order"],
            max_index_to_use=30,
            quiet=quiet,
        )
    print("nuc_ceatures complete")
    check_feature_set(feature_sets)
    print("check_feature_set complete")
    if learn_options["gc_features"]:
        gc_above_10, gc_below_10, gc_count = gc_features(data, length_audit)
        feature_sets["gc_above_10"] = pd.DataFrame(gc_above_10)
        feature_sets["gc_below_10"] = pd.DataFrame(gc_below_10)
        feature_sets["gc_count"] = pd.DataFrame(gc_count)
    print("gc features complete")

    if learn_options["include_NGGX_interaction"]:
        feature_sets["NGGX"] = nggx_interaction_feature(data, pam_audit)
    print("ggx features complte")
    if learn_options["include_Tm"]:
        feature_sets["Tm"] = Tm_feature(data, pam_audit, learn_options=None)
    print("Tm complete")

    time_1 = time()
    if not quiet:
        print(
            f"\t\tElapsed time for constructing features is {time_1 - time_0:.2f} seconds"
        )

    check_feature_set(feature_sets)
    print("final feature check complte")

    return feature_sets

def parallel_featurize_data(data, learn_options, pam_audit=True, length_audit=True, quiet=True, num_threads=1) -> dict:
    if num_threads > 1:
        dflist = np.array_split(data, num_threads)
        partial_fd = partial(featurize_data,learn_options=learn_options, pam_audit=pam_audit, length_audit=length_audit, quiet=quiet )
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
        return featurize_data(data=data, learn_options=learn_options, pam_audit=pam_audit, length_audit=length_audit, quiet=quiet)


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
                    "# of individuals do not match up across feature sets"
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


def get_all_order_nuc_features(
    data,
    feature_sets,
    learn_options,
    maxorder,
    max_index_to_use,
    prefix="",
    quiet=False,
):
    for order in range(1, maxorder + 1):
        if not quiet:
            print(f"\t\tconstructing order {order} features")
        nuc_features_pd, nuc_features_pi = apply_nucleotide_features(
            data,
            order,
            include_pos_independent=True,
            max_index_to_use=max_index_to_use,
            prefix=prefix,
        )
        feature_sets[f"{prefix}_nuc_pd_Order{order:d}"] = nuc_features_pd
        if learn_options["include_pi_nuc_feat"]:
            feature_sets[f"{prefix}_nuc_pi_Order{order:d}"] = nuc_features_pi
        check_feature_set(feature_sets)

        if not quiet:
            print("\t\t\t\t\t\t\tdone")


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


def apply_nucleotide_features(
    seq_data_frame, order, include_pos_independent, max_index_to_use, prefix=""
):
    if include_pos_independent:
        feat_pd = seq_data_frame.apply(
            nucleotide_features, args=(order, max_index_to_use, prefix, "pos_dependent")
        )
        feat_pi = seq_data_frame.apply(
            nucleotide_features,
            args=(order, max_index_to_use, prefix, "pos_independent"),
        )
        if np.any(np.isnan(feat_pd)):
            raise AssertionError(
                "nans here can arise from sequences of different lengths"
            )
        if np.any(np.isnan(feat_pi)):
            raise AssertionError(
                "nans here can arise from sequences of different lengths"
            )
    else:
        feat_pd = seq_data_frame.apply(
            nucleotide_features, args=(order, max_index_to_use, prefix, "pos_dependent")
        )
        if np.any(np.isnan(feat_pd)):
            raise AssertionError("found nan in feat_pd")
        feat_pi = None
    return feat_pd, feat_pi


def get_alphabet(order, raw_alphabet=("A", "T", "C", "G")):
    alphabet = ["".join(i) for i in product(raw_alphabet, repeat=order)]
    return alphabet


def nucleotide_features(
    s,
    order,
    max_index_to_use,
    prefix="",
    feature_type="all",
    raw_alphabet=("A", "T", "C", "G"),
):
    """
    compute position-specific order-mer features for the 4-letter alphabet
    (e.g. for a sequence of length 30, there are 30*4 single nucleotide features
          and (30-1)*4^2=464 double nucleotide features
    """
    if feature_type not in ["all", "pos_independent", "pos_dependent"]:
        raise AssertionError("Unknown feature_type")
    if max_index_to_use <= len(s):
        max_index_to_use = len(s)

    if max_index_to_use is not None:
        s = s[:max_index_to_use]
    # s = s[:30] #cut-off at thirty to clean up extra data that they accidentally left in,
    # nd were instructed to ignore in this way
    alphabet: List[str] = get_alphabet(order, raw_alphabet=raw_alphabet)
    features_pos_dependent = np.zeros(len(alphabet) * (len(s) - (order - 1)))
    features_pos_independent = np.zeros(np.power(len(raw_alphabet), order))

    index_dependent: List[str] = []
    index_independent: List[str] = []

    for position in range(0, len(s) - order + 1, 1):
        for let in alphabet:
            index_dependent.append(f"{prefix}{let}_{position:d}")

    for let in alphabet:
        index_independent.append(f"{prefix}{let}")

    for position in range(0, len(s) - order + 1, 1):
        nucl: object = s[position: position + order]
        features_pos_dependent[alphabet.index(nucl) + (position * len(alphabet))] = 1.0
        features_pos_independent[alphabet.index(nucl)] += 1.0

        # this is to check that the labels in the pd df actually match the nucl and position
        if (
            index_dependent[alphabet.index(nucl) + (position * len(alphabet))]
            != f"{prefix}{nucl}_{position:d}"
        ):
            raise AssertionError()
        if index_independent[alphabet.index(nucl)] != f"{prefix}{nucl}":
            raise AssertionError()

    if np.any(np.isnan(features_pos_dependent)):
        raise Exception("found nan features in features_pos_dependent")
    if np.any(np.isnan(features_pos_independent)):
        raise Exception("found nan features in features_pos_independent")

    if feature_type in ("all", "pos_independent"):
        if feature_type == "all":
            res = pd.Series(features_pos_dependent, index=index_dependent).append(
                pd.Series(features_pos_independent, index=index_independent)
            )
            if np.any(np.isnan(res.values)):
                raise AssertionError()
        else:
            res = pd.Series(features_pos_independent, index=index_independent)
            if np.any(np.isnan(res.values)):
                raise AssertionError()
    else:
        res = pd.Series(features_pos_dependent, index=index_dependent)

    if np.any(np.isnan(res.values)):
        raise AssertionError()
    return res


def nucleotide_features_dictionary(prefix=""):
    seqname = ["-4", "-3", "-2", "-1"]
    seqname.extend([str(i) for i in range(1, 21)])
    seqname.extend(["N", "G", "G", "+1", "+2", "+3"])

    orders = [1, 2, 3]
    sequence = 30
    feature_names_dep = []
    feature_names_indep = []
    index_dependent = []
    index_independent = []

    for order in orders:
        raw_alphabet = ["A", "T", "C", "G"]
        alphabet = ["".join(i) for i in product(raw_alphabet, repeat=order)]
        features_pos_dependent = np.zeros(len(alphabet) * (sequence - (order - 1)))
        features_pos_independent = np.zeros(np.power(len(raw_alphabet), order))

        index_dependent.extend(
            [
                f"{prefix}_pd.Order{order}_P{i}"
                for i in range(len(features_pos_dependent))
            ]
        )
        index_independent.extend(
            [
                f"{prefix}_pi.Order{order}_P{i}"
                for i in range(len(features_pos_independent))
            ]
        )

        for pos in range(sequence - (order - 1)):
            for letter in alphabet:
                feature_names_dep.append(f"{letter}_{seqname[pos]}")

        for letter in alphabet:
            feature_names_indep.append(f"{letter}")

        if len(feature_names_indep) != len(index_independent):
            raise AssertionError()
        if len(feature_names_dep) != len(index_dependent):
            raise AssertionError()

    index_all = index_dependent + index_independent
    feature_all = feature_names_dep + feature_names_indep

    return dict(zip(index_all, feature_all))
