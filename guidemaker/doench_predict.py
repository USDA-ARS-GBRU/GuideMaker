"""doench_predict.py. Simplified module to run the model 'V3_model_nopos' from Doench et al. 2016 for on-target
   scoring.

For use in Guidemaker https://guidemaker.org.
Adam Rivers, Unites States Department of Agriculture, Agricultural Research Service

The core code, https://github.com/MicrosoftResearch/Azimuth, is in Python2 and does not run well given changes to packages.
Miles Smith worked on porting to Python3 in this repo: https://github.com/milescsmith/Azimuth, including a new branch
that used Poetry to build. The work is not complete.

This work is derivative of that BSD 3-clause, Modified licensed work. The key changes are:
 1. Much of the code needed for tasks other thant prediction of the V3_model_nopos was removed.
 2. The Calculation of NGGX features was re-written. A bug that prevented scaling to thousands guides efficiently.
 3. the Pickle model and scikit-learn were replaced with an Onnx model ('V3_model_nopos.onnx"), and medadata file
    ("V3_model_nopos_options.json") and onnxruntime for better persistence, security, and performance.

Reference:

John G. Doench*, Nicolo Fusi*, Meagan Sullender*, Mudra Hegde*, Emma W. Vaimberg*, Katherine F. Donovan, Ian Smith,
Zuzana Tothova, Craig Wilen , Robert Orchard , Herbert W. Virgin, Jennifer Listgarten*, David E. Root.
Optimized sgRNA design to maximize activity and minimize off-target effects for genetic screens with CRISPR-Cas9.
Nature Biotechnology Jan 2016, doi:10.1038/nbt.3437.
"""

import pandas as pd
from guidemaker.doench_featurization import featurize_data, parallel_featurize_data
from typing import List, Optional
import numpy as np
from pathlib import Path
import onnxruntime as rt
import json
import logging
import os

logger = logging.getLogger(__name__)

DIR = os.path.dirname(os.path.abspath(__file__))  # This is your Project Root


#MODEL = "guidemaker/data/V3_model_nopos.onnx"
MODEL = os.path.join(DIR,"data/V3_model_nopos.onnx")
#MODEL_META = "guidemaker/data/V3_model_nopos_options.json"
MODEL_META = os.path.join(DIR,"data/V3_model_nopos_options.json")

def concatenate_feature_sets(feature_sets: dict, keys: List[str] = None) -> tuple:
    """ Combine features

    Given a dictionary of sets of features, each in a pd.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set

    Args:
        feature_sets (dict): a Dict of feature sets as pandas DataFrames

    Returns:
        (tuple: inputs(numpy.ndarray), dim (tuple)
    """
    if feature_sets == {}:
        raise AssertionError("no feature sets present")
    if keys is None:
        keys = list(feature_sets.keys())
    feature_1 = feature_sets[keys[0]].shape[0]
    for assemblage in feature_sets:
        feature_2 = feature_sets[assemblage].shape[0]
        if feature_1 != feature_2:
            raise AssertionError(
                f"not same # items for features {keys[0]} and {assemblage}"
            )
    num_sets = feature_sets[keys[0]].shape[0]
    inputs = np.zeros((num_sets, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for assemblage in keys:
        inputs_set = feature_sets[assemblage].values
        dim[assemblage] = inputs_set.shape[1]
        dimsum = dimsum + dim[assemblage]
        inputs = np.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[assemblage].columns.tolist())
    return inputs, dim, dimsum, feature_names


def predict(
    seq: np.ndarray,
    model_file: Optional[Path] = MODEL,
    model_metadata: Optional[Path] = MODEL_META,
    pam_audit: bool = True,
    length_audit: bool = False,
    num_threads: int = 1
) -> np.array:
    """Pedicts regressions scored from sequences.

    Args:
        seq (numpy.ndarray) numpy array of 30 nt sequences with 25 nt of guide, NGG pam in 25:27 and the following 2 nts.
        model_file (str): file path of the onnx Boosted Gradien Regressor model file without position data
        model_metadata (str): file path of the json model parameters metadata file.
        pam_audit (bool): check PAM of each sequence.
        length_audit(bool) : check length of each sequence.

    Returns:
        (numpy.array): An array with regression values.

     """
    if not isinstance(seq, np.ndarray):
        raise AssertionError("Please ensure seq is a numpy array")
    if len(seq[0]) <= 0:
        raise AssertionError("Make sure that seq is not empty")
    if not isinstance(seq[0], str):
        raise AssertionError(
            f"Please ensure input sequences are in string format, i.e. 'AGAG' "
            f"rather than ['A' 'G' 'A' 'G'] or alternate representations"
        )
    # Open onnx runtime session
    sess = rt.InferenceSession(model_file)
    with open(model_metadata, "r") as f:
        learn_options = json.load(f)
    x_df = pd.DataFrame(
        columns=["30mer", "Strand"],
        data=list(zip(seq, ["NA" for x in range(len(seq))])),
    )

    feature_sets = parallel_featurize_data(
        x_df,
        learn_options,
        pam_audit=pam_audit,
        length_audit=length_audit,
        num_threads=num_threads
    )
    inputs, *_ = concatenate_feature_sets(feature_sets)
    preds = sess.run(None, {'float_input': inputs.astype(np.float32)})[0]
    return preds
