"""cfd_score_calculator.py This is a modified version of the CDF score calculator in Doench et al. (2016) for use in
Guidemaker (https://guidemaker.org)
Adam Rivers, USDA Agricultural Research Service

We score only the CFD for off targets with a NGG site, we do not collect these non-matching PAM off targets Guidemaker.
For this reason we omit the PAM scoring portion of CFD. For that reason we omit the pam scoring
part of the doench et al. (2016) script.  Results are identical for all off-targets that are scored.

Very few off trargets with non-pam matching sites would interact with targets in a small geneome (The highest scoring
non-Pam,NGT, has a score of 0.3). Additionally we require all our guides have a distance of at least 2 by default so
any off targets would have a score below the 0.2 threshold most people use.

We also modified the script to score pam sites longer than 20 by ignoring the 5' end past 20 and for shorter pam's by
only scoring the sites present.
"""
import json
from typing import Tuple, Dict
import logging

logger = logging.getLogger(__name__)

def get_mm_pam_scores() -> Tuple[Dict, Dict]:
    """load json file of mismatch scores and PAM scores

    Returns:
        (tuple):dict of mismatch scores, dict of pam scores

    """
    try:
        with open("guidemaker/data/cfd_data.json") as dat:
            scores = json.load(dat)
        mm_s = scores['mm']
        pam_s = scores['pam']
        return mm_s, pam_s
    except (FileNotFoundError, IOError):
        raise Exception("Could not find file with reference mismatch scores and PAM scores")


def check_len(wt: str, off: str) -> int:
    """Verify the lengths of guide and off target match returning the length

    Args:
        wt: the guide type guide sequence
        off: the off target sequence

    Returns:
        (int): the length of the data

    """
    wtl = len(wt)
    offl = len(off)
    assert (wtl == offl), "The lengths wt and off differ: wt = {}, off = {}".format(str(wtl), str(offl))
    return wtl


def calc_cfd(wt: str, off: str, mm_scores=None) -> float:
    """Calculate the CFD score using precalculated weights

    Args:
        wt: wild-type gRNA sequence, excluding the PAM Cas9 site
        off: off target sequence, excluding the PAM Cas9 site

    Returns:
        (float): CDF score of the pair

    """
    guidelen = check_len(wt, off)
    if mm_scores is None:
        mm_scores, _ = get_mm_pam_scores()
    score = 1.
    off = off.upper().replace('T', 'U')
    wt = wt.upper().replace('T', 'U')
    s_list = list(off)
    wt_list = list(wt)
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    for i, sl in enumerate(s_list):
        if (guidelen - 20 - i) <= 0:
            if wt_list[i] != sl:
                key = 'r' + wt_list[i] + ':d' + basecomp[sl] + ',' + str(20 + i + 1 - guidelen)
                score *= mm_scores[key]
    return score
