"""GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems."""

import importlib_resources
from ._version import get_versions
from .core import *
from .cli import *
from .definitions import *
from .cfd_score_calculator import *
from .doench_predict import *
from .doench_featurization import *



CONFIG_PATH = importlib_resources.files('guidemaker') / 'data/config_default.yaml'
APP_PARAMETER_IMG = importlib_resources.files('guidemaker') / 'data/parameters.png'
APP_EXPERIMENT_FILE = importlib_resources.files('guidemaker') / 'data/PooledCRISPRExperiments.md'
WEB_APP = importlib_resources.files('guidemaker') / 'data/app.py'

__all__ = ["core", "cli", "definitions", "doench_featurization", "doench_predict", "cfd_score_calculator"]

__version__ = get_versions()['version']
del get_versions
