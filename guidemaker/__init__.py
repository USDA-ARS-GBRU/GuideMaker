"""GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems."""
import pkg_resources
from ._version import get_versions
from .core import *
from .cli import *
from .definitions import *

CONFIG_PATH = pkg_resources.resource_filename('guidemaker', 'data/config_default.yaml')
APP_PARAMETER_IMG = pkg_resources.resource_filename('guidemaker', 'data/parameters.png')
APP_EXPERIMENT_FILE = pkg_resources.resource_filename('guidemaker', 'data/PooledCRISPRExperiments.md')
WEB_APP = pkg_resources.resource_filename('guidemaker','data/app.py')

__all__ = ["core", "cli", "definitions", "doench_featurization", "doench_predict", "cfd_score_calculator"]

__version__ = get_versions()['version']
del get_versions
