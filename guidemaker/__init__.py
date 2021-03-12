"""Guidemaker: globally design gRNAs for any CRISPR-Cas system in any small genome

"""
import pkg_resources


CONFIG_PATH = pkg_resources.resource_filename('guidemaker', 'data/config_default.yaml')
APP_PARAMETER_IMG = pkg_resources.resource_filename('guidemaker', 'data/parameters.png')
APP_EXPERIMENT_FILE = pkg_resources.resource_filename('guidemaker', 'data/PooledCRISPRExperiments.md')



from ._version import get_versions
from .core import *
from .cli import *

__all__ = ["core", "cli"]


__version__ = get_versions()['version']
del get_versions
