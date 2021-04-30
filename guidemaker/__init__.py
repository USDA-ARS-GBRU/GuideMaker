"""Guidemaker: globally design gRNAs for any CRISPR-Cas system in any small genome."""
import pkg_resources
from ._version import get_versions
from .core import *
from .cli import *
from .definitions import *
from .app import *

CONFIG_PATH = pkg_resources.resource_filename('guidemaker', 'data/config_default.yaml')
APP_PARAMETER_IMG = pkg_resources.resource_filename('guidemaker', 'data/parameters.png')
APP_EXPERIMENT_FILE = pkg_resources.resource_filename('guidemaker', 'data/PooledCRISPRExperiments.md')

__all__ = ["core", "cli", "app","definitions"]

__version__ = get_versions()['version']
del get_versions
