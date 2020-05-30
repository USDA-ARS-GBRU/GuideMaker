"""Guidemaker: globally design gRNAs for any CRISPR-Cas system in any genome

"""

from .core import *
from .cli import *

__all__ = ["core", "cli"]

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
