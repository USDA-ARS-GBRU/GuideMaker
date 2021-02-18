"""Guidemaker: globally design gRNAs for any CRISPR-Cas system in any genome

"""


from ._version import get_versions
from .core import *
from .cli import *
from .getsize import *

__all__ = ["core", "cli"]


__version__ = get_versions()['version']
del get_versions
