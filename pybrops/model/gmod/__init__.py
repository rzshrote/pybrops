"""
Module containing genomic prediction models.
"""

__all__ = [
    "GenomicModel",
    "NonlinearGenomicModel",
    "LinearGenomicModel",
    "CoancestryLinearGenomicModel",
    "AdditiveLinearGenomicModel",
    "AdditiveDominanceLinearGenomicModel",
    "AdditiveDominanceEpistaticLinearGenomicModel",
    "DenseLinearGenomicModel",
    "DenseAdditiveLinearGenomicModel",
]

# imports are order dependent!!!

# utilities

# abstract interfaces
from pybrops.model.gmod import GenomicModel
from pybrops.model.gmod import NonlinearGenomicModel
from pybrops.model.gmod import LinearGenomicModel
from pybrops.model.gmod import CoancestryLinearGenomicModel
from pybrops.model.gmod import AdditiveLinearGenomicModel
from pybrops.model.gmod import AdditiveDominanceLinearGenomicModel
from pybrops.model.gmod import AdditiveDominanceEpistaticLinearGenomicModel

# concrete class implementations
from pybrops.model.gmod import DenseLinearGenomicModel
from pybrops.model.gmod import DenseAdditiveLinearGenomicModel
