"""
Module containing genomic prediction models.
"""
# imports are order dependent!!!

# utilities

# abstract classes
# level 1
from . import GenomicModel

# level 2
from . import LinearGenomicModel
from . import NonlinearGenomicModel

# concrete classes
# level 3
from . import AdditiveLinearGenomicModel
