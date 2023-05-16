"""
Module defining phenotype dataframes.
"""

__all__ = [
    "PhenotypeDataFrame",
    "DictPhenotypeDataFrame"
]

# abstract interfaces
from pybrops.popgen.ptdf import PhenotypeDataFrame

# concrete implementations
from pybrops.popgen.ptdf import DictPhenotypeDataFrame
