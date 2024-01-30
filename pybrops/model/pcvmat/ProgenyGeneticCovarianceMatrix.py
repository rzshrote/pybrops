"""
Module defining interfaces and error checking routines for matrices storing progeny genetic variance-covariance estimates.
"""

from abc import ABCMeta
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.SquareTaxaSquareTraitMatrix import SquareTaxaSquareTraitMatrix

class ProgenyGeneticCovarianceMatrix(
        SquareTaxaSquareTraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrices storing progeny covariances for specific crosses.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of genetic covariances from a genomic model.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for ProgenyGeneticCovarianceMatrix.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(ProgenyGeneticCovarianceMatrix, self).__init__(**kwargs)

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
