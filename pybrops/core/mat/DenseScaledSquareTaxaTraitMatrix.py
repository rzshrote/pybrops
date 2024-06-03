"""
Module implementing a dense, scaled matrix with taxa axes that are square and
a trait axis which is not square, and associated error checking routines.
"""

import copy
from numbers import Real
from typing import Optional, Union
import numpy
from pybrops.core.mat.DenseScaledMatrix import DenseScaledMatrix
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.core.mat.ScaledSquareTaxaTraitMatrix import ScaledSquareTaxaTraitMatrix


class DenseScaledSquareTaxaTraitMatrix(
        DenseSquareTaxaTraitMatrix,
        DenseScaledMatrix,
        ScaledSquareTaxaTraitMatrix,
    ):
    """
    A concrete class for dense, scaled matrices with taxa axes that are square 
    and a trait axis which is not square.

    The purpose of this abstract class is to merge the following implementations
    and interfaces:

        1. DenseSquareTaxaTraitMatrix (implementation)
        2. DenseScaledMatrix (implementation)
        3. ScaledSquareTaxaTraitMatrix (interface)
    """

    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self,
            mat: numpy.ndarray, 
            location: Union[numpy.ndarray,Real] = 0.0, 
            scale: Union[numpy.ndarray,Real] = 1.0, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseScaledSquareTaxaTraitMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            Matrix values to store.

        location : numpy.ndarray, Real
            An array of shape ``(t,)`` containing locations for each trait.
            If ``Real``, then the provided location is used for each trait.

        scale : numpy.ndarray, Real
            An array of shape ``(t,)`` containing scales for each trait.
            If ``Real``, then the provided scale is used for each trait.

        taxa : numpy.ndarray, None
            An array of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.

        taxa_grp : numpy.ndarray, None
            An array of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.

        trait : numpy.ndarray, None
            An array of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.

        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call DenseSquareTaxaTraitMatrix constructor
        super(DenseScaledSquareTaxaTraitMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )
        # set location and scale
        self.location = location
        self.scale = scale

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A copy of the DenseScaledSquareTaxaTraitMatrix.
        """
        return self.__class__(
            mat = copy.copy(self.mat),
            location = copy.copy(self.location),
            scale = copy.copy(self.scale),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            trait = copy.copy(self.trait),
        )
    
    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict, None, default = None
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A deep copy of the DenseScaledSquareTaxaTraitMatrix.
        """
        return self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            location = copy.deepcopy(self.location, memo),
            scale = copy.deepcopy(self.scale, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            trait = copy.deepcopy(self.trait, memo),
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a shallow copy of the DenseScaledSquareTaxaTraitMatrix.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A shallow copy of the original DenseScaledSquareTaxaTraitMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a deep copy of the DenseScaledSquareTaxaTraitMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A deep copy of the original DenseScaledSquareTaxaTraitMatrix.
        """
        return copy.deepcopy(self, memo)

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseScaledSquareTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseScaledSquareTaxaTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseScaledSquareTaxaTraitMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseScaledSquareTaxaTraitMatrix.__name__,
                type(v).__name__
            )
        )
