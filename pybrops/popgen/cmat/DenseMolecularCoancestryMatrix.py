"""
Module providing a dense coancestry matrix implementation for identity by state
and associated error checking routines.
"""

from typing import Any, Optional
import numpy

from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class DenseMolecularCoancestryMatrix(DenseCoancestryMatrix):
    """
    A concrete class for a dense coancestry matrix calculated using molecular
    coancestry (identity by state). Coancestry matrices are square.

    The purpose of this concrete class is to implement functionality for:
        1) Dense coancestry matrix value calculation.
        2) Dense coancestry matrix value access.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseMolecularCoancestryMatrix class.

        Parameters
        ----------
        mat : numpy.ndarray
        taxa : numpy.ndarray
        taxa_grp : numpy.ndarray
        """
        super(DenseMolecularCoancestryMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################


    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    # methods already inherited from DenseSquareTaxaMatrix and DenseCoancestryMatrix
    # ################### Sorting Methods ####################
    # def lexsort(self, keys = None, axis = None):
    #     """
    #     Performn an indirect stable sort using a tuple of keys.

    #     Parameters
    #     ----------
    #     keys : tuple, None
    #         A tuple of columns to be sorted. The last column is the primary
    #         sort key. If None, sort using vrnt_chrgrp as primary key, and
    #         vrnt_phypos as secondary key.
    #     axis : int
    #         Ignored by this class

    #     Returns
    #     -------
    #     indices : numpy.ndarray
    #         Array of indices that sort the keys.
    #     """
    #     # if no keys were provided, set a default
    #     if keys is None:
    #         keys = (self._taxa, self._taxa_grp)

    #     # remove keys that are None
    #     keys = tuple(k for k in keys if k is not None)

    #     # check for errors
    #     if len(keys) == 0:
    #         raise RuntimeError("cannot lexsort: no default keys")

    #     # get indices
    #     indices = numpy.lexsort(keys)

    #     # return indices
    #     return indices

    # def reorder(self, indices, axis = None):
    #     """
    #     Reorder the CoancestryMatrix along both axes, maintaining symmetry.

    #     Parameters
    #     ----------
    #     indices : numpy.ndarray
    #         Indices of where to place elements.
    #     axis : int
    #         Ignored by this class.
    #     """
    #     # reorder along axis 0, then along axis 1
    #     self._mat[indices,:]
    #     self._mat[:,indices]

    # def sort(self, keys = None, axis = None):
    #     # reset taxa group metadata
    #     self.taxa_grp_name = None
    #     self.taxa_grp_stix = None
    #     self.taxa_grp_spix = None
    #     self.taxa_grp_len = None

    #     # get indices for sort
    #     indices = self.lexsort(keys, axis)

    #     # reorder internals
    #     self.reorder(indices, axis)

    # ################### Grouping Methods ###################
    # def group(self, axis = None):
    #     """

    #     """
    #     # sort matrix
    #     self.sort(axis = axis)

    #     # create metadata
    #     if self._taxa_grp is not None:
    #         # get unique taxa group names, starting indices, group lengths
    #         uniq = numpy.unique(self._taxa_grp, return_index = True, return_counts = True)
    #         # make assignments to instance data
    #         self._taxa_grp_name, self._taxa_grp_stix, self._taxa_grp_len = uniq
    #         # calculate stop indices
    #         self._taxa_grp_spix = self._taxa_grp_stix + self._taxa_grp_len

    # def is_grouped(self, axis = None):
    #     """
    #     Determine whether the CoancestryMatrix has been sorted and grouped.

    #     Returns
    #     -------
    #     grouped : bool
    #         True or False indicating whether the CoancestryMatrix has been
    #         sorted and grouped.
    #     """
    #     return (
    #         (self._taxa_grp_name is not None) and
    #         (self._taxa_grp_stix is not None) and
    #         (self._taxa_grp_spix is not None) and
    #         (self._taxa_grp_len is not None)
    #     )

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @classmethod
    def from_gmat(
            cls, 
            gmat: GenotypeMatrix, 
            **kwargs: dict
        ) -> 'DenseMolecularCoancestryMatrix':
        """
        Create a CoancestryMatrix from a GenotypeMatrix.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Input genotype matrix from which to calculate coancestry.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseMolecularCoancestryMatrix
            A dense molecular coancestry matrix.
        """
        ####################################################
        ### Calculate the coancestry matrix.

        # get ploidy level and reciprocol of number of loci
        ploidy = gmat.ploidy
        rnvrnt = 1.0 / gmat.nvrnt

        # check if we have subroutines to calculate coancestry
        if ploidy not in [1,2]:
            raise RuntimeError("Genotype ploidy level {0} not supported".format(ploidy))

        # get genotype matrix as a native integer type
        X = gmat.tacount(int)

        # declare pointer to matrix
        mat = None

        # test ploidy level and apply appropriate coancestry calculation
        if ploidy == 1:
            Y = 1 - X                                   # calculate complement to X
            mat = (2.0 * rnvrnt) * ((X @ X.T) + (Y @ Y.T))  # (2/m)(XX' + YY')
        elif ploidy == 2:
            X -= 1                                      # {-1,0,1} format
            mat = (1.0 + (rnvrnt * (X @ X.T)))          # (1+((1/m)XX'))
        ####################################################

        ####################################################
        ### Construct DenseMolecularCoancestryMatrix
        # copy taxa data if available
        taxa = gmat.taxa.copy() if gmat.taxa is not None else None
        taxa_grp = gmat.taxa_grp.copy() if gmat.taxa_grp is not None else None
        taxa_grp_name = gmat.taxa_grp_name.copy() if gmat.taxa_grp_name is not None else None
        taxa_grp_stix = gmat.taxa_grp_stix.copy() if gmat.taxa_grp_stix is not None else None
        taxa_grp_spix = gmat.taxa_grp_spix.copy() if gmat.taxa_grp_spix is not None else None
        taxa_grp_len = gmat.taxa_grp_len.copy() if gmat.taxa_grp_len is not None else None

        # construct basic object
        out = cls(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp
        )

        # apply taxa metadata
        out.taxa_grp_name = taxa_grp_name
        out.taxa_grp_stix = taxa_grp_stix
        out.taxa_grp_spix = taxa_grp_spix
        out.taxa_grp_len = taxa_grp_len

        # return matrix
        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseMolecularCoancestryMatrix(v: object) -> bool:
    """
    Determine whether an object is a DenseMolecularCoancestryMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseMolecularCoancestryMatrix object instance.
    """
    return isinstance(v, DenseMolecularCoancestryMatrix)

def check_is_DenseMolecularCoancestryMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseMolecularCoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseMolecularCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseMolecularCoancestryMatrix".format(vname))
