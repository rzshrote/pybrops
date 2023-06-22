"""
Module providing a dense coancestry matrix implementation using the VanRaden
method and associated error checking routines.
"""

__all__ = [
    "DenseVanRadenCoancestryMatrix",
    "check_is_DenseVanRadenCoancestryMatrix"
]

import numpy
from typing import Union

from pybrops.core.error.error_value_python import check_all_equal
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_type_numpy import check_ndarray_dtype
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_in_interval
from pybrops.core.error.error_value_numpy import check_ndarray_in_interval
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_object
from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import check_is_GenotypeMatrix

class DenseVanRadenCoancestryMatrix(DenseCoancestryMatrix):
    """
    A concrete class for a dense coancestry matrix calculated using the VanRaden
    method. Coancestry matrices are square.

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
            taxa: Union[numpy.ndarray,None] = None, 
            taxa_grp: Union[numpy.ndarray,None] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseVanRadenCoancestryMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array from which to construct genomic relationship matrix
        taxa : numpy.ndarray, None, default = None
            Names of taxa.
        taxa_grp : numpy.ndarray, None, default = None
            Taxa group assignments.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseVanRadenCoancestryMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Coancestry Data Properites ##############
    @DenseCoancestryMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        check_is_ndarray(value, "mat")
        check_all_equal(value.shape, "mat.shape")
        check_ndarray_dtype(value, "mat", numpy.float64)
        check_ndarray_ndim(value, "mat", 2)
        self._mat = value

    ################# Taxa Data Properites #################
    @DenseCoancestryMatrix.taxa.setter
    def taxa(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa")
            check_ndarray_dtype_is_object(value, "taxa")
            check_ndarray_ndim(value, "taxa", 1)
            check_ndarray_axis_len(value, "taxa", 0, self._mat.shape[0])
        self._taxa = value

    @DenseCoancestryMatrix.taxa_grp.setter
    def taxa_grp(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa_grp")
            check_ndarray_dtype(value, "taxa_grp", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp", 1)
            check_ndarray_axis_len(value, "taxa_grp", 0, self._mat.shape[0])
        self._taxa_grp = value

    ############### Taxa Metadata Properites ###############
    @DenseCoancestryMatrix.taxa_grp_name.setter
    def taxa_grp_name(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa_grp_name")
            check_ndarray_dtype(value, "taxa_grp_name", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_name", 1)
        self._taxa_grp_name = value

    @DenseCoancestryMatrix.taxa_grp_stix.setter
    def taxa_grp_stix(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa_grp_stix")
            check_ndarray_dtype(value, "taxa_grp_stix", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_stix", 1)
        self._taxa_grp_stix = value

    @DenseCoancestryMatrix.taxa_grp_spix.setter
    def taxa_grp_spix(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa_grp_spix")
            check_ndarray_dtype(value, "taxa_grp_spix", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_spix", 1)
        self._taxa_grp_spix = value

    @DenseCoancestryMatrix.taxa_grp_len.setter
    def taxa_grp_len(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa_grp_len")
            check_ndarray_dtype(value, "taxa_grp_len", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_len", 1)
        self._taxa_grp_len = value

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def from_gmat(
            cls, 
            gmat: GenotypeMatrix, 
            p_anc: Union[numpy.ndarray,float,None] = None, 
            **kwargs: dict
        ) -> 'DenseVanRadenCoancestryMatrix':
        """
        Create a DenseVanRadenCoancestryMatrix from a GenotypeMatrix
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Input genotype matrix from which to calculate the genomic relationship matrix.
        p_anc : numpy.ndarray, float, None, default = None
            Ancestral allele frequencies.
            If numpy.ndarray, the array must be of shape (p,) where ``p`` is the number of marker loci.
            If float, ancestral allele frequencies are assumed constant ``p_anc`` across all marker loci. 
            If None, ancestral allele frequencies are estimated from ``gmat``.

        Returns
        -------
        out : DenseVanRadenCoancestryMatrix
            Dense VanRaden genomic relationship matrix.
        """
        # check input datatypes
        check_is_GenotypeMatrix(gmat, "gmat")

        # if None, calculate allele frequencies from provided matrix.
        if p_anc is None:
            p_anc = gmat.afreq()                    # (n,p) -> (p,)

        # if p_anc is numpy.ndarray, check for correct shape, and value ranges 
        elif isinstance(p_anc, numpy.ndarray):
            # check values
            check_ndarray_ndim(p_anc, "p_anc", 1)
            check_ndarray_axis_len(p_anc, "p_anc", 0, gmat.nvrnt)
            check_ndarray_in_interval(p_anc, "p_anc", 0.0, 1.0)
        
        # if p_anc is float, check for correct value range
        elif isinstance(p_anc, float):
            # check values
            check_is_in_interval(p_anc, "p_anc", 0.0, 1.0)
            p_anc = numpy.repeat(p_anc, gmat.nvrnt) # scalar -> (p,)
        
        # raise error for incorrect types
        else:
            raise TypeError("p_anc must be of type 'numpy.ndarray', 'numbers.Number', or None")
        
        # multiply p_anc by ploidy level to get mean number of alleles we expect
        # (p,) -> (1,p)
        # (1,p) * scalar -> (1,p)
        M = p_anc[None,:] * float(gmat.ploidy)

        # get the genotype matrix as {0,1,2,...}
        X = gmat.tacount()
        
        # subtract mean from genotype matrix
        Z = X - M

        # calculate the scaling coefficient for ZZ'
        # 1 / (2p(1-p))
        G_scale = 1.0 / (float(gmat.ploidy) * p_anc.dot(1.0 - p_anc))

        # calculate the G matrix:
        # G = G_scale * ZZ'
        G = G_scale * Z.dot(Z.T)

        # create output
        out = cls(
            mat = G,
            taxa = gmat.taxa,
            taxa_grp = gmat.taxa_grp,
            **kwargs
        )

        # copy taxa metadata if available
        out.taxa_grp_name = gmat.taxa_grp_name
        out.taxa_grp_stix = gmat.taxa_grp_stix
        out.taxa_grp_spix = gmat.taxa_grp_spix
        out.taxa_grp_len = gmat.taxa_grp_len

        # return matrix
        return out



################################## Utilities ###################################
def check_is_DenseVanRadenCoancestryMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseVanRadenCoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseVanRadenCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseVanRadenCoancestryMatrix".format(vname))
