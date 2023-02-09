"""
Module providing a dense coancestry matrix implementation using the Yang et al. (2010)
method and associated error checking routines.
"""

import numpy
import numbers
from typing import Union
from typing import Any

from pybrops.core.error import check_all_equal
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_ndarray_axis_len
from pybrops.core.error import check_ndarray_dtype
from pybrops.core.error import check_ndarray_ndim
from pybrops.core.error import check_float_in_interval
from pybrops.core.error import check_ndarray_in_interval
from pybrops.core.error import check_ndarray_dtype_is_object
from pybrops.core.error import check_isinstance
from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import check_is_GenotypeMatrix

class DenseYangCoancestryMatrix(DenseCoancestryMatrix):
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
    def __init__(self, mat: numpy.ndarray, taxa: Union[numpy.ndarray,None] = None, taxa_grp: Union[numpy.ndarray,None] = None, **kwargs: dict):
        """
        Constructor for the concrete class DenseYangCoancestryMatrix.

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
        super(DenseYangCoancestryMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Taxa Data Properites #################
    def taxa():
        doc = "The taxa property."
        def fget(self):
            return self._taxa
        def fset(self, value):
            if value is not None:
                check_is_ndarray(value, "taxa")
                check_ndarray_dtype_is_object(value, "taxa")
                check_ndarray_ndim(value, "taxa", 1)
                check_ndarray_axis_len(value, "taxa", 0, self._mat.shape[0])
            self._taxa = value
        def fdel(self):
            del self._taxa
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    taxa = property(**taxa())

    def taxa_grp():
        doc = "The taxa_grp property."
        def fget(self):
            return self._taxa_grp
        def fset(self, value):
            if value is not None:
                check_is_ndarray(value, "taxa_grp")
                check_ndarray_dtype(value, "taxa_grp", numpy.int64)
                check_ndarray_ndim(value, "taxa_grp", 1)
                check_ndarray_axis_len(value, "taxa_grp", 0, self._mat.shape[0])
            self._taxa_grp = value
        def fdel(self):
            del self._taxa_grp
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    taxa_grp = property(**taxa_grp())

    ############### Taxa Metadata Properites ###############
    def taxa_grp_name():
        doc = "The taxa_grp_name property."
        def fget(self):
            return self._taxa_grp_name
        def fset(self, value):
            if value is not None:
                check_is_ndarray(value, "taxa_grp_name")
                check_ndarray_dtype(value, "taxa_grp_name", numpy.int64)
                check_ndarray_ndim(value, "taxa_grp_name", 1)
            self._taxa_grp_name = value
        def fdel(self):
            del self._taxa_grp_name
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    taxa_grp_name = property(**taxa_grp_name())

    def taxa_grp_stix():
        doc = "The taxa_grp_stix property."
        def fget(self):
            return self._taxa_grp_stix
        def fset(self, value):
            if value is not None:
                check_is_ndarray(value, "taxa_grp_stix")
                check_ndarray_dtype(value, "taxa_grp_stix", numpy.int64)
                check_ndarray_ndim(value, "taxa_grp_stix", 1)
            self._taxa_grp_stix = value
        def fdel(self):
            del self._taxa_grp_stix
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    taxa_grp_stix = property(**taxa_grp_stix())

    def taxa_grp_spix():
        doc = "The taxa_grp_spix property."
        def fget(self):
            return self._taxa_grp_spix
        def fset(self, value):
            if value is not None:
                check_is_ndarray(value, "taxa_grp_spix")
                check_ndarray_dtype(value, "taxa_grp_spix", numpy.int64)
                check_ndarray_ndim(value, "taxa_grp_spix", 1)
            self._taxa_grp_spix = value
        def fdel(self):
            del self._taxa_grp_spix
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    taxa_grp_spix = property(**taxa_grp_spix())

    def taxa_grp_len():
        doc = "The taxa_grp_len property."
        def fget(self):
            return self._taxa_grp_len
        def fset(self, value):
            if value is not None:
                check_is_ndarray(value, "taxa_grp_len")
                check_ndarray_dtype(value, "taxa_grp_len", numpy.int64)
                check_ndarray_ndim(value, "taxa_grp_len", 1)
            self._taxa_grp_len = value
        def fdel(self):
            del self._taxa_grp_len
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    taxa_grp_len = property(**taxa_grp_len())

    ############## Coancestry Data Properites ##############
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            check_is_ndarray(value, "mat")
            check_all_equal(value.shape, "mat.shape")
            check_ndarray_dtype(value, "mat", numpy.float64)
            check_ndarray_ndim(value, "mat", 2)
            self._mat = value
        def fdel(self):
            del self._mat
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    mat = property(**mat())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    @classmethod
    def from_gmat(cls, gmat: GenotypeMatrix, p_anc: Union[numpy.ndarray,float,None] = None, **kwargs: dict):
        """
        Create a dense genomic relationship matrix using methods from Yang et al. (2010) from a GenotypeMatrix.
        
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
        out : DenseYangCoancestryMatrix
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
            check_float_in_interval(p_anc, "p_anc", 0.0, 1.0)
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

        # calculate scaling factor as 1/sqrt(2p(1-p))
        Z_scale = 1.0 / numpy.sqrt(float(gmat.ploidy) * p_anc * (1.0 - p_anc))

        # scale matrix columns by their variances
        Z = Z * Z_scale

        # calculate the scaling coefficient for ZZ'
        # 1 / m
        G_scale = 1.0 / gmat.nvrnt

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



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseYangCoancestryMatrix(v: Any) -> bool:
    """
    Determine whether an object is a DenseYangCoancestryMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseYangCoancestryMatrix object instance.
    """
    return isinstance(v, DenseYangCoancestryMatrix)

def check_is_DenseYangCoancestryMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type DenseYangCoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseYangCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseYangCoancestryMatrix".format(vname))
