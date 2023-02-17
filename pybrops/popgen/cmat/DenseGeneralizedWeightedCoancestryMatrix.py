"""
Module providing an implementation for a dense generalized weighted genomic relationship matrix.
"""

import math
import numbers
from typing import Any, Union
import numpy
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_ndarray_in_interval
from pybrops.core.error import check_ndarray_ndim
from pybrops.core.error import check_ndarray_axis_len
from pybrops.core.error.error_value_python import check_number_in_interval
from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import check_is_GenotypeMatrix

class DenseGeneralizedWeightedCoancestryMatrix(DenseCoancestryMatrix):
    """
    docstring for DenseGeneralizedWeightedCoancestryMatrix.
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
        Constructor for DenseGeneralizedWeightedCoancestryMatrix.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseGeneralizedWeightedCoancestryMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )
    
    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def from_gmat(
            cls, 
            gmat: GenotypeMatrix, 
            mkrwt: Union[numpy.ndarray,numbers.Number,None] = None,
            afreq: Union[numpy.ndarray,numbers.Number,None] = None, 
            **kwargs: dict
        ) -> 'DenseGeneralizedWeightedCoancestryMatrix':
        """
        Construct a generalized weighted genomic relationship matrix.
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Input genotype matrix.
        mkrwt : numpy.ndarray, numbers.Number, None
            Marker weights to apply to the genotype matrix.
            If None, markers weights are assumed to be 1.
        afreq : numpy.ndarray, numbers.Number, None
            Allele frequency of the mean point.
            If None, allele frequencies are assumed to equal the allele frequency of the input GenotypeMatrix.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : DenseGeneralizedWeightedCoancestryMatrix
            A dense generalized weighted coancestry matrix.
        """
        # check input datatypes
        check_is_GenotypeMatrix(gmat, "gmat")

        # check mkrwt inputs
        if mkrwt is None:
            mkrwt = numpy.full((gmat.nvrnt,), 1.0, dtype = "float64")
        elif isinstance(mkrwt, numpy.ndarray):
            check_is_ndarray(mkrwt, "mkrwt")
            check_ndarray_ndim(mkrwt, "mkrwt", 1)
            check_ndarray_axis_len(mkrwt, "mkrwt", 0, gmat.nvrnt)
        elif isinstance(mkrwt, numbers.Number):
            check_number_in_interval(mkrwt, "mkrwt", 0.0, math.inf)
            mkrwt = numpy.full((gmat.nvrnt,), mkrwt, dtype = "float64")
        else:
            raise TypeError("variable 'mkrwt' must be of type 'numpy.ndarray', 'numbers.Number', or None")
            
        # check afreq inputs
        if afreq is None:
            afreq = gmat.afreq() # (n,p) -> (p,)
        elif isinstance(afreq, numpy.ndarray):
            check_ndarray_ndim(afreq, "afreq", 1)
            check_ndarray_axis_len(afreq, "afreq", 0, gmat.nvrnt)
            check_ndarray_in_interval(afreq, "afreq", 0.0, 1.0)
        elif isinstance(afreq, numbers.Number):
            check_number_in_interval(afreq, "afreq", 0.0, 1.0)
            afreq = numpy.full((gmat.nvrnt,), afreq, dtype = "float64")
        else:
            raise TypeError("variable 'afreq' must be of type 'numpy.ndarray', 'numbers.Number', or None")
        
        # multiply afreq by ploidy level to get the number of alleles on which to center
        # (p,) -> (1,p)
        # (1,p) * scalar -> (1,p)
        M = float(gmat.ploidy) * afreq[None,:]

        # get the genotype matrix as {0,1,2,...}
        # (n,p)
        X = gmat.tacount()
        
        # subtract mean from genotype matrix
        # (n,p) - (1,p) -> (n,p)
        Z = X - M

        # calculate the G matrix:
        # G = ZDZ'
        # (n,p) * (1,p) -> (n,p)
        # (n,p) @ (p,n) -> (n,n)
        G = (Z * mkrwt[None,:]).dot(Z.T)

        ### Construct DenseGeneralizedWeightedCoancestryMatrix
        # copy taxa data if available
        taxa = gmat.taxa.copy() if gmat.taxa is not None else None
        taxa_grp = gmat.taxa_grp.copy() if gmat.taxa_grp is not None else None
        taxa_grp_name = gmat.taxa_grp_name.copy() if gmat.taxa_grp_name is not None else None
        taxa_grp_stix = gmat.taxa_grp_stix.copy() if gmat.taxa_grp_stix is not None else None
        taxa_grp_spix = gmat.taxa_grp_spix.copy() if gmat.taxa_grp_spix is not None else None
        taxa_grp_len = gmat.taxa_grp_len.copy() if gmat.taxa_grp_len is not None else None

        # create output
        out = cls(
            mat = G,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
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
def is_DenseGeneralizedWeightedCoancestryMatrix(v: Any) -> bool:
    """
    Determine whether an object is a ``DenseGeneralizedWeightedCoancestryMatrix``.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``v`` is a ``DenseGeneralizedWeightedCoancestryMatrix`` object instance.
    """
    return isinstance(v, DenseGeneralizedWeightedCoancestryMatrix)

def check_is_DenseGeneralizedWeightedCoancestryMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``DenseGeneralizedWeightedCoancestryMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseGeneralizedWeightedCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseGeneralizedWeightedCoancestryMatrix".format(vname))
