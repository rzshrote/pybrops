"""
Module providing an implementation for a dense generalized weighted genomic relationship matrix.
"""

import math
import numbers
from typing import Union
import numpy
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_ndarray_in_interval
from pybrops.core.error import check_ndarray_ndim
from pybrops.core.error import check_ndarray_axis_len
from pybrops.core.error.error_value_python import check_float_in_interval, check_number_in_interval
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
        ):
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
        ):
        """
        Construct a generalized weighted genomic relationship matrix.
        
        Parameters
        ----------
        gmat : argtype
            argdesc
        
        Returns
        -------
        out : outtype
            outdesc
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
        # G = ZZ'
        G = Z.dot(Z.T)

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
