import numpy

from pybropt.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix

from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.error import check_ndarray_dtype
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_ndim
from pybropt.core.error import cond_check_ndarray_dtype
from pybropt.core.error import cond_check_ndarray_axis_len
from pybropt.core.error import check_all_equal
from pybropt.core.error import cond_check_ndarray_dtype_is_object

class DenseMolecularCoancestryMatrix(DenseCoancestryMatrix):
    """docstring for DenseMolecularCoancestryMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, **kwargs):
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

    ################# Taxa Data Properites #################
    def taxa():
        doc = "The taxa property."
        def fget(self):
            return self._taxa
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa")
            cond_check_ndarray_dtype_is_object(value, "taxa")
            cond_check_ndarray_ndim(value, "taxa", 1)
            cond_check_ndarray_axis_len(value, "taxa", 0, self._mat.shape[0])
            self._taxa = value
        def fdel(self):
            del self._taxa
        return locals()
    taxa = property(**taxa())

    def taxa_grp():
        doc = "The taxa_grp property."
        def fget(self):
            return self._taxa_grp
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp")
            cond_check_ndarray_dtype(value, "taxa_grp", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp", 1)
            cond_check_ndarray_axis_len(value, "taxa_grp", 0, self._mat.shape[0])
            self._taxa_grp = value
        def fdel(self):
            del self._taxa_grp
        return locals()
    taxa_grp = property(**taxa_grp())

    ############### Taxa Metadata Properites ###############
    def taxa_grp_name():
        doc = "The taxa_grp_name property."
        def fget(self):
            return self._taxa_grp_name
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_name")
            cond_check_ndarray_dtype(value, "taxa_grp_name", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_name", 1)
            self._taxa_grp_name = value
        def fdel(self):
            del self._taxa_grp_name
        return locals()
    taxa_grp_name = property(**taxa_grp_name())

    def taxa_grp_stix():
        doc = "The taxa_grp_stix property."
        def fget(self):
            return self._taxa_grp_stix
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_stix")
            cond_check_ndarray_dtype(value, "taxa_grp_stix", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_stix", 1)
            self._taxa_grp_stix = value
        def fdel(self):
            del self._taxa_grp_stix
        return locals()
    taxa_grp_stix = property(**taxa_grp_stix())

    def taxa_grp_spix():
        doc = "The taxa_grp_spix property."
        def fget(self):
            return self._taxa_grp_spix
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_spix")
            cond_check_ndarray_dtype(value, "taxa_grp_spix", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_spix", 1)
            self._taxa_grp_spix = value
        def fdel(self):
            del self._taxa_grp_spix
        return locals()
    taxa_grp_spix = property(**taxa_grp_spix())

    def taxa_grp_len():
        doc = "The taxa_grp_len property."
        def fget(self):
            return self._taxa_grp_len
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_len")
            cond_check_ndarray_dtype(value, "taxa_grp_len", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_len", 1)
            self._taxa_grp_len = value
        def fdel(self):
            del self._taxa_grp_len
        return locals()
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
        return locals()
    mat = property(**mat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys = None, axis = None):
        """
        Performn an indirect stable sort using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using vrnt_chrgrp as primary key, and
            vrnt_phypos as secondary key.
        axis : int
            Ignored by this class

        Returns
        -------
        indices : numpy.ndarray
            Array of indices that sort the keys.
        """
        # if no keys were provided, set a default
        if keys is None:
            keys = (self._taxa, self._taxa_grp)

        # remove keys that are None
        keys = tuple(k for k in keys if k is not None)

        # check for errors
        if len(keys) == 0:
            raise RuntimeError("cannot lexsort: no default keys")

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def reorder(self, indices, axis = None):
        """
        Reorder the CoancestryMatrix along both axes, maintaining symmetry.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        axis : int
            Ignored by this class.
        """
        # reorder along axis 0, then along axis 1
        self._mat[indices,:]
        self._mat[:,indices]

    def sort(self, keys = None, axis = None):
        # reset taxa group metadata
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

        # get indices for sort
        indices = self.lexsort(keys, axis)

        # reorder internals
        self.reorder(indices, axis)

    ################### Grouping Methods ###################
    def group(self, axis = None):
        """

        """
        # sort matrix
        self.sort(axis = axis)

        # create metadata
        if self._taxa_grp is not None:
            # get unique taxa group names, starting indices, group lengths
            uniq = numpy.unique(self._taxa_grp, return_index = True, return_counts = True)
            # make assignments to instance data
            self._taxa_grp_name, self._taxa_grp_stix, self._taxa_grp_len = uniq
            # calculate stop indices
            self._taxa_grp_spix = self._taxa_grp_stix + self._taxa_grp_len

    def is_grouped(self, axis = None):
        """
        Determine whether the CoancestryMatrix has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the CoancestryMatrix has been
            sorted and grouped.
        """
        return (
            (self._taxa_grp_name is not None) and
            (self._taxa_grp_stix is not None) and
            (self._taxa_grp_spix is not None) and
            (self._taxa_grp_len is not None)
        )

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @classmethod
    def from_gmat(cls, gmat):
        """
        Create a CoancestryMatrix from a GenotypeMatrix.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Input genotype matrix from which to calculate coancestry.
        kwargs : dict
            Additional keyword arguments.
        """
        ####################################################
        ### Calculate the coancestry matrix.

        # get ploidy level and reciprocol of number of loci
        ploidy = gmat.ploidy
        rnvrnt = 1.0 / gmat.nvrnt

        # check if we have subroutines to calculate coancestry
        if ploidy not in [1,2]:
            raise RuntimeError("Genotype ploidy level {0} not supported".format(ploidy))

        # get genotype matrix
        X = gmat.tacount()

        # declare pointer to matrix
        mat = None

        # test ploidy level and apply appropriate coancestry calculation
        if ploidy == 1:
            Y = 1 - X                                   # calculate complement to X
            mat = rnvrnt * ((X @ X.T) + (Y @ Y.T))      # (1/m)(XX' + YY')
        elif ploidy == 2:
            X -= 1                                      # {-1,0,1} format
            mat = 0.5 * (1.0 + (rnvrnt * (X @ X.T)))    # (1/2)(1+((1/m)XX'))
        ####################################################

        ####################################################
        ### Construct DenseMolecularCoancestryMatrix
        # copy taxa data if available
        taxa = numpy.object_(gmat.taxa) if gmat.taxa is not None else None
        taxa_grp = numpy.int64(gmat.taxa_grp) if gmat.taxa_grp is not None else None

        # construct basic object
        dmcmat = cls(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp
        )

        # copy taxa metadata if available
        dmcmat.taxa_grp_name = gmat.taxa_grp_name
        dmcmat.taxa_grp_stix = gmat.taxa_grp_stix
        dmcmat.taxa_grp_spix = gmat.taxa_grp_spix
        dmcmat.taxa_grp_len = gmat.taxa_grp_len

        # return matrix
        return dmcmat



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseMolecularCoancestryMatrix(v):
    """
    Determine whether an object is a DenseMolecularCoancestryMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseMolecularCoancestryMatrix object instance.
    """
    return isinstance(v, DenseMolecularCoancestryMatrix)

def check_is_DenseMolecularCoancestryMatrix(v, vname):
    """
    Check if object is of type DenseMolecularCoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseMolecularCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseMolecularCoancestryMatrix".format(vname))

def cond_check_is_DenseMolecularCoancestryMatrix(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DenseMolecularCoancestryMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        DenseMolecularCoancestryMatrix.
    """
    if cond(v):
        check_is_DenseMolecularCoancestryMatrix(v, vname)
