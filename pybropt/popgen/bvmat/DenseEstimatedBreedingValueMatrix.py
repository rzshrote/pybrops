import numpy
from . import DenseBreedingValueMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.error import check_ndarray_dtype
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_ndim

class DenseEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """docstring for DenseEstimatedBreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, trait = None, taxa = None, taxa_grp = None, **kwargs):
        """
        Parameters
        ----------
        mat : numpy.ndarray
            A float64 matrix of breeding values of shape (n, t).
        trait : numpy.ndarray
            A string matrix of trait names of shape (t,).
        taxa : numpy.ndarray
            A string_ matrix of taxa names of shape (n,).
        taxa_grp : numpy.ndarray
            An int64 matrix of taxa group labels of shape (n,).
        """
        super(DenseEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            trait = trait,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )
        # self.mat = mat            # already checked in super constructor
        self.trait = trait
        self.taxa = taxa
        self.taxa_grp = taxa_grp

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
            pybropt.util.cond_check_matrix_dtype(value, "taxa", numpy.string_)
            cond_check_ndarray_ndim(value, "taxa", 1)
            pybropt.util.cond_check_matrix_axis_len(value, "taxa", 0, self._mat.shape[0])
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
            pybropt.util.cond_check_matrix_dtype(value, "taxa_grp", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp", 1)
            pybropt.util.cond_check_matrix_axis_len(value, "taxa_grp", 0, self._gmat.geno.shape[0])
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
            pybropt.util.cond_check_matrix_dtype(value, "taxa_grp_name", numpy.int64)
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
            pybropt.util.cond_check_matrix_dtype(value, "taxa_grp_stix", numpy.int64)
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
            pybropt.util.cond_check_matrix_dtype(value, "taxa_grp_spix", numpy.int64)
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
            pybropt.util.cond_check_matrix_dtype(value, "taxa_grp_len", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_len", 1)
            self._taxa_grp_len = value
        def fdel(self):
            del self._taxa_grp_len
        return locals()
    taxa_grp_len = property(**taxa_grp_len())

    ################# Breeding Value Data ##################
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            check_is_ndarray(value, "mat")
            check_ndarray_dtype(value, "mat", numpy.float64)
            check_ndarray_ndim(value, "mat", 2)
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())

    def trait():
        doc = "The trait property."
        def fget(self):
            return self._trait
        def fset(self, value):
            cond_check_is_ndarray(value, "trait")
            pybropt.util.cond_check_matrix_dtype(value, "trait", numpy.string_)
            cond_check_ndarray_ndim(value, "trait", 1)
            pybropt.util.cond_check_matrix_axis_len(value, "trait", 0, self._mat.shape[1])
            self._trait = value
        def fdel(self):
            del self._trait
        return locals()
    trait = property(**trait())

    def ntrait():
        doc = "The ntrait property."
        def fget(self):
            return self.mat.shape[1]
        def fset(self, value):
            pybropt.util.error_readonly("ntrait")
        def fdel(self):
            pybropt.util.error_readonly("ntrait")
        return locals()
    ntrait = property(**ntrait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def axis_index(self, axis):
        """
        Return an index (unsigned) from a provided axis integer (signed)

        Parameters
        ----------
        axis : int
            Integer representation of the axis. Can be in range (-ndim,ndim).
            If outside this range, will raise an AxisError.

        Returns
        -------
        index : int
            Index representation of the axis. In range [0,ndim).
        """
        # get the number of axis
        naxes = self._mat.ndim

        # handle axis argument
        if (axis >= naxes) or (axis <= -naxes):
            raise numpy.AxisError("axis %s is out of bounds for array of dimension %s" % (axis, naxes))

        # modulo the axis number to get the axis (in the case of negative axis)
        axis %= naxes

        return axis

    def lexsort(self, keys = None, axis = -1):
        """
        Perform an indirect stable sort using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using vrnt_chrgrp as primary key, and
            vrnt_phypos as secondary key.
        axis : int
            The axis of the Matrix over which to sort values.

        Returns
        -------
        indices : numpy.ndarray
            Array of indices that sort the keys.
        """
        # transform axis number to an index
        axis = self.axis_index(axis)

        # if no keys were provided, set a default
        if keys is None:
            # assign default keys
            if axis == 0:
                keys = (self._taxa, self._taxa_grp)
            elif axis == 1:
                keys = (self._trait,)

            # remove keys that are None
            keys = tuple(k for k in keys if k is not None)

            # check for errors
            if len(keys) == 0:
                raise RuntimeError("cannot lexsort on axis %s: no default keys" % axis)
        else:
            l = self._mat.shape[axis]
            for i,k in enumerate(keys):
                if len(k) != l:
                    raise RuntimeError("cannot lexsort on axis %s: key %s is incompatible with axis length %s" % (axis, i, l))

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def reorder(self, indices, axis = -1):
        """
        Reorder the VariantMatrix.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        axis : int
            The axis over which to reorder values.

        """
        # transform axis number to an index
        axis = self.axis_index(axis)

        ########################################################################
        if axis == 0:                                           ### TAXA AXIS
            self._mat = self._mat[indices,:]                  # reorder mat array
            if self._taxa is not None:
                self._taxa = self._taxa[indices]                # reorder taxa array
            if self._taxa_grp is not None:
                self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array
        ########################################################################
        elif axis == 1:                                         ### LOCUS AXIS
            self._mat = self._mat[:,indices]                  # reorder mat array
            if self._trait is not None:
                self._trait = self._trait[indices]              # reorder trait array

    def sort(self, keys = None, axis = -1):
        """
        Reset metadata for corresponding axis: name, stix, spix, len.
        Sort the VariantMatrix using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using vrnt_chrgrp as primary key, and
            vrnt_phypos as secondary key.
        axis : int
            The axis over which to sort values.
        """
        # get axis
        axis = self.axis_index(axis)

        if axis == 0:
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
    def group(self, axis = -1):
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        # get axis index
        axis = self.axis_index(axis)

        # sort along taxa axis
        self.sort(axis = axis)

        if axis == 0:
            if self._taxa_grp is not None:
                # get unique taxa group names, starting indices, group lengths
                uniq = numpy.unique(self._taxa_grp, return_index = True, return_counts = True)
                # make assignments to instance data
                self._taxa_grp_name, self._taxa_grp_stix, self._taxa_grp_len = uniq
                # calculate stop indices
                self._taxa_grp_spix = self._taxa_grp_stix + self._taxa_grp_len

    def is_grouped(self, axis = -1):
        """
        Determine whether the Matrix has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        # convert axis to index
        axis = self.axis_index(axis)

        if axis == 1:
            return (
                (self._taxa_grp_name is not None) and
                (self._taxa_grp_stix is not None) and
                (self._taxa_grp_spix is not None) and
                (self._taxa_grp_len is not None)
            )
        return False



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseEstimatedBreedingValueMatrix(v):
    return isinstance(v, DenseEstimatedBreedingValueMatrix)

def check_is_DenseEstimatedBreedingValueMatrix(v, varname):
    if not isinstance(v, DenseEstimatedBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseEstimatedBreedingValueMatrix." % varname)

def cond_check_is_DenseEstimatedBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseEstimatedBreedingValueMatrix(v, varname)
