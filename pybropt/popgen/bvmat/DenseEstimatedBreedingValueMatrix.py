import copy
import numpy

from . import DenseBreedingValueMatrix

from pybropt.core.mat import get_axis

from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.error import check_ndarray_dtype
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_ndim
from pybropt.core.error import error_readonly
from pybropt.core.error import cond_check_ndarray_dtype
from pybropt.core.error import cond_check_ndarray_axis_len
from pybropt.core.error import cond_check_ndarray_dtype_is_object

class DenseEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """docstring for DenseEstimatedBreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, raw = None, trait = None, taxa = None, taxa_grp = None, **kwargs):
        """
        Parameters
        ----------
        mat : numpy.ndarray
            A float64 matrix of breeding values of shape (n, t).
        raw : numpy.ndarray
            A float64 matrix of raw phenotypic values of shape (r, n, t).
        trait : numpy.ndarray
            A string matrix of trait names of shape (t,).
        taxa : numpy.ndarray
            A string_ matrix of taxa names of shape (n,).
        taxa_grp : numpy.ndarray
            An int64 matrix of taxa group labels of shape (n,).
        """
        super(DenseEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

        # set instance data
        self.raw = raw
        self.trait = trait
        self.taxa = taxa
        self.taxa_grp = taxa_grp

        # set metadata
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

    def __copy__(self):
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        # construct new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            raw = copy.copy(self.raw),
            trait = copy.copy(self.trait),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
        )

        # copy metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        return out

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Matrix
        """
        # construct new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat),
            raw = copy.deepcopy(self.raw),
            trait = copy.deepcopy(self.trait),
            taxa = copy.deepcopy(self.taxa),
            taxa_grp = copy.deepcopy(self.taxa_grp),
        )

        # copy metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len)

        return out

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

    def ntaxa():
        doc = "The ntaxa property."
        def fget(self):
            return self._mat.shape[0]
        def fset(self, value):
            error_readonly("ntaxa")
        def fdel(self):
            error_readonly("ntaxa")
        return locals()
    ntaxa = property(**ntaxa())

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

    ################## Raw Phenotype Data ##################
    def raw():
        doc = "Raw phenotype matrix of shape (r, n, t). r = rep, n = indiv, t = trait"
        def fget(self):
            return self._raw
        def fset(self, value):
            cond_check_is_ndarray(value, "raw")
            cond_check_ndarray_dtype(value, "raw", numpy.float64)
            cond_check_ndarray_ndim(value, "raw", 3)
            cond_check_ndarray_axis_len(value, "raw", 1, self._mat.shape[0])
            cond_check_ndarray_axis_len(value, "raw", 2, self._mat.shape[1])
            self._raw = value
        def fdel(self):
            del self._raw
        return locals()
    raw = property(**raw())

    ###################### Trait Data ######################
    def trait():
        doc = "The trait property."
        def fget(self):
            return self._trait
        def fset(self, value):
            cond_check_is_ndarray(value, "trait")
            cond_check_ndarray_dtype_is_object(value, "trait")
            cond_check_ndarray_ndim(value, "trait", 1)
            cond_check_ndarray_axis_len(value, "trait", 0, self._mat.shape[1])
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
            error_readonly("ntrait")
        def fdel(self):
            error_readonly("ntrait")
        return locals()
    ntrait = property(**ntrait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############# Matrix element manipulation ##############
    def append(self, values, axis, raw = None, trait = None, taxa = None, taxa_grp = None, **kwargs):
        """
        Append values to the matrix. Cannot add additional locations to 'raw'.

        Parameters
        ----------
        values : numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        if axis == 0:
            if (self._raw is not None) and (raw is None):
                raise RuntimeError("cannot append: raw argument is required")
            if (self._taxa is not None) and (taxa is None):
                raise RuntimeError("cannot append: taxa argument is required")
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise RuntimeError("cannot append: taxa_grp argument is required")
        elif axis == 1:
            if (self._raw is not None) and (raw is None):
                raise RuntimeError("cannot append: raw argument is required")
            if (self._trait is not None) and (trait is None):
                raise RuntimeError("cannto append: trait argument is required")

        # append values
        self._mat = numpy.append(self._mat, values, axis = axis)
        if axis == 0:   # taxa axis
            if (self._raw is not None) and (raw is not None):
                self._raw = numpy.append(self._raw, raw, axis = 1)
            if (self._taxa is not None) and (taxa is not None):
                self._taxa = numpy.append(self._taxa, taxa, axis = 0)
            if (self._taxa_grp is not None) and (taxa_grp is not None):
                self._taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)
        elif axis == 1:
            if (self._raw is not None) and (raw is not None):
                self._raw = numpy.append(self._raw, raw, axis = 2)
            if (self._trait is not None) and (trait is not None):
                self._trait = numpy.append(self._trait, trait, axis = 0)

    def delete(self, obj, axis, **kwargs):
        """
        Delete sub-arrays along an axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to delete the subarray defined by obj.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # delete values
        self._mat = numpy.delete(self._mat, obj, axis = axis)
        if axis == 0:   # taxa axis
            if self._raw is not None:
                self._raw = numpy.delete(self._raw, obj, axis = 1)
            if self._taxa is not None:
                self._taxa = numpy.delete(self._taxa, obj, axis = 0)
            if self._taxa_grp is not None:
                self._taxa_grp = numpy.delete(self._taxa_grp, obj, axis = 0)
        if axis == 2:   # trait axis
            if self._raw is not None:
                self._raw = numpy.delete(self._raw, obj, axis = 2)
            if self._trait is not None:
                self._trait = numpy.delete(self._trait, obj, axis = 0)

    def insert(self, obj, values, axis, raw = None, trait = None, taxa = None, taxa_grp = None, **kwargs):
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        if axis == 0:
            if (self._raw is not None) and (raw is None):
                raise RuntimeError("cannot append: raw argument is required")
            if (self._taxa is not None) and (taxa is None):
                raise RuntimeError("cannot append: taxa argument is required")
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise RuntimeError("cannot append: taxa_grp argument is required")
        elif axis == 1:
            if (self._raw is not None) and (raw is None):
                raise RuntimeError("cannot append: raw argument is required")
            if (self._trait is not None) and (trait is None):
                raise RuntimeError("cannto append: trait argument is required")

        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = axis)
        if axis == 0:   # taxa axis
            if (self._raw is not None) and (raw is not None):
                self._raw = numpy.insert(self._raw, obj, raw, axis = 1)
            if (self._taxa is not None) and (taxa is not None):
                self._taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
            if (self._taxa_grp is not None) and (taxa_grp is not None):
                self._taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)
        elif axis == 1: # trait axis
            if (self._raw is not None) and (raw is not None):
                self._raw = numpy.insert(self._raw, obj, raw, axis = 2)
            if (self._trait is not None) and (trait is not None):
                self._trait = numpy.insert(self._trait, obj, trait, axis = 0)

    def select(self, obj, axis, **kwargs):
        """
        Select certain values from the GenotypeMatrix.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices where values are selected.
        axis : int
            The axis along which values are selected.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # initialize to null pointers
        mat_sel = None
        raw_sel = None
        taxa_sel = None
        taxa_grp_sel = None
        trait_sel = None

        # custom selections along each axis
        if axis == 0:   # taxa axis
            mat_sel = self._mat[obj,:]
            if self._raw is not None:
                raw_sel = self._raw[:,obj,:]
            if self._taxa is not None:
                taxa_sel = self._taxa[obj]
            if self._taxa_grp is not None:
                taxa_grp_sel = self._taxa_grp[obj]
            trait_sel = self._trait
        elif axis == 1: # trait axis
            mat_sel = self._mat[:,obj]
            if self._raw is not None:
                raw_sel = self._raw[:,:,obj]
            taxa_sel = self._taxa
            taxa_grp_sel = self._taxa_grp
            if self._trait is not None:
                trait_sel = self._trait[obj]

        # select elements
        dbvmat = DenseEstimatedBreedingValueMatrix(
            mat = mat_sel,
            raw = raw_sel,
            trait = trait_sel,
            taxa = taxa_sel,
            taxa_grp = taxa_grp_sel
        )

        return dbvmat

    ################### Sorting Methods ####################
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
        axis = get_axis(axis, self._mat.ndim)

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
        axis = get_axis(axis, self._mat.ndim)

        ########################################################################
        if axis == 0:                                           ### TAXA AXIS
            self._mat = self._mat[indices,:]                  # reorder mat array
            if self._raw is not None:
                self._raw = self._raw[:,indices,:]
            if self._taxa is not None:
                self._taxa = self._taxa[indices]                # reorder taxa array
            if self._taxa_grp is not None:
                self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array
        ########################################################################
        elif axis == 1:                                         ### LOCUS AXIS
            self._mat = self._mat[:,indices]                  # reorder mat array
            if self._raw is not None:
                self._raw = self._raw[:,:,indices]
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
        axis = get_axis(axis, self._mat.ndim)

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
        axis = get_axis(axis, self._mat.ndim)

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
        axis = get_axis(axis, self._mat.ndim)

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
