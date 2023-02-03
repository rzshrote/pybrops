import copy
from typing import Sequence, Union
import numpy
import h5py

from pybrops.core.error.error_attr_python import check_is_iterable

from . import DenseBreedingValueMatrix

from pybrops.core.mat import get_axis

from pybrops.core.error import check_file_exists
from pybrops.core.error import check_group_in_hdf5
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_is_str
from pybrops.core.error import check_ndarray_ndim
from pybrops.core.error import check_ndarray_dtype
from pybrops.core.error import cond_check_is_ndarray
from pybrops.core.error import cond_check_ndarray_axis_len
from pybrops.core.error import cond_check_ndarray_dtype
from pybrops.core.error import cond_check_ndarray_dtype_is_object
from pybrops.core.error import cond_check_ndarray_ndim
from pybrops.core.error import error_readonly
from pybrops.core.util import save_dict_to_hdf5

class DenseEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """docstring for DenseEstimatedBreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat: numpy.ndarray, raw = None, taxa = None, taxa_grp = None, trait = None, **kwargs: dict):
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
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.trait = trait

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

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values: Union[Matrix,numpy.ndarray], axis = -1, raw = None, taxa = None, taxa_grp = None, trait = None, **kwargs: dict):
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : Matrix or numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are appended.
        raw : numpy.ndarray
            A float64 matrix of raw phenotypic values of shape (r, n, t).
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseEstimatedBreedingValueMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseEstimatedBreedingValueMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait breeding values to adjoin to the Matrix.
            If values is a DenseEstimatedBreedingValueMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a Matrix extract Matrix.mat values
        if is_DenseEstimatedBreedingValueMatrix(values):
            if raw is None:
                raw = values.raw
            if trait is None:
                trait = values.trait
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseEstimatedBreedingValueMatrix or numpy.ndarray")

        # perform error checks
        if axis == 0:
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot adjoin: raw argument is required")
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot adjoin: taxa_grp argument is required")
        elif axis == 1:
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot adjoin: raw argument is required")
            if (self._trait is not None) and (trait is None):
                raise TypeError("cannot adjoin: trait argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # adjoin values
        values = numpy.append(self._mat, values, axis = axis)
        if axis == 0:
            if self._raw is not None:
                raw = numpy.append(self._raw, raw, axis = 1)
            if self._taxa is not None:
                taxa = numpy.append(self._taxa, taxa, axis = 0)
            if self._taxa_grp is not None:
                taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)
        elif axis == 1:
            if self._raw is not None:
                raw = numpy.append(self._raw, raw, axis = 1)
            if self._trait is not None:
                trait = numpy.append(self._trait, trait, axis = 0)

        # create new output
        out = self.__class__(
            mat = values,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait
        )

        return out

    def adjoin_taxa(self, values: Union[Matrix,numpy.ndarray], raw = None, taxa = None, taxa_grp = None, **kwargs: dict):
        """
        Add additional elements to the end of the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        raw : numpy.ndarray
            A float64 matrix of raw phenotypic values of shape (r, n, t).
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        return self.adjoin(
            values = values,
            axis = 0,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def adjoin_trait(self, values: Union[Matrix,numpy.ndarray], raw = None, trait = None, **kwargs: dict):
        """
        Add additional elements to the end of the Matrix along the trait axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        raw : numpy.ndarray
            A float64 matrix of raw phenotypic values of shape (r, n, t).
        trait : numpy.ndarray
            Taxa names to adjoin to the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        return self.adjoin(
            values = values,
            axis = 1,
            raw = raw,
            trait = trait,
            **kwargs
        )

    def delete(self, obj: Union[int,slice,Sequence], axis = -1, **kwargs: dict):
        """
        Delete sub-arrays along an axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to delete the subarray defined by obj.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # get values
        mat = self._mat
        raw = self._raw
        taxa = self._taxa
        taxa_grp = self._taxa_grp
        trait = self._trait

        # delete values
        mat = numpy.delete(mat, obj, axis = axis)
        if axis == 0:
            if raw is not None:
                raw = numpy.delete(raw, obj, axis = 1)
            if taxa is not None:
                taxa = numpy.delete(taxa, obj, axis = 0)
            if taxa_grp is not None:
                taxa_grp = numpy.delete(taxa_grp, obj, axis = 0)
        elif axis == 1:
            if raw is not None:
                raw = numpy.delete(raw, obj, axis = 2)
            if trait is not None:
                trait = numpy.delete(trait, obj, axis = 0)

        # create new output
        out = self.__class__(
            mat = mat,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

        return out

    def delete_taxa(self, obj: Union[int,slice,Sequence], **kwargs: dict):
        """
        Delete sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.delete(
            obj = obj,
            axis = 0,
            **kwargs
        )

    def delete_trait(self, obj: Union[int,slice,Sequence], **kwargs: dict):
        """
        Delete sub-arrays along the trait axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.delete(
            obj = obj,
            axis = 1,
            **kwargs
        )

    def insert(self, obj: Union[int,slice,Sequence], values: Union[Matrix,numpy.ndarray], axis = -1, raw = None, taxa = None, taxa_grp = None, trait = None, **kwargs: dict):
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        raw : numpy.ndarray
            A float64 matrix of raw phenotypic values of shape (r, n, t).
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait breeding values to insert to the Matrix.
            If values is a DenseEstimatedBreedingValueMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DenseEstimatedBreedingValueMatrix(values):
            if raw is None:
                raw = values.raw
            if trait is None:
                trait = values.trait
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseEstimatedBreedingValueMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if axis == 0:
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot insert: raw argument is required")
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot insert: taxa_grp argument is required")
        elif axis == 1:
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot insert: raw argument is required")
            if (self._trait is not None) and (trait is None):
                raise TypeError("cannot insert: trait argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # insert values
        values = numpy.insert(self._mat, obj, values, axis = axis)
        if axis == 0:
            if self._raw is not None:
                raw = numpy.insert(self._raw, obj, raw, axis = 1)
            if self._taxa is not None:
                taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
            if self._taxa_grp is not None:
                taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)
        elif axis == 1:
            if self._raw is not None:
                raw = numpy.insert(self._raw, obj, raw, axis = 2)
            if self._trait is not None:
                trait = numpy.insert(self._trait, obj, trait, axis = 0)

        # create output
        out = self.__class__(
            mat = values,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

        return out

        # if given a Matrix extract Matrix.mat values
        if is_Matrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type Matrix or numpy.ndarray")

        # append values
        mat = numpy.insert(self._mat, obj, values, axis)

        # create new output
        out = self.__class__(mat = mat)

        return out

    def insert_taxa(self, obj: Union[int,slice,Sequence], values: Union[Matrix,numpy.ndarray], raw = None, taxa = None, taxa_grp = None, **kwargs: dict):
        """
        Insert values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.insert(
            obj = obj,
            values = values,
            axis = 0,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def insert_trait(self, obj: Union[int,slice,Sequence], values, raw = None, trait = None, **kwargs: dict):
        """
        Insert values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        trait : numpy.ndarray
            Taxa names to insert into the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.insert(
            obj = obj,
            values = values,
            raw = raw,
            trait = trait,
            **kwargs
        )

    def select(self, indices, axis = -1, **kwargs: dict):
        """
        Select certain values from the matrix.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        axis : int
            The axis along which values are selected.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # get values
        mat = self._mat
        raw = self._raw
        taxa = self._taxa
        taxa_grp = self._taxa_grp
        trait = self._trait

        # select values
        mat = numpy.take(self._mat, indices, axis)
        if axis == 0:
            if raw is not None:
                raw = numpy.take(raw, indices, axis = 1)
            if taxa is not None:
                taxa = numpy.take(taxa, indices, axis = 0)
            if taxa_grp is not None:
                taxa_grp = numpy.take(taxa_grp, indices, axis = 0)
        elif axis == 1:
            if raw is not None:
                raw = numpy.take(raw, indices, axis = 2)
            if trait is not None:
                trait = numpy.take(trait, indices, axis = 0)

        # create new output
        out = self.__class__(
            mat = mat,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

        return out

    def select_taxa(self, indices, **kwargs: dict):
        """
        Select certain values from the Matrix along the taxa axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        return self.select(
            indices = indices,
            axis = 0,
            **kwargs
        )

    def select_trait(self, indices, **kwargs: dict):
        """
        Select certain values from the Matrix along the trait axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        return self.select(
            indices = indices,
            axis = 1,
            **kwargs
        )

    @staticmethod
    def concat(mats, axis = -1, **kwargs: dict):
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        # ensure that we have an iterable object
        check_is_iterable(mats, "mats")

        # get length of mats
        mats_len = len(mats)

        # ensure that we have an array_like of length >= 1
        if mats_len <= 0:
            raise ValueError("need at least one Matrix to concatenate")

        # ensure that all items in mats are DenseEstimatedBreedingValueMatrix
        for i in range(mats_len):
            check_is_DenseEstimatedBreedingValueMatrix(mats[i], "mats[{0}]".format(i))

        # get first matrix
        mats0 = mats[0]

        # get axis
        axis = get_axis(axis, mats0.mat.ndim)

        # extract tuples of shape parameters for Matrix
        ntaxa_t, ntrait_t = zip(*[m.mat.shape for m in mats])

        # extract first Matrix shape parameters
        ntaxa, ntrait = mats0.mat.shape

        # create matrix lists
        mat_l = [m.mat for m in mats]
        raw_l = None
        taxa_l = None
        taxa_grp_l = None
        trait_l = None

        # check shapes and add to list
        if axis == 0:                                           # concatenate additional taxa
            # check matrix shapes
            if any(e != ntrait for e in ntrait_t):              # raise error if any have different trait number
                raise ValueError("Matrix shapes do not all align along axis 2 (trait axis)")
            # add taxa related attributes to lists
            if mats0.raw is not None:
                raw_l = [m.raw for m in mats]
                if any(e is None for e in raw_l):
                    raise ValueError("cannot concat: raw needed for all Matrix in list")
            if mats0.taxa is not None:                          # populate taxa_l
                taxa_l = [numpy.object_([None]*m.ntaxa) if m.taxa is None else m.taxa for m in mats]
            if mats0.taxa_grp is not None:
                taxa_grp_l = [m.taxa_grp for m in mats]
                if any(e is None for e in taxa_grp_l):
                    raise ValueError("cannot concat: taxa_grp needed for all Matrix in list")
        elif axis == 1:                                         # concatenate additional traits
            # check matrix shapes
            if any(e != ntaxa for e in ntaxa_t):                # raise error if any have different taxa number
                raise ValueError("Matrix shapes do not all align along axis 1 (taxa axis)")
            # add loci related attributes to lists
            if mats0.raw is not None:
                raw_l = [m.raw for m in mats]
                if any(e is None for e in raw_l):
                    raise ValueError("cannot concat: raw needed for all Matrix in list")
            if mats0.trait is not None:
                trait_l = [m.trait for m in mats]
                if any(e is None for e in trait_l):
                    raise ValueError("cannot concat: trait needed for all Matrix in list")

        # concatenate everything
        mat = numpy.concatenate(mat_l, axis = axis)
        raw = mats0.raw if raw_l is None else numpy.concatenate(raw_l, axis = axis+1)
        taxa = mats0.taxa if taxa_l is None else numpy.concatenate(taxa_l, axis = 0)
        taxa_grp = mats0.taxa_grp if taxa_grp_l is None else numpy.concatenate(taxa_grp_l, axis = 0)
        trait = mats0.trait if trait_l is None else numpy.concatenate(trait_l, axis = 0)

        # concatenate everything and put into new DensePhasedGenotypeVariantMatrix
        out = DenseEstimatedBreedingValueMatrix(
            mat = mat,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

        return out

    @staticmethod
    def concat_taxa(mats: Sequence, **kwargs: dict):
        """
        Concatenate list of Matrix together along the taxa axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        return self.concat(
            mats = mats,
            axis = 0,
            **kwargs
        )

    @staticmethod
    def concat_trait(mats: Sequence, **kwargs: dict):
        """
        Concatenate list of Matrix together along the trait axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        return self.concat(
            mats = mats,
            axis = 1,
            **kwargs
        )

    ######### Matrix element in-place-manipulation #########
    def append(self, values: Union[Matrix,numpy.ndarray], axis = -1, raw = None, taxa = None, taxa_grp = None, trait = None, **kwargs: dict):
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

        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DenseEstimatedBreedingValueMatrix(values):
            if raw is None:
                raw = values.raw
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            if trait is None:
                trait = values.trait
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseEstimatedBreedingValueMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if axis == 0:
            if self._mat.shape[1] != values.shape[1]:
                raise ValueError("Matrix shapes do not all align along axis 1 (trait axis)")
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot append: raw argument is required")
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot append: taxa_grp argument is required")
        elif axis == 1:
            if self._mat.shape[0] != values.shape[0]:
                raise ValueError("Matrix shapes do not all align along axis 0 (taxa axis)")
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot append: raw argument is required")
            if (self._trait is not None) and (trait is None):
                raise TypeError("cannot append: trait argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # append values
        self._mat = numpy.append(self._mat, values, axis = axis)
        if axis == 0:
            # set fields
            if self._raw is not None:
                self._raw = numpy.append(self._raw, raw, axis = 1)
            if self._taxa is not None:
                self._taxa = numpy.append(self._taxa, taxa, axis = 0)
            if self._taxa_grp is not None:
                self._taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)
            # reset metadata
            self._taxa_grp_len = None
            self._taxa_grp_name = None
            self._taxa_grp_stix = None
            self._taxa_grp_spix = None
        elif axis == 0:
            # set fields
            if self._raw is not None:
                self._raw = numpy.append(self._raw, raw, axis = 2)
            if self._trait is not None:
                self._trait = numpy.append(self._trait, trait, axis = 0)

    def append_taxa(self, values: Union[Matrix,numpy.ndarray], raw = None, taxa = None, taxa_grp = None, **kwargs: dict):
        """
        Append values to the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.append(
            values = values,
            axis = 0,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def append_trait(self, values: Union[Matrix,numpy.ndarray], raw = None, trait = None, **kwargs: dict):
        """
        Append values to the Matrix along the trait axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        trait : numpy.ndarray
            Taxa names to append to the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.append(
            values = values,
            axis = 0,
            raw = raw,
            trait = trait,
            **kwargs
        )

    def remove(self, obj: Union[int,slice,Sequence], axis, **kwargs: dict):
        """
        Remove sub-arrays along an axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
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
            # reset metadata
            self._taxa_grp_len = None
            self._taxa_grp_name = None
            self._taxa_grp_stix = None
            self._taxa_grp_spix = None
        if axis == 1:   # trait axis
            if self._raw is not None:
                self._raw = numpy.delete(self._raw, obj, axis = 2)
            if self._trait is not None:
                self._trait = numpy.delete(self._trait, obj, axis = 0)

    def remove_taxa(self, obj: Union[int,slice,Sequence], **kwargs: dict):
        """
        Remove sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.remove(
            obj = obj,
            axis = 0,
            **kwargs
        )

    def remove_trait(self, obj: Union[int,slice,Sequence], **kwargs: dict):
        """
        Remove sub-arrays along the trait axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.remove(
            obj = obj,
            axis = 1,
            **kwargs
        )

    def incorp(self, obj: Union[int,slice,Sequence], values, axis = -1, raw = None, taxa = None, taxa_grp = None, trait = None, **kwargs: dict):
        """
        Incorporate values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : array_like
            Values to incorporate into the matrix.
        axis : int
            The axis along which values are incorporated.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DenseEstimatedBreedingValueMatrix(values):
            if raw is None:
                raw = values.raw
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            if trait is None:
                trait = values.trait
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DensePhasedGenotypeVariantMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if axis == 0:
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot incorp: raw argument is required")
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot incorp: taxa_grp argument is required")
        elif axis == 1:
            if (self._raw is not None) and (raw is None):
                raise TypeError("cannot incorp: raw argument is required")
            if (self._trait is not None) and (trait is None):
                raise TypeError("cannot incorp: trait argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = axis)
        if axis == 0:
            if self._raw is not None:
                self._raw = numpy.insert(self._raw, obj, raw, axis = 1)
            if self._taxa is not None:
                self._taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
            if self._taxa_grp is not None:
                self._taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)
            # reset metadata
            self._taxa_grp_len = None
            self._taxa_grp_name = None
            self._taxa_grp_stix = None
            self._taxa_grp_spix = None
        elif axis == 1:
            if self._raw is not None:
                self._raw = numpy.insert(self._raw, obj, raw, axis = 2)
            if self._trait is not None:
                self._trait = numpy.insert(self._trait, obj, trait, axis = 0)

    def incorp_taxa(self, obj: Union[int,slice,Sequence], values, raw = None, taxa = None, taxa_grp = None, **kwargs: dict):
        """
        Incorporate values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        taxa : numpy.ndarray
            Taxa names to incorporate into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.incorp(
            obj = obj,
            values = values,
            axis = 0,
            raw = raw,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def incorp_trait(self, obj, values, raw = None, trait = None, **kwargs: dict):
        """
        Incorporate values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        trait : numpy.ndarray
            Taxa names to incorporate into the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.incorp(
            obj = obj,
            values = values,
            axis = 1,
            raw = raw,
            trait = trait,
            **kwargs
        )

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
        axis = get_axis(axis, self._mat.ndim)                   # transform axis number to an index
        emess = None                                            # error message

        if keys is None:                                        # if no keys were provided, set a default
            if axis == 0:                                       # taxa axis
                keys = (self._taxa, self._taxa_grp)             # taxa default keys
                emess = "taxa, taxa_grp are None"               # taxa error message
            elif axis == 1:                                     # trait axis
                keys = (self._trait)                            # trait default keys
                emess = "trait is None"                         # loci error message
            keys = tuple(k for k in keys if k is not None)      # remove None keys
            if len(keys) == 0:                                  # raise error if needed
                raise TypeError("cannot lexsort on axis {0}: {1}".format(axis, emess))
        else:
            l = self._mat.shape[axis]
            for i,k in enumerate(keys):
                if len(k) != l:
                    raise TypeError("cannot lexsort on axis %s: key %s is incompatible with axis length %s" % (axis, i, l))

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def lexsort_taxa(self, keys: Union[tuple,numpy.ndarray,None] = None, **kwargs: dict):
        """
        Perform an indirect stable sort using a sequence of keys along the taxa
        axis.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        return self.lexsort(
            keys = keys,
            axis = 0,
            **kwargs
        )

    def lexsort_trait(self, keys, **kwargs: dict):
        """
        Perform an indirect stable sort using a sequence of keys along the trait
        axis.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        return self.lexsort(
            keys = keys,
            axis = 1,
            **kwargs
        )

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
        # TODO: move reset of sort metadata to here
        # transform axis number to an index
        axis = get_axis(axis, self._mat.ndim)

        ########################################################################
        if axis == 0:                                           ### TAXA AXIS
            self._mat = self._mat[indices,:]                    # reorder mat array
            if self._raw is not None:
                self._raw = self._raw[:,indices,:]              # reorder raw array
            if self._taxa is not None:
                self._taxa = self._taxa[indices]                # reorder taxa array
            if self._taxa_grp is not None:
                self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array
        ########################################################################
        elif axis == 1:                                         ### LOCUS AXIS
            self._mat = self._mat[:,indices]                    # reorder mat array
            if self._raw is not None:
                self._raw = self._raw[:,:,indices]              # reorder raw array
            if self._trait is not None:
                self._trait = self._trait[indices]              # reorder trait array

    def reorder_taxa(self, indices, **kwargs: dict):
        """
        Reorder elements of the Matrix along the taxa axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.reorder(
            indices = indices,
            axis = 0,
            **kwargs
        )

    def reorder_trait(self, indices, **kwargs: dict):
        """
        Reorder elements of the Matrix along the trait axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.reorder(
            indices = indices,
            axis = 1,
            **kwargs
        )

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

        # TODO: move this section to reorder function
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

    def sort_taxa(self, keys: Union[tuple,numpy.ndarray,None] = None, **kwargs: dict):
        """
        Sort slements of the Matrix along the taxa axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.
        """
        self.sort(
            keys = keys,
            axis = 0,
            **kwargs
        )

    def sort_trait(self, keys, **kwargs: dict):
        """
        Sort slements of the Matrix along the trait axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.
        """
        self.sort(
            keys = keys,
            axis = 1,
            **kwargs
        )

    ################### Grouping Methods ###################
    def group(self, axis = -1):
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        # get axis index
        axis = get_axis(axis, self._mat.ndim)

        if axis == 1:
            raise ValueError("cannot group along axis 1 (trait axis): not groupable")

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

    def group_taxa(self, **kwargs: dict):
        """
        Sort the Matrix along the taxa axis, then populate grouping indices for
        the taxa axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.
        """
        self.group(
            axis = 0,
            **kwargs
        )

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

        if axis == 0:
            return (
                (self._taxa_grp_name is not None) and
                (self._taxa_grp_stix is not None) and
                (self._taxa_grp_spix is not None) and
                (self._taxa_grp_len is not None)
            )
        return False

    def is_grouped_taxa(self, **kwargs: dict):
        """
        Determine whether the Matrix has been sorted and grouped along the taxa
        axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        return self.is_grouped(
            axis = 0,
            **kwargs
        )

    ################### Matrix File I/O ####################
    @staticmethod
    def from_hdf5(filename, groupname = None):
        """
        Read GenotypeMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is read from base HDF5 group.

        Returns
        -------
        gmat : GenotypeMatrix
            A genotype matrix read from file.
        """
        check_file_exists(filename)                             # check file exists
        h5file = h5py.File(filename, "r")                       # open HDF5 in read only
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            check_group_in_hdf5(groupname, h5file, filename)    # check that group exists
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### check that we have all required fields
        required_fields = ["mat"]                               # all required arguments
        for field in required_fields:                           # for each required field
            fieldname = groupname + field                       # concatenate base groupname and field
            check_group_in_hdf5(fieldname, h5file, filename)    # check that group exists
        ######################################################### read data
        data_dict = {                                           # output dictionary
            "mat": None,
            "raw": None,
            "taxa": None,
            "taxa_grp": None,
            "trait": None
        }
        for field in data_dict.keys():                          # for each field
            fieldname = groupname + field                       # concatenate base groupname and field
            if fieldname in h5file:                             # if the field exists in the HDF5 file
                data_dict[field] = h5file[fieldname][:]         # read array
        ######################################################### read conclusion
        h5file.close()                                          # close file
        data_dict["taxa"] = numpy.object_(                      # convert taxa strings from byte to utf-8
            [s.decode("utf-8") for s in data_dict["taxa"]]
        )
        data_dict["trait"] = numpy.object_(                     # convert trait string from byte to utf-8
            [s.decode("utf-8") for s in data_dict["trait"]]
        )
        ######################################################### create object
        debvmat = DenseEstimatedBreedingValueMatrix(**data_dict)# create object from read data
        return debvmat

    def to_hdf5(self, filename, groupname = None):
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is written to the base HDF5 group.
        """
        h5file = h5py.File(filename, "a")                       # open HDF5 in write mode
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### populate HDF5 file
        data_dict = {                                           # data dictionary
            "mat": self.mat,
            "raw": self.raw,
            "taxa": self.taxa,
            "taxa_grp": self.taxa_grp,
            "trait": self.trait
        }
        save_dict_to_hdf5(h5file, groupname, data_dict)         # save data
        ######################################################### write conclusion
        h5file.close()                                          # close the file



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
