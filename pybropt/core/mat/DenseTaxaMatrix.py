import copy
import numpy

from . import DenseMutableMatrix
from . import TaxaMatrix

from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_iterable
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_dtype_is_int8
from pybropt.core.error import check_ndarray_is_2d
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_axis_len
from pybropt.core.error import cond_check_ndarray_dtype_is_bool
from pybropt.core.error import cond_check_ndarray_dtype_is_int64
from pybropt.core.error import cond_check_ndarray_dtype_is_float64
from pybropt.core.error import cond_check_ndarray_dtype_is_object
from pybropt.core.error import cond_check_ndarray_ndim
from pybropt.core.error import error_readonly
from pybropt.core.error import generic_check_isinstance
from pybropt.core.mat import get_axis

class DenseTaxaMatrix(DenseMutableMatrix,TaxaMatrix):
    """docstring for DenseTaxaMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, **kwargs):
        super(DenseTaxaMatrix, self).__init__(
            mat = mat,
            **kwargs
        )
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

    #################### Matrix copying ####################
    def __copy__(self):
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp)
        )

        # copy taxa metadata
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
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Taxa Data Properites #################
    def taxa():
        doc = "Taxa label property."
        def fget(self):
            """Get taxa label array"""
            return self._taxa
        def fset(self, value):
            """Set taxa label array"""
            cond_check_is_ndarray(value, "taxa")
            cond_check_ndarray_dtype_is_object(value, "taxa")
            cond_check_ndarray_ndim(value, "taxa", 1)
            cond_check_ndarray_axis_len(value, "taxa", 0, self.ntaxa)
            self._taxa = value
        def fdel(self):
            """Delete taxa label array"""
            del self._taxa
        return locals()
    taxa = property(**taxa())

    def taxa_grp():
        doc = "Taxa group label property."
        def fget(self):
            """Get taxa group label array"""
            return self._taxa_grp
        def fset(self, value):
            """Set taxa group label array"""
            cond_check_is_ndarray(value, "taxa_grp")
            cond_check_ndarray_dtype_is_int64(value, "taxa_grp")
            cond_check_ndarray_ndim(value, "taxa_grp", 1)
            cond_check_ndarray_axis_len(value, "taxa_grp", 0, self.ntaxa)
            self._taxa_grp = value
        def fdel(self):
            """Delete taxa group label array"""
            del self._taxa_grp
        return locals()
    taxa_grp = property(**taxa_grp())

    ############### Taxa Metadata Properites ###############
    def ntaxa():
        doc = "Number of taxa property."
        def fget(self):
            """Get number of taxa"""
            return self._mat.shape[self.taxa_axis]
        def fset(self, value):
            """Set number of taxa"""
            error_readonly("ntaxa")
        def fdel(self):
            """Delete number of taxa"""
            error_readonly("ntaxa")
        return locals()
    ntaxa = property(**ntaxa())

    def taxa_axis():
        doc = "Axis along which taxa are stored property."
        def fget(self):
            """Get taxa axis number"""
            return 0
        def fset(self, value):
            """Set taxa axis number"""
            error_readonly("taxa_axis")
        def fdel(self):
            """Delete taxa axis number"""
            error_readonly("taxa_axis")
        return locals()
    taxa_axis = property(**taxa_axis())

    def taxa_grp_name():
        doc = "Taxa group name property."
        def fget(self):
            """Get taxa group name array"""
            return self._taxa_grp_name
        def fset(self, value):
            """Set taxa group name array"""
            cond_check_is_ndarray(value, "taxa_grp_name")
            cond_check_ndarray_dtype_is_int64(value, "taxa_grp_name")
            cond_check_ndarray_ndim(value, "taxa_grp_name", 1)
            self._taxa_grp_name = value
        def fdel(self):
            """Delete taxa group array"""
            del self._taxa_grp_name
        return locals()
    taxa_grp_name = property(**taxa_grp_name())

    def taxa_grp_stix():
        doc = "Taxa group start index property."
        def fget(self):
            """Get taxa group start indices array"""
            return self._taxa_grp_stix
        def fset(self, value):
            """Set taxa group start indices array"""
            cond_check_is_ndarray(value, "taxa_grp_stix")
            cond_check_ndarray_dtype_is_int64(value, "taxa_grp_stix")
            cond_check_ndarray_ndim(value, "taxa_grp_stix", 1)
            self._taxa_grp_stix = value
        def fdel(self):
            """Delete taxa group start indices array"""
            del self._taxa_grp_stix
        return locals()
    taxa_grp_stix = property(**taxa_grp_stix())

    def taxa_grp_spix():
        doc = "Taxa group stop index property."
        def fget(self):
            """Get taxa group stop indices array"""
            return self._taxa_grp_spix
        def fset(self, value):
            """Set taxa group stop indices array"""
            cond_check_is_ndarray(value, "taxa_grp_spix")
            cond_check_ndarray_dtype_is_int64(value, "taxa_grp_spix")
            cond_check_ndarray_ndim(value, "taxa_grp_spix", 1)
            self._taxa_grp_spix = value
        def fdel(self):
            """Delete taxa group stop indices array"""
            del self._taxa_grp_spix
        return locals()
    taxa_grp_spix = property(**taxa_grp_spix())

    def taxa_grp_len():
        doc = "Taxa group length property."
        def fget(self):
            """Get taxa group length array"""
            return self._taxa_grp_len
        def fset(self, value):
            """Set taxa group length array"""
            cond_check_is_ndarray(value, "taxa_grp_len")
            cond_check_ndarray_dtype_is_int64(value, "taxa_grp_len")
            cond_check_ndarray_ndim(value, "taxa_grp_len", 1)
            self._taxa_grp_len = value
        def fdel(self):
            """Delete taxa group length array"""
            del self._taxa_grp_len
        return locals()
    taxa_grp_len = property(**taxa_grp_len())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values, axis = -1, taxa = None, taxa_grp = None, **kwargs):
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DensePhasedGenotypeMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedGenotypeMatrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.adjoin_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

        return out

    def adjoin_taxa(self, values, taxa = None, taxa_grp = None, **kwargs):
        """
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        **kwargs
            Additional keyword arguments.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("cannot adjoin: 'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")   # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot adjoin: 'taxa_grp' argument is required")

        # adjoin values
        values = numpy.append(self._mat, values, axis = self.taxa_axis)
        if self._taxa is not None:
            taxa = numpy.append(self._taxa, taxa, axis = 0)
        if self._taxa_grp is not None:
            taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)

        out = self.__class__(
            mat = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    def delete(self, obj, axis = -1, **kwargs):
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

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.delete_taxa(
                obj = obj,
                **kwargs
            )
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_taxa(self, obj, **kwargs):
        """
        Delete sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp

        # delete values
        mat = numpy.delete(mat, obj, axis = self.taxa_axis)
        if taxa is not None:
            taxa = numpy.delete(taxa, obj, axis = 0)
        if taxa_grp is not None:
            taxa_grp = numpy.delete(taxa_grp, obj, axis = 0)

        out = self.__class__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    def insert(self, obj, values, axis = -1, taxa = None, taxa_grp = None, **kwargs):
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : DenseHaplotypeMatrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : DenseHaplotypeMatrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.insert_taxa(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    def insert_taxa(self, obj, values, taxa = None, taxa_grp = None, **kwargs):
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
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot insert: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")   # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot insert: 'taxa_grp' argument is required")

        # insert values
        values = numpy.insert(self._mat, obj, values, axis = self.taxa_axis)
        if self._taxa is not None:
            taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
        if self._taxa_grp is not None:
            taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)

        # create output
        out = self.__class__(
            mat = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    def select(self, indices, axis = -1, **kwargs):
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
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.select_taxa(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_taxa(self, indices, **kwargs):
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
        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp

        # select values
        mat = numpy.take(mat, indices, axis = self.taxa_axis)
        if taxa is not None:
            taxa = numpy.take(taxa, indices, axis = 0)
        if taxa_grp is not None:
            taxa_grp = numpy.take(taxa_grp, indices, axis = 0)

        out = self.__class__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    @classmethod
    def concat(cls, mats, axis = -1, **kwargs):
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : array_like of matrices
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : DenseHaplotypeMatrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].taxa_axis:
            out = cls.concat_taxa(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_taxa(cls, mats, **kwargs):
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
        # ensure that we have an iterable object
        check_is_iterable(mats, "mats")

        # ensure that we have an array_like of length >= 1
        if len(mats) <= 0:
            raise ValueError("need at least one Matrix to concatenate")

        # ensure that all items in mats are DenseHaplotypeMatrix
        for i,v in enumerate(mats):
            generic_check_isinstance(v, "mats[{0}]".format(i), cls)

        # make sure dimensions are all identical to first element in mats
        if any(m.mat_ndim != mats[0].mat_ndim for m in mats):
            raise ValueError("cannot concat: not all matrices have the same number of dimensions")

        # extract tuple of shapes for testing compatibility
        shape_t = tuple(zip(*[m.mat.shape for m in mats]))

        # test matrix compatibility (same axis length along non-taxa axes)
        for i,v in enumerate(shape_t):                              # for each index,tuple in shape_t
            if (i != mats[0].taxa_axis) and any(l != v[0] for l in v): # if not the taxa axis AND axis lengths are different
                raise ValueError("cannot concat: matrix shapes do not all align along axis {0}".format(i))

        # create matrix lists
        mat_ls = [m.mat for m in mats]
        taxa_ls = [m.taxa for m in mats]
        taxa_grp_ls = [m.taxa_grp for m in mats]

        # process/error check taxa_ls
        if all(e is None for e in taxa_ls):                             # if all elements are None
            taxa_ls = None                                              # replace list with None
        else:                                                           # else at least one matrix does not have a taxa array
            for i,v in enumerate(taxa_ls):                              # for each index,element in taxa_ls
                if v is None:                                           # if element is None
                    ntaxa = shape_t[mats[0].taxa_axis][i]               # get number of taxa
                    taxa_ls[i] = numpy.empty(ntaxa, dtype = "object")   # replace with array of None

        # process/error check taxa_grp_ls
        if all(e is None for e in taxa_grp_ls):         # if all elements are None
            taxa_grp_ls = None                          # replace list with None
        elif any(e is None for e in taxa_grp_ls):       # else if any elements are None
            raise ValueError("cannot concat: 'taxa_grp' needed for all Matrix in list")

        # concatenate mat, taxa, taxa_grp items
        mat = numpy.concatenate(mat_ls, axis = mats[0].taxa_axis)
        taxa = None if taxa_ls is None else numpy.concatenate(taxa_ls, axis = 0)
        taxa_grp = None if taxa_grp_ls is None else numpy.concatenate(taxa_grp_ls, axis = 0)

        # TODO: decide if first element in list is good source of information
        # concatenate everything and put into new DenseHaplotypeMatrix
        # use first element as source of variant data
        out = cls(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis = -1, taxa = None, taxa_grp = None, **kwargs):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DenseHaplotypeMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.append_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

    def append_taxa(self, values, taxa = None, taxa_grp = None, **kwargs):
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
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot append: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object") # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot append: 'taxa_grp' argument is required")

        # append values
        self._mat = numpy.append(self._mat, values, axis = self.taxa_axis)
        if self._taxa is not None:
            self._taxa = numpy.append(self._taxa, taxa, axis = 0)
        if self._taxa_grp is not None:
            self._taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)

        # reset metadata
        self._taxa_grp_len = None
        self._taxa_grp_name = None
        self._taxa_grp_stix = None
        self._taxa_grp_spix = None

    def remove(self, obj, axis = -1, **kwargs):
        """
        Remove sub-arrays along an axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to remove the subarray defined by obj.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.taxa_axis:
            self.remove_taxa(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def remove_taxa(self, obj, **kwargs):
        """
        Remove sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        # delete values
        self._mat = numpy.delete(self._mat, obj, axis = self.taxa_axis)

        if self._taxa is not None:
            self._taxa = numpy.delete(self._taxa, obj, axis = 0)
        if self._taxa_grp is not None:
            self._taxa_grp = numpy.delete(self._taxa_grp, obj, axis = 0)

        # reset metadata
        self._taxa_grp_len = None
        self._taxa_grp_name = None
        self._taxa_grp_stix = None
        self._taxa_grp_spix = None

    def incorp(self, obj, values, axis = -1, taxa = None, taxa_grp = None, **kwargs):
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
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.taxa_axis:
            self.incorp(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    def incorp_taxa(self, obj, values, taxa = None, taxa_grp = None, **kwargs):
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
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseHaplotypeMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")  # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot incorp: 'taxa_grp' argument is required")

        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = self.taxa_axis)

        if self._taxa is not None:
            self._taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
        if self._taxa_grp is not None:
            self._taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)

        # reset metadata
        self._taxa_grp_len = None
        self._taxa_grp_name = None
        self._taxa_grp_stix = None
        self._taxa_grp_spix = None

    ################### Sorting Methods ####################
    def lexsort(self, keys = None, axis = -1, **kwargs):
        """
        Perform an indirect stable sort using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key.
        axis : int
            The axis of the Matrix over which to sort values.

        Returns
        -------
        indices : numpy.ndarray
            Array of indices that sort the keys.
        """
        axis = get_axis(axis, self.mat_ndim)    # transform axis number to an index
        indices = None                          # declare variable

        # dispatch to correct function
        if axis == self.taxa_axis:
            self.lexsort_taxa(keys = keys, **kwargs)
        else:
            raise ValueError("cannot lexsort along axis {0}".format(axis))

        return indices

    def lexsort_taxa(self, keys = None, **kwargs):
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
        # default error message
        emess = "no available keys to sort"

        # if no keys were provided, set a default
        if keys is None:
            keys = (self._taxa, self._taxa_grp)         # taxa default keys
            emess = "taxa, taxa_grp are None"           # taxa error message

        # remove None keys
        keys = tuple(k for k in keys if k is not None)

        # raise error if no keys remain
        if len(keys) == 0:
            raise ValueError("cannot lexsort on axis {0} (taxa axis): {1}".format(self.taxa_axis, emess))

        # raise error if keys are of incompatible length
        if any(len(k) != self.ntaxa for k in keys):
            emess = "keys are not all length {0}".format(self.ntaxa)
            raise ValueError("cannot lexsort on axis {0} (taxa axis): {1}".format(self.taxa_axis, emess))

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def reorder(self, indices, axis = -1, **kwargs):
        """
        Reorder the VariantMatrix.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        axis : int
            The axis over which to reorder values.

        """
        axis = get_axis(axis, self.mat_ndim)                   # transform axis number to an index

        if axis == self.taxa_axis:
            self.reorder(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    def reorder_taxa(self, indices, **kwargs):
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
        # build a tuple to slice the matrix
        ix = tuple(indices if i == self.taxa_axis else slice(None) for i in range(self.mat_ndim))

        # reorder arrays
        self._mat = self._mat[ix]

        if self._taxa is not None:
            self._taxa = self._taxa[indices]                # reorder taxa array
        if self._taxa_grp is not None:
            self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array

    def sort(self, keys = None, axis = -1, **kwargs):
        """
        Reset metadata for corresponding axis: name, stix, spix, len.
        Sort the VariantMatrix using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key.
        axis : int
            The axis over which to sort values.
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.sort_taxa(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    def sort_taxa(self, keys = None, **kwargs):
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
        # reset taxa group metadata
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

        # get indices for sort
        indices = self.lexsort_taxa(keys, **kwargs)

        # reorder internals
        self.reorder_taxa(indices, **kwargs)

    ################### Grouping Methods ###################
    def group(self, axis = -1, **kwargs):
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.group_taxa(**kwargs)
        else:
            raise ValueError("cannot group along axis {0}".format(axis))

    def group_taxa(self, **kwargs):
        """
        Sort the Matrix along the taxa axis, then populate grouping indices for
        the taxa axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.
        """
        # sort taxa using default keys
        self.sort_taxa()

        if self._taxa_grp is not None:
            # get unique taxa group names, starting indices, group lengths
            uniq = numpy.unique(self._taxa_grp, return_index = True, return_counts = True)
            # make assignments to instance data
            self._taxa_grp_name, self._taxa_grp_stix, self._taxa_grp_len = uniq
            # calculate stop indices
            self._taxa_grp_spix = self._taxa_grp_stix + self._taxa_grp_len

    def is_grouped(self, axis = -1, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        axis = get_axis(axis, self.mat_ndim)    # transform axis number to an index
        grouped = False                         # default output

        if axis == self.taxa_axis:
            grouped = self.is_grouped_taxa(**kwargs)
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped

    def is_grouped_taxa(self, **kwargs):
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
        return (
            (self._taxa_grp_name is not None) and
            (self._taxa_grp_stix is not None) and
            (self._taxa_grp_spix is not None) and
            (self._taxa_grp_len is not None)
        )



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseTaxaMatrix(v):
    """
    Determine whether an object is a DenseTaxaMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseTaxaMatrix object instance.
    """
    return isinstance(v, DenseTaxaMatrix)

def check_is_DenseTaxaMatrix(v, varname):
    """
    Check if object is of type DenseTaxaMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_DenseTaxaMatrix(v):
        raise TypeError("'{0}' must be a DenseTaxaMatrix".format(varname))

def cond_check_is_DenseTaxaMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DenseTaxaMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        DenseTaxaMatrix.
    """
    if cond(v):
        check_is_DenseTaxaMatrix(v, varname)