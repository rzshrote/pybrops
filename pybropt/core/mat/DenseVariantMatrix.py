import copy
import numpy

from pybropt.core.mat.DenseMutableMatrix import DenseMutableMatrix
from pybropt.core.mat.VariantMatrix import VariantMatrix

from pybropt.core.error import check_is_int
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
from pybropt.core.mat.util import get_axis
from pybropt.core.error import check_is_iterable

class DenseVariantMatrix(DenseMutableMatrix,VariantMatrix):
    """docstring for DenseVariantMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Parameters
        ----------
        mat : numpy.ndarray
            An int8 haplotype matrix. Must be {0,1,2} format.
        ploidy : int
            The ploidy represented by the haplotype matrix. This only represents
            ploidy of the reproductive habit. If the organism represented is an
            allopolyploid (e.g. hexaploid wheat), the ploidy is 2 since it
            reproduces in a diploid manner.
        # TODO: Add a mat_format option to store as {0,1,2}, {-1,0,1}, etc.
        """
        super(DenseVariantMatrix, self).__init__(
            mat = mat,
            **kwargs
        )
        # error check and setting values using properties
        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_name = vrnt_name
        self.vrnt_genpos = vrnt_genpos
        self.vrnt_xoprob = vrnt_xoprob
        self.vrnt_hapgrp = vrnt_hapgrp
        self.vrnt_hapalt = vrnt_hapalt
        self.vrnt_hapref = vrnt_hapref
        self.vrnt_mask = vrnt_mask

        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len = None

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
            vrnt_chrgrp = copy.copy(self.vrnt_chrgrp),
            vrnt_phypos = copy.copy(self.vrnt_phypos),
            vrnt_name = copy.copy(self.vrnt_name),
            vrnt_genpos = copy.copy(self.vrnt_genpos),
            vrnt_xoprob = copy.copy(self.vrnt_xoprob),
            vrnt_hapgrp = copy.copy(self.vrnt_hapgrp),
            vrnt_hapalt = copy.copy(self.vrnt_hapalt),
            vrnt_hapref = copy.copy(self.vrnt_hapref),
            vrnt_mask = copy.copy(self.vrnt_mask)
        )

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.copy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.copy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.copy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len = copy.copy(self.vrnt_chrgrp_len)

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
            vrnt_chrgrp = copy.deepcopy(self.vrnt_chrgrp, memo),
            vrnt_phypos = copy.deepcopy(self.vrnt_phypos, memo),
            vrnt_name = copy.deepcopy(self.vrnt_name, memo),
            vrnt_genpos = copy.deepcopy(self.vrnt_genpos, memo),
            vrnt_xoprob = copy.deepcopy(self.vrnt_xoprob, memo),
            vrnt_hapgrp = copy.deepcopy(self.vrnt_hapgrp, memo),
            vrnt_hapalt = copy.deepcopy(self.vrnt_hapalt, memo),
            vrnt_hapref = copy.deepcopy(self.vrnt_hapref, memo),
            vrnt_mask = copy.deepcopy(self.vrnt_mask, memo)
        )

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name, memo)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix, memo)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix, memo)
        out.vrnt_chrgrp_len = copy.deepcopy(self.vrnt_chrgrp_len, memo)

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Variant Data Properites ################
    def vrnt_chrgrp():
        doc = "Variant chromosome group label property."
        def fget(self):
            """Get variant chromosome group lable array"""
            return self._vrnt_chrgrp
        def fset(self, value):
            """Set variant chromosome group lable array"""
            cond_check_is_ndarray(value, "vrnt_chrgrp")
            cond_check_ndarray_dtype_is_int64(value, "vrnt_chrgrp")
            cond_check_ndarray_ndim(value, "vrnt_chrgrp", 1)
            cond_check_ndarray_axis_len(value, "vrnt_chrgrp", 0, self.nvrnt)
            self._vrnt_chrgrp = value
        def fdel(self):
            """Delete variant chromosome group lable array"""
            del self._vrnt_chrgrp
        return locals()
    vrnt_chrgrp = property(**vrnt_chrgrp())

    def vrnt_phypos():
        doc = "Variant physical position property."
        def fget(self):
            """Get variant physical position array"""
            return self._vrnt_phypos
        def fset(self, value):
            """Set variant physical position array"""
            cond_check_is_ndarray(value, "vrnt_phypos")
            cond_check_ndarray_dtype_is_int64(value, "vrnt_phypos")
            cond_check_ndarray_ndim(value, "vrnt_phypos", 1)
            cond_check_ndarray_axis_len(value, "vrnt_phypos", 0, self.nvrnt)
            self._vrnt_phypos = value
        def fdel(self):
            """Delete variant physical position array"""
            del self._vrnt_phypos
        return locals()
    vrnt_phypos = property(**vrnt_phypos())

    def vrnt_name():
        doc = "Variant name property."
        def fget(self):
            """Get variant name array"""
            return self._vrnt_name
        def fset(self, value):
            """Set variant name array"""
            cond_check_is_ndarray(value, "vrnt_name")
            cond_check_ndarray_dtype_is_object(value, "vrnt_name")
            cond_check_ndarray_ndim(value, "vrnt_name", 1)
            cond_check_ndarray_axis_len(value, "vrnt_name", 0, self.nvrnt)
            self._vrnt_name = value
        def fdel(self):
            """Delete variant name array"""
            del self._vrnt_name
        return locals()
    vrnt_name = property(**vrnt_name())

    def vrnt_genpos():
        doc = "Variant genetic position property."
        def fget(self):
            """Get variant genetic position array"""
            return self._vrnt_genpos
        def fset(self, value):
            """Set variant genetic position array"""
            cond_check_is_ndarray(value, "vrnt_genpos")
            cond_check_ndarray_dtype_is_float64(value, "vrnt_genpos")
            cond_check_ndarray_ndim(value, "vrnt_genpos", 1)
            cond_check_ndarray_axis_len(value, "vrnt_genpos", 0, self.nvrnt)
            self._vrnt_genpos = value
        def fdel(self):
            """Delete variant genetic position array"""
            del self._vrnt_genpos
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    def vrnt_xoprob():
        doc = "Variant crossover sequential probability property."
        def fget(self):
            """Get variant crossover sequential probability array"""
            return self._vrnt_xoprob
        def fset(self, value):
            """Set variant crossover sequential probability array"""
            cond_check_is_ndarray(value, "vrnt_xoprob")
            cond_check_ndarray_dtype_is_float64(value, "vrnt_xoprob")
            cond_check_ndarray_ndim(value, "vrnt_xoprob", 1)
            cond_check_ndarray_axis_len(value, "vrnt_xoprob", 0, self.nvrnt)
            self._vrnt_xoprob = value
        def fdel(self):
            """Delete variant crossover sequential probability array"""
            del self._vrnt_xoprob
        return locals()
    vrnt_xoprob = property(**vrnt_xoprob())

    def vrnt_hapgrp():
        doc = "Variant haplotype group label property."
        def fget(self):
            """Get variant haplotype group label array"""
            return self._vrnt_hapgrp
        def fset(self, value):
            """Set variant haplotype group label array"""
            cond_check_is_ndarray(value, "vrnt_hapgrp")
            cond_check_ndarray_dtype_is_int64(value, "vrnt_hapgrp")
            cond_check_ndarray_ndim(value, "vrnt_hapgrp", 1)
            cond_check_ndarray_axis_len(value, "vrnt_hapgrp", 0, self.nvrnt)
            self._vrnt_hapgrp = value
        def fdel(self):
            """Delete variant haplotype group label array"""
            del self._vrnt_hapgrp
        return locals()
    vrnt_hapgrp = property(**vrnt_hapgrp())

    def vrnt_hapalt():
        doc = "Variant haplotype sequence property."
        def fget(self):
            """Get variant haplotype sequence"""
            return self._vrnt_hapalt
        def fset(self, value):
            """Set variant haplotype sequence"""
            cond_check_is_ndarray(value, "vrnt_hapalt")
            cond_check_ndarray_dtype_is_object(value, "vrnt_hapalt")
            cond_check_ndarray_ndim(value, "vrnt_hapalt", 1)
            cond_check_ndarray_axis_len(value, "vrnt_hapalt", 0, self.nvrnt)
            self._vrnt_hapalt = value
        def fdel(self):
            """Delete variant haplotype sequence"""
            del self._vrnt_hapalt
        return locals()
    vrnt_hapalt = property(**vrnt_hapalt())

    def vrnt_hapref():
        doc = "Variant reference haplotype sequence property."
        def fget(self):
            """Get variant reference haplotype sequence"""
            return self._vrnt_hapref
        def fset(self, value):
            """Set variant reference haplotype sequence"""
            cond_check_is_ndarray(value, "vrnt_hapref")
            cond_check_ndarray_dtype_is_object(value, "vrnt_hapref")
            cond_check_ndarray_ndim(value, "vrnt_hapref", 1)
            cond_check_ndarray_axis_len(value, "vrnt_hapref", 0, self.nvrnt)
            self._vrnt_hapref = value
        def fdel(self):
            """Delete variant reference haplotype sequence"""
            del self._vrnt_hapref
        return locals()
    vrnt_hapref = property(**vrnt_hapref())

    def vrnt_mask():
        doc = "Variant mask property."
        def fget(self):
            """Get variant mask"""
            return self._vrnt_mask
        def fset(self, value):
            """Set variant mask"""
            cond_check_is_ndarray(value, "vrnt_mask")
            cond_check_ndarray_dtype_is_bool(value, "vrnt_mask")
            cond_check_ndarray_ndim(value, "vrnt_mask", 1)
            cond_check_ndarray_axis_len(value, "vrnt_mask", 0, self.nvrnt)
            self._vrnt_mask = value
        def fdel(self):
            """Delete variant mask"""
            del self._vrnt_mask
        return locals()
    vrnt_mask = property(**vrnt_mask())

    ############# Variant Metadata Properites ##############
    def nvrnt():
        doc = "Number of variants property."
        def fget(self):
            """Get number of variants"""
            return self._mat.shape[self.vrnt_axis]
        def fset(self, value):
            """Set number of variants"""
            error_readonly("nvrnt")
        def fdel(self):
            """Delete number of variants"""
            error_readonly("nvrnt")
        return locals()
    nvrnt = property(**nvrnt())

    def vrnt_axis():
        doc = "Axis along which variants are stored property."
        def fget(self):
            """Get variant axis"""
            return 0
        def fset(self, value):
            """Set variant axis"""
            error_readonly("vrnt_axis")
        def fdel(self):
            """Delete variant axis"""
            error_readonly("vrnt_axis")
        return locals()
    vrnt_axis = property(**vrnt_axis())

    def vrnt_chrgrp_name():
        doc = "Variant chromosome group names property."
        def fget(self):
            """Get variant chromosome group name array"""
            return self._vrnt_chrgrp_name
        def fset(self, value):
            """Set variant chromosome group name array"""
            cond_check_is_ndarray(value, "vrnt_chrgrp_name")
            cond_check_ndarray_dtype_is_int64(value, "vrnt_chrgrp_name")
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_name", 1)
            self._vrnt_chrgrp_name = value
        def fdel(self):
            """Delete variant chromosome group name array"""
            del self._vrnt_chrgrp_name
        return locals()
    vrnt_chrgrp_name = property(**vrnt_chrgrp_name())

    def vrnt_chrgrp_stix():
        doc = "Variant chromosome group start indices property."
        def fget(self):
            """Get variant chromosome group start indices array"""
            return self._vrnt_chrgrp_stix
        def fset(self, value):
            """Set variant chromosome group start indices array"""
            cond_check_is_ndarray(value, "vrnt_chrgrp_stix")
            cond_check_ndarray_dtype_is_int64(value, "vrnt_chrgrp_stix")
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_stix", 1)
            self._vrnt_chrgrp_stix = value
        def fdel(self):
            """Delete variant chromosome group start indices array"""
            del self._vrnt_chrgrp_stix
        return locals()
    vrnt_chrgrp_stix = property(**vrnt_chrgrp_stix())

    def vrnt_chrgrp_spix():
        doc = "Variant chromosome group stop indices property."
        def fget(self):
            """Get variant chromosome group stop indices array"""
            return self._vrnt_chrgrp_spix
        def fset(self, value):
            """Set variant chromosome group stop indices array"""
            cond_check_is_ndarray(value, "vrnt_chrgrp_spix")
            cond_check_ndarray_dtype_is_int64(value, "vrnt_chrgrp_spix")
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_spix", 1)
            self._vrnt_chrgrp_spix = value
        def fdel(self):
            """Delete variant chromosome group stop indices array"""
            del self._vrnt_chrgrp_spix
        return locals()
    vrnt_chrgrp_spix = property(**vrnt_chrgrp_spix())

    def vrnt_chrgrp_len():
        doc = "Variant chromosome group length property."
        def fget(self):
            """Get variant chromosome group length array"""
            return self._vrnt_chrgrp_len
        def fset(self, value):
            """Set variant chromosome group length array"""
            cond_check_is_ndarray(value, "vrnt_chrgrp_len")
            cond_check_ndarray_dtype_is_int64(value, "vrnt_chrgrp_len")
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_len", 1)
            self._vrnt_chrgrp_len = value
        def fdel(self):
            """Delete variant chromosome group length array"""
            del self._vrnt_chrgrp_len
        return locals()
    vrnt_chrgrp_len = property(**vrnt_chrgrp_len())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values, axis = -1, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DensePhasedGenotypeMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        kwargs : dict
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
        if axis == self.vrnt_axis:
            out = self.adjoin_vrnt(
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

        return out

    def adjoin_vrnt(self, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Add additional elements to the end of the Matrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
        vrnt_hapalt : numpy.ndarray
            Variant haplotype sequence.
        vrnt_hapref : numpy.ndarray
            Variant haplotype reference sequence.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if isinstance(values, self.__class__):
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_hapalt is None:
                vrnt_hapalt = values.vrnt_hapalt
            if vrnt_hapref is None:
                vrnt_hapref = values.vrnt_hapref
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("cannot adjoin: 'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.vrnt_axis) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))
        if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
            raise TypeError("cannot adjoin: 'vrnt_chrgrp' argument is required")
        if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
            raise TypeError("cannot adjoin: 'vrnt_phypos' argument is required")
        if (self._vrnt_name is not None) and (vrnt_name is None):
            vrnt_name = numpy.empty(values.shape[self.vrnt_axis], dtype = "object")  # fill with None
        if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
            raise TypeError("cannot adjoin: 'vrnt_genpos' argument is required")
        if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
            raise TypeError("cannot adjoin: 'vrnt_xoprob' argument is required")
        if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
            raise TypeError("cannot adjoin: 'vrnt_hapgrp' argument is required")
        if (self._vrnt_hapalt is not None) and (vrnt_hapalt is None):
            raise TypeError("cannot adjoin: 'vrnt_hapalt' argument is required")
        if (self._vrnt_hapref is not None) and (vrnt_hapref is None):
            raise TypeError("cannot adjoin: 'vrnt_hapref' argument is required")
        if (self._vrnt_mask is not None) and (vrnt_mask is None):
            raise TypeError("cannot adjoin: 'vrnt_mask' argument is required")

        # adjoin values
        values = numpy.append(self._mat, values, axis = self.vrnt_axis)
        if self._vrnt_chrgrp is not None:
            vrnt_chrgrp = numpy.append(self._vrnt_chrgrp, vrnt_chrgrp, axis = 0)
        if self._vrnt_phypos is not None:
            vrnt_phypos = numpy.append(self._vrnt_phypos, vrnt_phypos, axis = 0)
        if self._vrnt_name is not None:
            vrnt_name = numpy.append(self._vrnt_name, vrnt_name, axis = 0)
        if self._vrnt_genpos is not None:
            vrnt_genpos = numpy.append(self._vrnt_genpos, vrnt_genpos, axis = 0)
        if self._vrnt_xoprob is not None:
            vrnt_xoprob = numpy.append(self._vrnt_xoprob, vrnt_xoprob, axis = 0)
        if self._vrnt_hapgrp is not None:
            vrnt_hapgrp = numpy.append(self._vrnt_hapgrp, vrnt_hapgrp, axis = 0)
        if self._vrnt_hapalt is not None:
            vrnt_hapalt = numpy.append(self._vrnt_hapalt, vrnt_hapalt, axis = 0)
        if self._vrnt_hapref is not None:
            vrnt_hapref = numpy.append(self._vrnt_hapref, vrnt_hapref, axis = 0)
        if self._vrnt_mask is not None:
            vrnt_mask = numpy.append(self._vrnt_mask, vrnt_mask, axis = 0)

        out = self.__class__(
            mat = values,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
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
        kwargs : dict
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
        if axis == self.vrnt_axis:
            out = self.delete_vrnt(
                obj = obj,
                **kwargs
            )
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_vrnt(self, obj, **kwargs):
        """
        Delete sub-arrays along the variant axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get values
        mat = self._mat
        vrnt_chrgrp = self._vrnt_chrgrp
        vrnt_phypos = self._vrnt_phypos
        vrnt_name = self._vrnt_name
        vrnt_genpos = self._vrnt_genpos
        vrnt_xoprob = self._vrnt_xoprob
        vrnt_hapgrp = self._vrnt_hapgrp
        vrnt_hapalt = self._vrnt_hapalt
        vrnt_hapref = self._vrnt_hapref
        vrnt_mask = self._vrnt_mask

        # delete values
        mat = numpy.delete(mat, obj, axis = self.vrnt_axis)
        if vrnt_chrgrp is not None:
            vrnt_chrgrp = numpy.delete(vrnt_chrgrp, obj, axis = 0)
        if vrnt_phypos is not None:
            vrnt_phypos = numpy.delete(vrnt_phypos, obj, axis = 0)
        if vrnt_name is not None:
            vrnt_name = numpy.delete(vrnt_name, obj, axis = 0)
        if vrnt_genpos is not None:
            vrnt_genpos = numpy.delete(vrnt_genpos, obj, axis = 0)
        if vrnt_xoprob is not None:
            vrnt_xoprob = numpy.delete(vrnt_xoprob, obj, axis = 0)
        if vrnt_hapgrp is not None:
            vrnt_hapgrp = numpy.delete(vrnt_hapgrp, obj, axis = 0)
        if vrnt_hapalt is not None:
            vrnt_hapalt = numpy.delete(vrnt_hapalt, obj, axis = 0)
        if vrnt_hapref is not None:
            vrnt_hapref = numpy.delete(vrnt_hapref, obj, axis = 0)
        if vrnt_mask is not None:
            vrnt_mask = numpy.delete(vrnt_mask, obj, axis = 0)

        out = self.__class__(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        return out

    def insert(self, obj, values, axis = -1, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : DenseVariantMatrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
            If values is a DenseVariantMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseVariantMatrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.vrnt_axis:
            out = self.insert_vrnt(
                obj = obj,
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    def insert_vrnt(self, obj, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Insert values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
        vrnt_hapalt : numpy.ndarray
        vrnt_hapref : numpy.ndarray
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if isinstance(values, self.__class__):
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_hapalt is None:
                vrnt_hapalt = values.vrnt_hapalt
            if vrnt_hapref is None:
                vrnt_hapref = values.vrnt_hapref
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.vrnt_axis) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))
        if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
            raise TypeError("cannot insert: 'vrnt_chrgrp' argument is required")
        if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
            raise TypeError("cannot insert: 'vrnt_phypos' argument is required")
        if (self._vrnt_name is not None) and (vrnt_name is None):
            vrnt_name = numpy.empty(values.shape[self.vrnt_axis], dtype = "object") # fill with None
        if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
            raise TypeError("cannot insert: 'vrnt_genpos' argument is required")
        if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
            raise TypeError("cannot insert: vrnt_xoprob argument is required")
        if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
            raise TypeError("cannot insert: vrnt_hapgrp argument is required")
        if (self._vrnt_mask is not None) and (vrnt_mask is None):
            raise TypeError("cannot insert: vrnt_mask argument is required")

        # insert values
        values = numpy.insert(self._mat, obj, values, axis = self.vrnt_axis)
        if self._vrnt_chrgrp is not None:
            vrnt_chrgrp = numpy.insert(self._vrnt_chrgrp, obj, vrnt_chrgrp, axis = 0)
        if self._vrnt_phypos is not None:
            vrnt_phypos = numpy.insert(self._vrnt_phypos, obj, vrnt_phypos, axis = 0)
        if self._vrnt_name is not None:
            vrnt_name = numpy.insert(self._vrnt_name, obj, vrnt_name, axis = 0)
        if self._vrnt_genpos is not None:
            vrnt_genpos = numpy.insert(self._vrnt_genpos, obj, vrnt_genpos, axis = 0)
        if self._vrnt_xoprob is not None:
            vrnt_xoprob = numpy.insert(self._vrnt_xoprob, obj, vrnt_xoprob, axis = 0)
        if self._vrnt_hapgrp is not None:
            vrnt_hapgrp = numpy.insert(self._vrnt_hapgrp, obj, vrnt_hapgrp, axis = 0)
        if self._vrnt_hapalt is not None:
            vrnt_hapalt = numpy.insert(self._vrnt_hapalt, obj, vrnt_hapalt, axis = 0)
        if self._vrnt_hapref is not None:
            vrnt_hapref = numpy.insert(self._vrnt_hapref, obj, vrnt_hapref, axis = 0)
        if self._vrnt_mask is not None:
            vrnt_mask = numpy.insert(self._vrnt_mask, obj, vrnt_mask, axis = 0)

        # create output
        out = self.__class__(
            mat = values,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
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
        kwargs : dict
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
        if axis == self.vrnt_axis:
            out = self.select_vrnt(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_vrnt(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the variant axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        # get values
        mat = self._mat
        vrnt_chrgrp = self._vrnt_chrgrp
        vrnt_phypos = self._vrnt_phypos
        vrnt_name = self._vrnt_name
        vrnt_genpos = self._vrnt_genpos
        vrnt_xoprob = self._vrnt_xoprob
        vrnt_hapgrp = self._vrnt_hapgrp
        vrnt_hapalt = self._vrnt_hapalt
        vrnt_hapref = self._vrnt_hapref
        vrnt_mask = self._vrnt_mask

        # select values
        mat = numpy.take(mat, indices, axis = self.vrnt_axis)
        if vrnt_chrgrp is not None:
            vrnt_chrgrp = numpy.take(vrnt_chrgrp, indices, axis = 0)
        if vrnt_phypos is not None:
            vrnt_phypos = numpy.take(vrnt_phypos, indices, axis = 0)
        if vrnt_name is not None:
            vrnt_name = numpy.take(vrnt_name, indices, axis = 0)
        if vrnt_genpos is not None:
            vrnt_genpos = numpy.take(vrnt_genpos, indices, axis = 0)
        if vrnt_xoprob is not None:
            vrnt_xoprob = numpy.take(vrnt_xoprob, indices, axis = 0)
        if vrnt_hapgrp is not None:
            vrnt_hapgrp = numpy.take(vrnt_hapgrp, indices, axis = 0)
        if vrnt_hapalt is not None:
            vrnt_hapalt = numpy.take(vrnt_hapalt, indices, axis = 0)
        if vrnt_hapref is not None:
            vrnt_hapref = numpy.take(vrnt_hapref, indices, axis = 0)
        if vrnt_mask is not None:
            vrnt_mask = numpy.take(vrnt_mask, indices, axis = 0)

        out = self.__class__(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
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
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : DenseVariantMatrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].vrnt_axis:
            out = cls.concat_vrnt(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_vrnt(cls, mats, **kwargs):
        """
        Concatenate list of Matrix together along the variant axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
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

        # ensure that all items in mats are DenseVariantMatrix
        for i,v in enumerate(mats):
            generic_check_isinstance(v, "mats[{0}]".format(i), cls)

        # make sure dimensions are all identical to first element in mats
        if any(m.mat_ndim != mats[0].mat_ndim for m in mats):
            raise ValueError("cannot concat: not all matrices have the same number of dimensions")

        # extract tuple of shapes for testing compatibility
        shape_t = tuple(zip(*[m.mat.shape for m in mats]))

        # test matrix compatibility (same axis length along non-variant axes)
        for i,v in enumerate(shape_t):                              # for each index,tuple in shape_t
            if (i != mats[0].vrnt_axis) and any(l != v[0] for l in v): # if not the variant axis AND axis lengths are different
                raise ValueError("cannot concat: matrix shapes do not all align along axis {0}".format(i))

        # create matrix lists
        mat_ls = [m.mat for m in mats]
        vrnt_chrgrp_ls = [m.vrnt_chrgrp for m in mats]
        vrnt_phypos_ls = [m.vrnt_phypos for m in mats]
        vrnt_name_ls = [m.vrnt_name for m in mats]
        vrnt_genpos_ls = [m.vrnt_genpos for m in mats]
        vrnt_xoprob_ls = [m.vrnt_xoprob for m in mats]
        vrnt_hapgrp_ls = [m.vrnt_hapgrp for m in mats]
        vrnt_hapalt_ls = [m.vrnt_hapalt for m in mats]
        vrnt_hapref_ls = [m.vrnt_hapref for m in mats]
        vrnt_mask_ls = [m.vrnt_mask for m in mats]

        # process/error check vrnt_chrgrp_ls
        if all(e is None for e in vrnt_chrgrp_ls):
            vrnt_chrgrp_ls = None
        elif any(e is None for e in vrnt_chrgrp_ls):
            raise ValueError("cannot concat: 'vrnt_chrgrp' needed for all Matrix in list")

        # process/error check vrnt_phypos_ls
        if all(e is None for e in vrnt_phypos_ls):
            vrnt_phypos_ls = None
        elif any(e is None for e in vrnt_phypos_ls):
            raise ValueError("cannot concat: 'vrnt_phypos' needed for all Matrix in list")

        # process/error check vrnt_name_ls
        if all(e is None for e in vrnt_name_ls):
            vrnt_name_ls = None
        else:                                                               # else at least one matrix does not have a vrnt_name array
            for i,v in enumerate(vrnt_name_ls):                             # for each index,element in vrnt_name_ls
                if v is None:                                               # if element is None
                    nvrnt = shape_t[mats[0].vrnt_axis][i]                   # get number of variants
                    vrnt_name_ls[i] = numpy.empty(nvrnt, dtype = "object")  # replace with array of None

        # process/error check vrnt_genpos_ls
        if all(e is None for e in vrnt_genpos_ls):
            vrnt_genpos_ls = None
        elif any(e is None for e in vrnt_genpos_ls):
            raise ValueError("cannot concat: 'vrnt_genpos' needed for all Matrix in list")

        # process/error check vrnt_xoprob_ls
        if all(e is None for e in vrnt_xoprob_ls):
            vrnt_xoprob_ls = None
        elif any(e is None for e in vrnt_xoprob_ls):
            raise ValueError("cannot concat: 'vrnt_xoprob' needed for all Matrix in list")

        # process/error check vrnt_hapgrp_ls
        if all(e is None for e in vrnt_hapgrp_ls):
            vrnt_hapgrp_ls = None
        elif any(e is None for e in vrnt_hapgrp_ls):
            raise ValueError("cannot concat: 'vrnt_hapgrp' needed for all Matrix in list")

        # process/error check vrnt_hapalt_ls
        if all(e is None for e in vrnt_hapalt_ls):
            vrnt_hapalt_ls = None
        elif any(e is None for e in vrnt_hapalt_ls):
            raise ValueError("cannot concat: 'vrnt_hapalt' needed for all Matrix in list")

        # process/error check vrnt_hapref_ls
        if all(e is None for e in vrnt_hapref_ls):
            vrnt_hapref_ls = None
        elif any(e is None for e in vrnt_hapref_ls):
            raise ValueError("cannot concat: 'vrnt_hapref' needed for all Matrix in list")

        # process/error check vrnt_mask_ls
        if all(e is None for e in vrnt_mask_ls):
            vrnt_mask_ls = None
        elif any(e is None for e in vrnt_mask_ls):
            raise ValueError("cannot concat: 'vrnt_mask' needed for all Matrix in list")

        # concatenate mat, variant items
        mat = numpy.concatenate(mat_ls, axis = mats[0].vrnt_axis)
        vrnt_chrgrp = None if vrnt_chrgrp_ls is None else numpy.concatenate(vrnt_chrgrp_ls, axis = 0)
        vrnt_phypos = None if vrnt_phypos_ls is None else numpy.concatenate(vrnt_phypos_ls, axis = 0)
        vrnt_name = None if vrnt_name_ls is None else numpy.concatenate(vrnt_name_ls, axis = 0)
        vrnt_genpos = None if vrnt_genpos_ls is None else numpy.concatenate(vrnt_genpos_ls, axis = 0)
        vrnt_xoprob = None if vrnt_xoprob_ls is None else numpy.concatenate(vrnt_xoprob_ls, axis = 0)
        vrnt_hapgrp = None if vrnt_hapgrp_ls is None else numpy.concatenate(vrnt_hapgrp_ls, axis = 0)
        vrnt_hapalt = None if vrnt_hapalt_ls is None else numpy.concatenate(vrnt_hapalt_ls, axis = 0)
        vrnt_hapref = None if vrnt_hapref_ls is None else numpy.concatenate(vrnt_hapref_ls, axis = 0)
        vrnt_mask = None if vrnt_mask_ls is None else numpy.concatenate(vrnt_mask_ls, axis = 0)

        # TODO: decide if first element in list is good source of information
        # concatenate everything and put into new DenseVariantMatrix
        # use first element as source of variant data
        out = cls(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        return out

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis = -1, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DenseVariantMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.vrnt_axis:
            self.append_vrnt(
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

    def append_vrnt(self, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Append values to the Matrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to append to the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to append to the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to append to the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to append to the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to append to the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to append to the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to append to the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if isinstance(values, self.__class__):
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_hapalt is None:
                vrnt_hapalt = values.vrnt_hapalt
            if vrnt_hapref is None:
                vrnt_hapref = values.vrnt_hapref
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot append: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.vrnt_axis) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))
        if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
            raise TypeError("cannot append: 'vrnt_chrgrp' argument is required")
        if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
            raise TypeError("cannot append: 'vrnt_phypos' argument is required")
        if (self._vrnt_name is not None) and (vrnt_name is None):
            vrnt_name = numpy.empty(values.shape[self.vrnt_axis], dtype = "object")  # fill with None
        if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
            raise TypeError("cannot append: 'vrnt_genpos' argument is required")
        if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
            raise TypeError("cannot append: 'vrnt_xoprob' argument is required")
        if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
            raise TypeError("cannot append: 'vrnt_hapgrp' argument is required")
        if (self._vrnt_hapalt is not None) and (vrnt_hapalt is None):
            raise TypeError("cannot append: 'vrnt_hapalt' argument is required")
        if (self._vrnt_hapref is not None) and (vrnt_hapref is None):
            raise TypeError("cannot append: 'vrnt_hapref' argument is required")
        if (self._vrnt_mask is not None) and (vrnt_mask is None):
            raise TypeError("cannot append: 'vrnt_mask' argument is required")

        # append values
        self._mat = numpy.append(self._mat, values, axis = self.vrnt_axis)

        # set fields
        if self._vrnt_chrgrp is not None:
            self._vrnt_chrgrp = numpy.append(self._vrnt_chrgrp, vrnt_chrgrp, axis = 0)
        if self._vrnt_phypos is not None:
            self._vrnt_phypos = numpy.append(self._vrnt_phypos, vrnt_phypos, axis = 0)
        if self._vrnt_name is not None:
            self._vrnt_name = numpy.append(self._vrnt_name, vrnt_name, axis = 0)
        if self._vrnt_genpos is not None:
            self._vrnt_genpos = numpy.append(self._vrnt_genpos, vrnt_genpos, axis = 0)
        if self._vrnt_xoprob is not None:
            self._vrnt_xoprob = numpy.append(self._vrnt_xoprob, vrnt_xoprob, axis = 0)
        if self._vrnt_hapgrp is not None:
            self._vrnt_hapgrp = numpy.append(self._vrnt_hapgrp, vrnt_hapgrp, axis = 0)
        if self._vrnt_hapalt is not None:
            self._vrnt_hapalt = numpy.append(self._vrnt_hapalt, vrnt_hapalt, axis = 0)
        if self._vrnt_hapref is not None:
            self._vrnt_hapref = numpy.append(self._vrnt_hapref, vrnt_hapref, axis = 0)
        if self._vrnt_mask is not None:
            self._vrnt_mask = numpy.append(self._vrnt_mask, vrnt_mask, axis = 0)

        # reset metadata
        self._vrnt_chrgrp_len = None
        self._vrnt_chrgrp_name = None
        self._vrnt_chrgrp_stix = None
        self._vrnt_chrgrp_spix = None

    def remove(self, obj, axis = -1, **kwargs):
        """
        Remove sub-arrays along an axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to remove the subarray defined by obj.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.vrnt_axis:
            self.remove_vrnt(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def remove_vrnt(self, obj, **kwargs):
        """
        Remove sub-arrays along the variant axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # delete values
        self._mat = numpy.delete(self._mat, obj, axis = self.vrnt_axis)

        if self._vrnt_chrgrp is not None:
            self._vrnt_chrgrp = numpy.delete(self._vrnt_chrgrp, obj, axis = 0)
        if self._vrnt_phypos is not None:
            self._vrnt_phypos = numpy.delete(self._vrnt_phypos, obj, axis = 0)
        if self._vrnt_name is not None:
            self._vrnt_name = numpy.delete(self._vrnt_name, obj, axis = 0)
        if self._vrnt_genpos is not None:
            self._vrnt_genpos = numpy.delete(self._vrnt_genpos, obj, axis = 0)
        if self._vrnt_xoprob is not None:
            self._vrnt_xoprob = numpy.delete(self._vrnt_xoprob, obj, axis = 0)
        if self._vrnt_hapgrp is not None:
            self._vrnt_hapgrp = numpy.delete(self._vrnt_hapgrp, obj, axis = 0)
        if self._vrnt_hapalt is not None:
            self._vrnt_hapalt = numpy.delete(self._vrnt_hapalt, obj, axis = 0)
        if self._vrnt_hapref is not None:
            self._vrnt_hapref = numpy.delete(self._vrnt_hapref, obj, axis = 0)
        if self._vrnt_mask is not None:
            self._vrnt_mask = numpy.delete(self._vrnt_mask, obj, axis = 0)

        # reset metadata
        self._vrnt_chrgrp_name = None
        self._vrnt_chrgrp_stix = None
        self._vrnt_chrgrp_spix = None
        self._vrnt_chrgrp_len = None

    def incorp(self, obj, values, axis = -1, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
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
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.vrnt_axis:
            self.incorp(
                obj = obj,
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    def incorp_vrnt(self, obj, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Incorporate values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to incorporate into the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to incorporate into the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to incorporate into the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to incorporate into the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to incorporate into the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to incorporate into the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to incorporate into the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        # if given a DenseVariantMatrix extract *.mat values
        if isinstance(values, self.__class__):
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_hapalt is None:
                vrnt_hapalt = values.vrnt_hapalt
            if vrnt_hapref is None:
                vrnt_hapref = values.vrnt_hapref
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.vrnt_axis) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))
        if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
            raise TypeError("cannot incorp: 'vrnt_chrgrp' argument is required")
        if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
            raise TypeError("cannot incorp: 'vrnt_phypos' argument is required")
        if (self._vrnt_name is not None) and (vrnt_name is None):
            vrnt_name = numpy.empty(values.shape[self.vrnt_axis], dtype = "object") # fill with None
        if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
            raise TypeError("cannot incorp: 'vrnt_genpos' argument is required")
        if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
            raise TypeError("cannot incorp: 'vrnt_xoprob' argument is required")
        if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
            raise TypeError("cannot incorp: 'vrnt_hapgrp' argument is required")
        if (self._vrnt_hapalt is not None) and (vrnt_hapalt is None):
            raise TypeError("cannot incorp: 'vrnt_hapalt' argument is required")
        if (self._vrnt_hapref is not None) and (vrnt_hapref is None):
            raise TypeError("cannot incorp: 'vrnt_hapref' argument is required")
        if (self._vrnt_mask is not None) and (vrnt_mask is None):
            raise TypeError("cannot incorp: 'vrnt_mask' argument is required")

        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = self.vrnt_axis)
        if self._vrnt_chrgrp is not None:
            self._vrnt_chrgrp = numpy.insert(self._vrnt_chrgrp, obj, vrnt_chrgrp, axis = 0)
        if self._vrnt_phypos is not None:
            self._vrnt_phypos = numpy.insert(self._vrnt_phypos, obj, vrnt_phypos, axis = 0)
        if self._vrnt_name is not None:
            self._vrnt_name = numpy.insert(self._vrnt_name, obj, vrnt_name, axis = 0)
        if self._vrnt_genpos is not None:
            self._vrnt_genpos = numpy.insert(self._vrnt_genpos, obj, vrnt_genpos, axis = 0)
        if self._vrnt_xoprob is not None:
            self._vrnt_xoprob = numpy.insert(self._vrnt_xoprob, obj, vrnt_xoprob, axis = 0)
        if self._vrnt_hapgrp is not None:
            self._vrnt_hapgrp = numpy.insert(self._vrnt_hapgrp, obj, vrnt_hapgrp, axis = 0)
        if self._vrnt_hapalt is not None:
            self._vrnt_hapalt = numpy.insert(self._vrnt_hapalt, obj, vrnt_hapalt, axis = 0)
        if self._vrnt_hapref is not None:
            self._vrnt_hapref = numpy.insert(self._vrnt_hapref, obj, vrnt_hapref, axis = 0)
        if self._vrnt_mask is not None:
            self._vrnt_mask = numpy.insert(self._vrnt_mask, obj, vrnt_mask, axis = 0)

        # reset metadata
        self._vrnt_chrgrp_name = None
        self._vrnt_chrgrp_stix = None
        self._vrnt_chrgrp_spix = None
        self._vrnt_chrgrp_len = None

    ################### Sorting Methods ####################
    def lexsort(self, keys = None, axis = -1, **kwargs):
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
        axis = get_axis(axis, self.mat_ndim)    # transform axis number to an index
        indices = None                          # declare variable

        # dispatch to correct function
        if axis == self.vrnt_axis:
            self.lexsort_vrnt(keys = keys, **kwargs)
        else:
            raise ValueError("cannot lexsort along axis {0}".format(axis))

        return indices

    def lexsort_vrnt(self, keys = None, **kwargs):
        """
        Perform an indirect stable sort using a sequence of keys along the
        variant axis.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
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
            keys = (self._vrnt_phypos, self._vrnt_chrgrp)   # loci default keys
            emess = "vrnt_phypos, vrnt_chrgrp are None"     # loci error message

        # remove None keys
        keys = tuple(k for k in keys if k is not None)

        # raise error if no keys remain
        if len(keys) == 0:
            raise ValueError("cannot lexsort on axis {0} (variant axis): {1}".format(self.vrnt_axis, emess))

        # raise error if keys are of incompatible length
        if any(len(k) != self.nvrnt for k in keys):
            emess = "keys are not all length {0}".format(self.nvrnt)
            raise ValueError("cannot lexsort on axis {0} (variant axis): {1}".format(self.vrnt_axis, emess))

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

        if axis == self.vrnt_axis:
            self.reorder(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    def reorder_vrnt(self, indices, **kwargs):
        """
        Reorder elements of the Matrix along the variant axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # build a tuple to slice the matrix
        ix = tuple(indices if i == self.vrnt_axis else slice(None) for i in range(self.mat_ndim))

        # reorder arrays
        self._mat = self._mat[ix]

        if self._vrnt_chrgrp is not None:
            self._vrnt_chrgrp = self._vrnt_chrgrp[indices]  # reorder chromosome group array
        if self._vrnt_phypos is not None:
            self._vrnt_phypos = self._vrnt_phypos[indices]  # reorder physical position array
        if self._vrnt_name is not None:
            self._vrnt_name = self._vrnt_name[indices]      # reorder marker name array
        if self._vrnt_genpos is not None:
            self._vrnt_genpos = self._vrnt_genpos[indices]  # reorder map position array
        if self._vrnt_xoprob is not None:
            self._vrnt_xoprob = self._vrnt_xoprob[indices]  # reorder crossover probability array
        if self._vrnt_hapgrp is not None:
            self._vrnt_hapgrp = self._vrnt_hapgrp[indices]  # reorder haplotype group array
        if self._vrnt_hapalt is not None:
            self._vrnt_hapalt = self._vrnt_hapalt[indices]  # reorder haplotype group array
        if self._vrnt_hapref is not None:
            self._vrnt_hapref = self._vrnt_hapref[indices]  # reorder haplotype group array
        if self._vrnt_mask is not None:
            self._vrnt_mask = self._vrnt_mask[indices]      # reorder variant mask array

    def sort(self, keys = None, axis = -1, **kwargs):
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
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.vrnt_axis:
            self.sort_vrnt(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    def sort_vrnt(self, keys = None, **kwargs):
        """
        Sort slements of the Matrix along the variant axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        # reset variant group metadata
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len = None

        # get indices for sort
        indices = self.lexsort_vrnt(keys, **kwargs)

        # reorder internals
        self.reorder_vrnt(indices, **kwargs)

    ################### Grouping Methods ###################
    def group(self, axis = -1, **kwargs):
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.vrnt_axis:
            self.group_vrnt(**kwargs)
        else:
            raise ValueError("cannot group along axis {0}".format(axis))

    def group_vrnt(self, **kwargs):
        """
        Sort the Matrix along the variant axis, then populate grouping indices
        for the variant axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # sort variants using default keys
        self.sort_vrnt()

        if self._vrnt_chrgrp is not None:
            # get unique chromosome group names, starting indices, group lengths
            uniq = numpy.unique(self._vrnt_chrgrp, return_index = True, return_counts = True)
            # make assignments to instance data
            self._vrnt_chrgrp_name, self._vrnt_chrgrp_stix, self._vrnt_chrgrp_len = uniq
            # calculate stop indices
            self._vrnt_chrgrp_spix = self._vrnt_chrgrp_stix + self._vrnt_chrgrp_len

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

        if axis == self.vrnt_axis:
            grouped = self.is_grouped_vrnt(**kwargs)
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped

    def is_grouped_vrnt(self, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped along the
        variant axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        return (
            (self._vrnt_chrgrp_name is not None) and
            (self._vrnt_chrgrp_stix is not None) and
            (self._vrnt_chrgrp_spix is not None) and
            (self._vrnt_chrgrp_len is not None)
        )



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseVariantMatrix(v):
    """
    Determine whether an object is a DenseVariantMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseVariantMatrix object instance.
    """
    return isinstance(v, DenseVariantMatrix)

def check_is_DenseVariantMatrix(v, varname):
    """
    Check if object is of type DenseVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_DenseVariantMatrix(v):
        raise TypeError("'{0}' must be a DenseVariantMatrix".format(varname))

def cond_check_is_DenseVariantMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DenseVariantMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        DenseVariantMatrix.
    """
    if cond(v):
        check_is_DenseVariantMatrix(v, varname)
