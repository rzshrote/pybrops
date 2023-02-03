# TODO: finish writing this class
from . import DenseHaplotypeMatrix

from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_ndarray_dtype_is_int8
from pybrops.core.error import check_ndarray_is_3d
from pybrops.core.error import error_readonly
from pybrops.core.mat import PhasedMatrix

# TODO: full implementation
class DensePhasedHaplotypeMatrix(DenseHaplotypeMatrix,PhasedMatrix):
    """docstring for DensePhasedHaplotypeMatrix ."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs: dict):
        """
        DensePhasedHaplotypeMatrix constructor

        Parameters
        ----------
        haplo : numpy.ndarray
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DensePhasedHaplotypeMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            ploidy = mat.shape[0],
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Haplotype Data Properites ##############
    def mat():
        doc = "Raw haplotype matrix."
        def fget(self):
            return self._mat
        def fset(self, value):
            check_is_ndarray(value, "mat")              # check ndarray
            check_ndarray_dtype_is_int8(value, "mat")   # check dtype
            check_ndarray_is_3d(value, "mat")           # check shape
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())

    ############## General matrix properties ###############
    def ploidy():
        doc = "Ploidy of the haplotypes represented by the matrix."
        def fget(self):
            return self._mat.shape[self.phase_axis]
        def fset(self, value):
            error_readonly("ploidy")
        def fdel(self):
            error_readonly("ploidy")
        return locals()
    ploidy = property(**ploidy())

    ############## Phase Metadata Properites ###############
    def nphase():
        doc = "Number of chromosome phases represented by the matrix."
        def fget(self):
            """Get number of phases"""
            return self._mat.shape[self.phase_axis]
        def fset(self, value):
            """Set number of phases"""
            error_readonly("nphase")
        def fdel(self):
            """Delete number of phases"""
            error_readonly("nphase")
        return locals()
    nphase = property(**nphase())

    def phase_axis():
        doc = "Axis along which phases are stored property."
        def fget(self):
            """Get phase axis number"""
            return 0
        def fset(self, value):
            """Set phase axis number"""
            error_readonly("phase_axis")
        def fdel(self):
            """Delete phase axis number"""
            error_readonly("phase_axis")
        return locals()
    phase_axis = property(**phase_axis())

    ############### Taxa Metadata Properites ###############
    def taxa_axis():
        doc = "Axis along which taxa are stored property."
        def fget(self):
            """Get taxa axis number"""
            return 1
        def fset(self, value):
            """Set taxa axis number"""
            error_readonly("taxa_axis")
        def fdel(self):
            """Delete taxa axis number"""
            error_readonly("taxa_axis")
        return locals()
    taxa_axis = property(**taxa_axis())

    ############# Variant Metadata Properites ##############
    def vrnt_axis():
        doc = "Axis along which variants are stored property."
        def fget(self):
            """Get variant axis"""
            return 2
        def fset(self, value):
            """Set variant axis"""
            error_readonly("vrnt_axis")
        def fdel(self):
            """Delete variant axis"""
            error_readonly("vrnt_axis")
        return locals()
    vrnt_axis = property(**vrnt_axis())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs: dict):
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
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
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
        if axis == self.phase_axis:
            out = self.adjoin_phase(
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
            out = self.adjoin_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
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

    def adjoin_phase(self, values: Union[Matrix,numpy.ndarray], **kwargs: dict):
        """
        Adjoin values along the phase axis.

        Parameters
        ----------
        values : Matrix or numpy.ndarray
            Values to adjoin along the phase axis.
        **kwargs
            Additional keyword arguments.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DenseHaplotypeMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("cannot adjoin: 'values' must be of type DenseHaplotypeMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.phase_axis) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))

        # adjoin values
        values = numpy.append(self._mat, values, axis = self.phase_axis)

        out = self.__class__(
            mat = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            ploidy = values.shape[self.phase_axis],
            **kwargs
        )

        return out

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
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.delete_phase(obj = obj, **kwargs)
        elif axis == self.taxa_axis:
            out = self.delete_taxa(obj = obj, **kwargs)
        elif axis == self.vrnt_axis:
            out = self.delete_vrnt(obj = obj, **kwargs)
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_phase(self, obj: Union[int,slice,Sequence], **kwargs: dict):
        """
        Delete sub-arrays along the phase axis.

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
        # get values
        mat = self._mat

        # delete values
        mat = numpy.delete(mat, obj, axis = self.phase_axis)

        out = self.__class__(
            mat = mat,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            ploidy = mat.shape[self.phase_axis],
            **kwargs
        )

        return out

    def insert(self, obj: Union[int,slice,Sequence], values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs: dict):
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
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
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
        if axis == self.phase_axis:
            out = self.insert_phase(
                obj = obj,
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
            out = self.insert_taxa(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
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

    def insert_phase(self, obj: Union[int,slice,Sequence], values: Union[Matrix,numpy.ndarray], **kwargs: dict):
        """
        Insert values along the phase axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # if given a DenseHaplotypeMatrix extract *.mat values
        if is_DenseHaplotypeMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseHaplotypeMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot insert: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.phase_axis) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))

        # insert values
        values = numpy.insert(self._mat, obj, values, axis = self.phase_axis)

        # create output
        out = self.__class__(
            mat = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            ploidy = values.shape[self.phase_axis],
            **kwargs
        )

        return out

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
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.select_phase(indices = indices, **kwargs)
        elif axis == self.taxa_axis:
            out = self.select_taxa(indices = indices, **kwargs)
        elif axis == self.vrnt_axis:
            out = self.select_vrnt(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_phase(self, indices, **kwargs: dict):
        """
        Select certain values from the Matrix along the phase axis.

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

        # select values
        mat = numpy.take(mat, indices, axis = self.phase_axis)

        out = self.__class__(
            mat = mat,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            ploidy = mat.shape[self.phase_axis],
            **kwargs
        )

        return out

    @staticmethod
    def concat(mats, axis = -1, **kwargs: dict):
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
        if axis == mats[0].phase_axis:
            out = DenseHaplotypeMatrix.concat_phase(mats, **kwargs)
        elif axis == mats[0].taxa_axis:
            out = DenseHaplotypeMatrix.concat_taxa(mats, **kwargs)
        elif axis == mats[0].vrnt_axis:
            out = DenseHaplotypeMatrix.concat_vrnt(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @staticmethod
    def concat_phase(mats: Sequence, **kwargs: dict):
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
            check_is_DenseHaplotypeMatrix(v, "mats[{0}]".format(i))

        # make sure dimensions are all identical to first element in mats
        if any(m.mat_ndim != mats[0].mat_ndim for m in mats):
            raise ValueError("cannot concat: not all matrices have the same number of dimensions")

        # extract tuple of shapes for testing compatibility
        shape_t = tuple(zip(*[m.mat.shape for m in mats]))

        # test matrix compatibility (same axis length along non-taxa axes)
        for i,v in enumerate(shape_t):                              # for each index,tuple in shape_t
            if (i != mats[0].phase_axis) and any(l != v[0] for l in v): # if not the taxa axis AND axis lengths are different
                raise ValueError("cannot concat: matrix shapes do not all align along axis {0}".format(i))

        # create matrix lists
        mat_ls = [m.mat for m in mats]

        # concatenate mat, taxa, taxa_grp items
        mat = numpy.concatenate(mat_ls, axis = mats[0].taxa_axis)

        # TODO: decide if first element in list is good source of information
        # concatenate everything and put into new DenseHaplotypeMatrix
        # use first element as source of variant data
        out = DenseHaplotypeMatrix(
            mat = mat,
            taxa = mats[0].taxa,
            taxa_grp = mats[0].taxa_grp,
            vrnt_chrgrp = mats[0].vrnt_chrgrp,
            vrnt_phypos = mats[0].vrnt_phypos,
            vrnt_name = mats[0].vrnt_name,
            vrnt_genpos = mats[0].vrnt_genpos,
            vrnt_xoprob = mats[0].vrnt_xoprob,
            vrnt_hapgrp = mats[0].vrnt_hapgrp,
            vrnt_hapalt = mats[0].vrnt_hapalt,
            vrnt_hapref = mats[0].vrnt_hapref,
            vrnt_mask = mats[0].vrnt_mask,
            ploidy = mat.shape[self.phase_axis],
            **kwargs
        )

        # copy taxa metadata from source
        out.taxa_grp_name = mats[0].taxa_grp_name
        out.taxa_grp_stix = mats[0].taxa_grp_stix
        out.taxa_grp_spix = mats[0].taxa_grp_spix
        out.taxa_grp_len = mats[0].taxa_grp_len

        # copy variant metadata from source
        out.vrnt_chrgrp_name = mats[0].vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = mats[0].vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = mats[0].vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = mats[0].vrnt_chrgrp_len

        return out

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs: dict):
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
        if axis == self.phase_axis:
            self.append_phase(values, **kwargs)
        elif axis == self.taxa_axis:
            self.append_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
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

    def append_phase(self, values: Union[Matrix,numpy.ndarray], **kwargs: dict):
        """
        Append values to the Matrix along the phase axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        **kwargs
            Additional keyword arguments.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DenseHaplotypeMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseHaplotypeMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot append: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))

        # append values
        self._mat = numpy.append(self._mat, values, axis = self.phase_axis)

    def remove(self, obj: Union[int,slice,Sequence], axis = -1, **kwargs: dict):
        """
        Remove sub-arrays along an axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to remove the subarray defined by obj.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.phase_axis:
            self.remove_phase(obj = obj, **kwargs)
        elif axis == self.taxa_axis:
            self.remove_taxa(obj = obj, **kwargs)
        elif axis == self.vrnt_axis:
            self.remove_vrnt(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def remove_phase(self, obj: Union[int,slice,Sequence], **kwargs: dict):
        """
        Remove sub-arrays along the phase axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        # delete values
        self._mat = numpy.delete(self._mat, obj, axis = self.phase_axis)

    def incorp(self, obj, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs: dict):
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

        if axis == self.phase_axis:
            self.incorp(
                obj = obj,
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
            self.incorp(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
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

    def incorp_phase(self, obj, values: Union[Matrix,numpy.ndarray], **kwargs: dict):
        """
        Incorporate values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        **kwargs
            Additional keyword arguments.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DenseHaplotypeMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseHaplotypeMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.phase_axis) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))

        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = self.phase_axis)

    ############## Matrix summary statistics ###############
    def thcount(self, dtype = None):
        """
        Haplotype count of the non-zero haplotype within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, h) containing haplotype counts of the
            haplotype coded as 1 for all 'n' individuals, for all 'h'
            haplotypes.
        """
        out = self._mat.sum(0)              # take sum across the phase axis (0)
        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype
        return out

    def thfreq(self, dtype = None):
        """
        Haplotype frequency of the non-zero haplotype within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, h) containing haplotype frequencies of
            the haplotype coded as 1 for all 'n' individuals, for all 'h'
            haplotypes.
        """
        rnphase = 1.0 / self._mat.shape[0]  # take the reciprocal of the number of phases = 1 / nphase
        out = rnphase * self._mat.sum(0)    # take sum across the phase axis (0) and divide by nphase
        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype
        return out

    def hcount(self, dtype = None):
        """
        Haplotype count of the non-zero haplotype across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype counts of the
            haplotype coded as 1 for all 'h' haplotypes.
        """
        out = self._mat.sum((0,1))      # take sum across the phase (0) and taxa (1) axes
        if dtype is not None:           # if dtype is specified
            dtype = numpy.dtype(dtype)  # ensure conversion to dtype class
            if out.dtype != dtype:      # if output dtype and desired are different
                out = dtype.type(out)   # convert to correct dtype
        return out

    def hfreq(self, dtype = None):
        """
        Haplotype frequency of the non-zero haplotype across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype frequencies of
            the haplotype coded as 1 for all 'h' haplotypes.
        """
        s = self._mat.shape                     # get matrix shape
        denom = (s[0] * s[1])                   # get nphase * ntaxa
        rnphase = 1.0 / denom                   # take 1 / (nphase * ntaxa)
        out = rnphase * self._mat.sum((0,1))    # take sum across the phase axis (0) and divide by nphase
        if dtype is not None:                   # if dtype is specified
            dtype = numpy.dtype(dtype)          # ensure conversion to dtype class
            if out.dtype != dtype:              # if output dtype and desired are different
                out = dtype.type(out)           # convert to correct dtype
        return out

    def mhf(self, dtype = None):
        """
        Minor haplotype frequency across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype frequencies for
            the minor haplotype.
        """
        out = self.hfreq(dtype)     # get haplotype frequencies
        mask = out > 0.5            # create mask of haplotype frequencies > 0.5
        out[mask] = 1.0 - out[mask] # take 1 - haplotype frequency
        return out

    def mehe(self):
        """
        Mean expected heterozygosity across all taxa.

        Returns
        -------
        mehe : numpy.float64
            A 64-bit floating point representing the mean expected
            heterozygosity.
        """
        p = self.hfreq()            # get haplotype frequency (p)
        out = (p * (1.0 - p)).sum() # take p*(1-p) across all loci and sum the products
        s = self._mat.shape         # get matrix shape
        out *= (s[0] / s[2])        # multiply summation by (nphase / nloci)
        return numpy.float64(out)

    def gtcount(self):
        """
        Gather haplotype counts across for homozygous major, heterozygous,
        homozygous minor all individuals.

        Returns
        -------
        out : numpy.ndarray
            An int64 array of shape (3, h) containing haplotype counts across
            all 'h' haplotypes.
            Rows are as follows:
                out[0] = count of '0' genotype across all loci
                out[1] = count of '1' genotype across all loci
                out[2] = count of '2' genotype across all loci
        """
        gt = self._mat.sum(0)           # get genotypes as {0, 1, 2}
        out = numpy.empty(              # allocate output array
            (3, self._mat.shape[2]),    # shape (3, h)
            dtype='int64'               # int64 dtype
        )

        # get genotype counts
        out[0] = (gt == 0).sum(0)       # homozygous 0/0
        out[1] = (gt == 1).sum(0)       # heterozygous 0/1, 1/0
        out[2] = (gt == 2).sum(0)       # homozygous 1/1

        return out

    def gtfreq(self):
        """
        Gather haplotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Returns
        -------
        out : numpy.ndarray
            An float64 array of shape (3, h) containing haplotype counts across
            all 'h' haplotypes.
            Rows are as follows:
                out[0] = frequency of '0' genotype across all loci
                out[1] = frequency of '1' genotype across all loci
                out[2] = frequency of '2' genotype across all loci
        """
        recip = 1.0 / self._mat.shape[1]    # get reciprocal of number of taxa
        out = recip * self.gtcount()        # calculate genotype frequencies
        return out

    ################### Matrix File I/O ####################
    @staticmethod
    def from_hdf5(filename, groupname):
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
        raise NotImplementedError("method is abstract")

    def to_hdf5(self, filename, groupname):
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
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhasedHaplotypeMatrix(v):
    return isinstance(v, DensePhasedHaplotypeMatrix)

def check_is_DensePhasedHaplotypeMatrix(v, varname):
    if not isinstance(v, DensePhasedHaplotypeMatrix):
        raise TypeError("'%s' must be a DensePhasedHaplotypeMatrix." % varname)

def cond_check_is_DensePhasedHaplotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DensePhasedHaplotypeMatrix(v, varname)
