from . import PhasedTaxaVariantMatrix
from . import DensePhasedMatrix
from . import DenseTaxaVariantMatrix

from . import get_axis
from pybropt.core.error import error_readonly

class DensePhasedTaxaVariantMatrix(DenseTaxaVariantMatrix,DensePhasedMatrix,PhasedTaxaVariantMatrix):
    """docstring for DensePhasedTaxaVariantMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Parameters
        ----------
        mat : numpy.ndarray
            An int8 haplotype matrix. Must be {0,1,2} format.
        """
        super(DensePhasedTaxaVariantMatrix, self).__init__(
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
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Phase Metadata Properites ###############
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
    def adjoin(self, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DensePhasedGenotypeMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
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

    def adjoin_phase(self, values, **kwargs):
        """
        Adjoin values along the phase axis.

        Parameters
        ----------
        values : Matrix or numpy.ndarray
            Values to adjoin along the phase axis.
        **kwargs
            Additional keyword arguments.
        """
        out = super(DensePhasedTaxaVariantMatrix, self).adjoin_phase(
            values = values,
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
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

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
        if axis == self.phase_axis:
            out = self.delete_phase(obj = obj, **kwargs)
        elif axis == self.taxa_axis:
            out = self.delete_taxa(obj = obj, **kwargs)
        elif axis == self.vrnt_axis:
            out = self.delete_vrnt(obj = obj, **kwargs)
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_phase(self, obj, **kwargs):
        """
        Delete sub-arrays along the phase axis.

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
        out = super(DensePhasedTaxaVariantMatrix, self).delete_vrnt(
            obj = obj,
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
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    def insert(self, obj, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
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

    def insert_phase(self, obj, values, **kwargs):
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
        out = super(DensePhasedTaxaVariantMatrix, self).insert_phase(
            obj = obj,
            values = values,
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
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

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
        if axis == self.phase_axis:
            out = self.select_phase(indices = indices, **kwargs)
        elif axis == self.taxa_axis:
            out = self.select_taxa(indices = indices, **kwargs)
        elif axis == self.vrnt_axis:
            out = self.select_vrnt(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_phase(self, indices, **kwargs):
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
        out = super(DensePhasedTaxaVariantMatrix, self).select_phase(
            indices = indices,
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
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

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
        if axis == mats[0].phase_axis:
            out = cls.concat_phase(mats, **kwargs)
        elif axis == mats[0].taxa_axis:
            out = cls.concat_taxa(mats, **kwargs)
        elif axis == mats[0].vrnt_axis:
            out = cls.concat_vrnt(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_phase(cls, mats, **kwargs):
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
        out = super(DensePhasedTaxaVariantMatrix, cls).concat_phase(
            mats = mats,
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
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = mats[0].taxa_grp_name
        out.taxa_grp_stix = mats[0].taxa_grp_stix
        out.taxa_grp_spix = mats[0].taxa_grp_spix
        out.taxa_grp_len = mats[0].taxa_grp_len
        out.vrnt_chrgrp_name = mats[0].vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = mats[0].vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = mats[0].vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = mats[0].vrnt_chrgrp_len

        return out

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
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
            self.append_phase(
                values = values,
                **kwargs
            )
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

        if axis == self.phase_axis:
            self.remove_phase(obj = obj, **kwargs)
        elif axis == self.taxa_axis:
            self.remove_taxa(obj = obj, **kwargs)
        elif axis == self.vrnt_axis:
            self.remove_vrnt(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def incorp(self, obj, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
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
            self.incorp_phase(
                obj = obj,
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
            self.incorp_taxa(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
            self.incorp_vrnt(
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

    ################### Sorting Methods ####################

    ################### Grouping Methods ###################
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

        if axis == self.phase_axis:
            grouped = False     # always results to false, even though axis is ungroupable
        elif axis == self.taxa_axis:
            grouped = self.is_grouped_taxa(**kwargs)
        elif axis == self.vrnt_axis:
            grouped = self.is_grouped_vrnt(**kwargs)
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped


################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhasedTaxaVariantMatrix(v):
    """
    Determine whether an object is a DensePhasedTaxaVariantMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DensePhasedTaxaVariantMatrix object instance.
    """
    return isinstance(v, DensePhasedTaxaVariantMatrix)

def check_is_DensePhasedTaxaVariantMatrix(v, varname):
    """
    Check if object is of type DensePhasedTaxaVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_DensePhasedTaxaVariantMatrix(v):
        raise TypeError("'{0}' must be a DensePhasedTaxaVariantMatrix".format(varname))

def cond_check_is_DensePhasedTaxaVariantMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DensePhasedTaxaVariantMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        DensePhasedTaxaVariantMatrix.
    """
    if cond(v):
        check_is_DensePhasedTaxaVariantMatrix(v, varname)
