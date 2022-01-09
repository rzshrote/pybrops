from . import GroupableMatrix

class VariantMatrix(GroupableMatrix):
    """
    An abstract class for matrix wrapper objects with variant metadata.

    The purpose of this abstract class is to provide base functionality for:
        1) variant metadata manipulation routines.
        2) variant manipulation routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        VariantMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(VariantMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Variant Data Properites ################
    def vrnt_chrgrp():
        doc = "Variant chromosome group label property."
        def fget(self):
            """Get variant chromosome group lable array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group lable array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group lable array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp = property(**vrnt_chrgrp())

    def vrnt_phypos():
        doc = "Variant physical position property."
        def fget(self):
            """Get variant physical position array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant physical position array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant physical position array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_phypos = property(**vrnt_phypos())

    def vrnt_name():
        doc = "Variant name property."
        def fget(self):
            """Get variant name array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant name array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant name array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_name = property(**vrnt_name())

    def vrnt_genpos():
        doc = "Variant genetic position property."
        def fget(self):
            """Get variant genetic position array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant genetic position array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant genetic position array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    def vrnt_xoprob():
        doc = "Variant crossover sequential probability property."
        def fget(self):
            """Get variant crossover sequential probability array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant crossover sequential probability array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant crossover sequential probability array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_xoprob = property(**vrnt_xoprob())

    def vrnt_hapgrp():
        doc = "Variant haplotype group label property."
        def fget(self):
            """Get variant haplotype group label array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant haplotype group label array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant haplotype group label array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_hapgrp = property(**vrnt_hapgrp())

    def vrnt_hapalt():
        doc = "Variant haplotype sequence property."
        def fget(self):
            """Get variant haplotype sequence"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant haplotype sequence"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant haplotype sequence"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_hapalt = property(**vrnt_hapalt())

    def vrnt_hapref():
        doc = "Variant reference haplotype sequence property."
        def fget(self):
            """Get variant reference haplotype sequence"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant reference haplotype sequence"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant reference haplotype sequence"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_hapref = property(**vrnt_hapref())

    def vrnt_mask():
        doc = "Variant mask property."
        def fget(self):
            """Get variant mask"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant mask"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant mask"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_mask = property(**vrnt_mask())

    ############# Variant Metadata Properites ##############
    def nvrnt():
        doc = "Number of variants property."
        def fget(self):
            """Get number of variants"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of variants"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of variants"""
            raise NotImplementedError("method is abstract")
        return locals()
    nvrnt = property(**nvrnt())

    def vrnt_axis():
        doc = "Axis along which variants are stored property."
        def fget(self):
            """Get variant axis"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant axis"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant axis"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_axis = property(**vrnt_axis())

    def vrnt_chrgrp_name():
        doc = "Variant chromosome group names property."
        def fget(self):
            """Get variant chromosome group name array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group name array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group name array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_name = property(**vrnt_chrgrp_name())

    def vrnt_chrgrp_stix():
        doc = "Variant chromosome group start indices property."
        def fget(self):
            """Get variant chromosome group start indices array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group start indices array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group start indices array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_stix = property(**vrnt_chrgrp_stix())

    def vrnt_chrgrp_spix():
        doc = "Variant chromosome group stop indices property."
        def fget(self):
            """Get variant chromosome group stop indices array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group stop indices array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group stop indices array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_spix = property(**vrnt_chrgrp_spix())

    def vrnt_chrgrp_len():
        doc = "Variant chromosome group length property."
        def fget(self):
            """Get variant chromosome group length array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group length array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group length array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_len = property(**vrnt_chrgrp_len())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin_vrnt(self, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs):
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
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def delete_vrnt(self, obj, **kwargs):
        """
        Delete sub-arrays along the variant axis.

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
        raise NotImplementedError("static method is abstract")

    def insert_vrnt(self, obj, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs):
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
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def select_vrnt(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the variant axis.

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
        raise NotImplementedError("method is abstract")

    @classmethod
    def concat_vrnt(cls, mats, **kwargs):
        """
        Concatenate list of Matrix together along the variant axis.

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
        raise NotImplementedError("static method is abstract")

    ######### Matrix element in-place-manipulation #########
    def append_vrnt(self, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs):
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
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def remove_vrnt(self, obj, **kwargs):
        """
        Remove sub-arrays along the variant axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def incorp_vrnt(self, obj, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs):
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
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Sorting Methods ####################
    def lexsort_vrnt(self, keys, **kwargs):
        """
        Perform an indirect stable sort using a sequence of keys along the
        variant axis.

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
        raise NotImplementedError("method is abstract")

    def reorder_vrnt(self, indices, **kwargs):
        """
        Reorder elements of the Matrix along the variant axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def sort_vrnt(self, keys, **kwargs):
        """
        Sort slements of the Matrix along the variant axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Grouping Methods ###################
    def group_vrnt(self, **kwargs):
        """
        Sort the Matrix along the variant axis, then populate grouping indices
        for the variant axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def is_grouped_vrnt(self, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped along the
        variant axis.

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
        raise NotImplementedError("method is abstract")

    ################# Interpolation Methods ################
    def interp_genpos(self, gmap, **kwargs):
        """
        Interpolate genetic map postions for variants using a GeneticMap

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def interp_xoprob(self, gmap, gmapfn, **kwargs):
        """
        Interpolate genetic map positions AND crossover probabilities between
        sequential markers using a GeneticMap and a GeneticMapFunction.

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        gmapfn : GeneticMapFunction
            A genetic map function from which to interpolate crossover
            probabilities for loci within the VariantMatrix.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################## Clustering Methods ##################
    def assign_hapgrp(self, k, **kwargs):
        """
        Assign haplotype groups using k-means clustering.

        Parameters
        ----------
        k : int, numpy.ndarray
            Number of haplotype groups to assign to each
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_VariantMatrix(v):
    """
    Determine whether an object is a VariantMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a VariantMatrix object instance.
    """
    return isinstance(v, VariantMatrix)

def check_is_VariantMatrix(v, varname):
    """
    Check if object is of type VariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_VariantMatrix(v):
        raise TypeError("'{0}' must be a VariantMatrix".format(varname))

def cond_check_is_VariantMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type VariantMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        VariantMatrix.
    """
    if cond(v):
        check_is_VariantMatrix(v, varname)
