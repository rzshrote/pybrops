"""
Module defining interfaces and associated error checking routines for matrices
with taxa metadata.
"""

from typing import Any
from pybrops.core.mat.GroupableMatrix import GroupableMatrix

class TaxaMatrix(GroupableMatrix):
    """
    An abstract class for matrix wrapper objects with taxa metadata.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix taxa metadata.
        2) Matrix taxa routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        TaxaMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(TaxaMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Taxa Data Properites #################
    def taxa():
        doc = "Taxa label property."
        def fget(self):
            """Get taxa label array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set taxa label array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete taxa label array"""
            raise NotImplementedError("method is abstract")
        return locals()
    taxa = property(**taxa())

    def taxa_grp():
        doc = "Taxa group label property."
        def fget(self):
            """Get taxa group label array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set taxa group label array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete taxa group label array"""
            raise NotImplementedError("method is abstract")
        return locals()
    taxa_grp = property(**taxa_grp())

    ############### Taxa Metadata Properites ###############
    def ntaxa():
        doc = "Number of taxa property."
        def fget(self):
            """Get number of taxa"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of taxa"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of taxa"""
            raise NotImplementedError("method is abstract")
        return locals()
    ntaxa = property(**ntaxa())

    def taxa_axis():
        doc = "Axis along which taxa are stored property."
        def fget(self):
            """Get taxa axis number"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set taxa axis number"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete taxa axis number"""
            raise NotImplementedError("method is abstract")
        return locals()
    taxa_axis = property(**taxa_axis())

    def taxa_grp_name():
        doc = "Taxa group name property."
        def fget(self):
            """Get taxa group name array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set taxa group name array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete taxa group array"""
            raise NotImplementedError("method is abstract")
        return locals()
    taxa_grp_name = property(**taxa_grp_name())

    def taxa_grp_stix():
        doc = "Taxa group start index property."
        def fget(self):
            """Get taxa group start indices array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set taxa group start indices array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete taxa group start indices array"""
            raise NotImplementedError("method is abstract")
        return locals()
    taxa_grp_stix = property(**taxa_grp_stix())

    def taxa_grp_spix():
        doc = "Taxa group stop index property."
        def fget(self):
            """Get taxa group stop indices array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set taxa group stop indices array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete taxa group stop indices array"""
            raise NotImplementedError("method is abstract")
        return locals()
    taxa_grp_spix = property(**taxa_grp_spix())

    def taxa_grp_len():
        doc = "Taxa group length property."
        def fget(self):
            """Get taxa group length array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set taxa group length array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete taxa group length array"""
            raise NotImplementedError("method is abstract")
        return locals()
    taxa_grp_len = property(**taxa_grp_len())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin_taxa(self, values, taxa, taxa_grp, **kwargs):
        """
        Add additional elements to the end of the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def delete_taxa(self, obj, **kwargs):
        """
        Delete sub-arrays along the taxa axis.

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
        raise NotImplementedError("static method is abstract")

    def insert_taxa(self, obj, values, taxa, taxa_grp, **kwargs):
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def select_taxa(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the taxa axis.

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
        raise NotImplementedError("method is abstract")

    @classmethod
    def concat_taxa(cls, mats, **kwargs):
        """
        Concatenate list of Matrix together along the taxa axis.

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
        raise NotImplementedError("static method is abstract")

    ######### Matrix element in-place-manipulation #########
    def append_taxa(self, values, taxa, taxa_grp, **kwargs):
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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def remove_taxa(self, obj, **kwargs):
        """
        Remove sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def incorp_taxa(self, obj, values, taxa, taxa_grp, **kwargs):
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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Sorting Methods ####################
    def lexsort_taxa(self, keys, **kwargs):
        """
        Perform an indirect stable sort using a sequence of keys along the taxa
        axis.

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
        raise NotImplementedError("method is abstract")

    def reorder_taxa(self, indices, **kwargs):
        """
        Reorder elements of the Matrix along the taxa axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def sort_taxa(self, keys, **kwargs):
        """
        Sort slements of the Matrix along the taxa axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Grouping Methods ###################
    def group_taxa(self, **kwargs):
        """
        Sort the Matrix along the taxa axis, then populate grouping indices for
        the taxa axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def is_grouped_taxa(self, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped along the taxa
        axis.

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
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_TaxaMatrix(v: Any) -> bool:
    """
    Determine whether an object is a TaxaMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a TaxaMatrix object instance.
    """
    return isinstance(v, TaxaMatrix)

def check_is_TaxaMatrix(v: Any, varname: str) -> None:
    """
    Check if object is of type TaxaMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TaxaMatrix):
        raise TypeError("'{0}' must be a TaxaMatrix".format(varname))
