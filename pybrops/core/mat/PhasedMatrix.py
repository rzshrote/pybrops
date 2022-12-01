"""
Module defining phased matrix interfaces and associated error checking routines.
"""

from typing import Any

from pybrops.core.mat.MutableMatrix import MutableMatrix

class PhasedMatrix(MutableMatrix):
    """
    An abstract class for phased matrix wrapper objects.

    A phased matrix is defined as a matrix with a third dimension. This
    interface mostly pertains to phased genotype matrices.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix phase manipulation routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        PhasedMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(PhasedMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Phase Metadata Properites ###############
    def nphase():
        doc = "Number of chromosome phases represented by the matrix."
        def fget(self):
            """Get number of phases"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of phases"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of phases"""
            raise NotImplementedError("method is abstract")
        return locals()
    nphase = property(**nphase())

    def phase_axis():
        doc = "Axis along which phases are stored property."
        def fget(self):
            """Get phase axis number"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set phase axis number"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete phase axis number"""
            raise NotImplementedError("method is abstract")
        return locals()
    phase_axis = property(**phase_axis())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin_phase(self, values, **kwargs):
        """
        Add additional elements to the end of the Matrix along the phase axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def delete_phase(self, obj, **kwargs):
        """
        Delete sub-arrays along the phase axis.

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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def select_phase(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the phase axis.

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
    def concat_phase(cls, mats, **kwargs):
        """
        Concatenate list of Matrix together along the phase axis.

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
    def append_phase(self, values, **kwargs):
        """
        Append values to the Matrix along the phase axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def remove_phase(self, obj, **kwargs):
        """
        Remove sub-arrays along the phase axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def incorp_phase(self, obj, values, **kwargs):
        """
        Incorporate values along the phase axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhasedMatrix(v: Any):
    """
    Determine whether an object is a PhasedMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PhasedMatrix object instance.
    """
    return isinstance(v, PhasedMatrix)

def check_is_PhasedMatrix(v: Any, vname: str):
    """
    Check if object is of type PhasedMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_PhasedMatrix(v):
        raise TypeError("'{0}' must be a PhasedMatrix".format(vname))
