"""
Module defining HDF5 I/O interfaces and associated error checking routines.
"""

class HDF5InputOutput:
    """
    Abstract class for defining HDF5 input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_hdf5`` - write an object to an HDF5 file.
    - ``from_hdf5`` - load an object from an HDF5 file.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class HDF5InputOutput.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(HDF5InputOutput, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Matrix File I/O ####################
    def to_hdf5(self, filename, groupname):
        """
        Write object to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which object data is stored.
            If None, object is written to the base HDF5 group.
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(cls, filename, groupname):
        """
        Read object from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which object data is stored.
            If None, object is read from base HDF5 group.

        Returns
        -------
        obj : cls
            A genotype matrix read from file.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_HDF5InputOutput(v):
    """
    Determine whether an object is a HDF5InputOutput.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a HDF5InputOutput object instance.
    """
    return isinstance(v, HDF5InputOutput)

def check_is_HDF5InputOutput(v, vname):
    """
    Check if object is of type HDF5InputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, HDF5InputOutput):
        raise TypeError("variable '{0}' must be a HDF5InputOutput".format(vname))

def cond_check_is_HDF5InputOutput(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type HDF5InputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a HDF5InputOutput.
    """
    if cond(v):
        check_is_HDF5InputOutput(v, vname)
