"""
Module defining comma separated value I/O interfaces and associated error
checking routines.
"""

from typing import Any


class CSVInputOutput:
    """
    Abstract class for defining CSV input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_csv`` - write an object to a csv file.
    - ``from_csv`` - load an object from a csv file.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class CSVInputOutput.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(CSVInputOutput, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ####################### File I/O #######################
    def to_csv(self, filename):
        """
        Write and object to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ####################### File I/O #######################
    @classmethod
    def from_csv(cls, filename):
        """
        Read object from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.

        Returns
        -------
        out : cls
            An object of the appropriate class type.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_CSVInputOutput(v: Any) -> bool:
    """
    Determine whether an object is a CSVInputOutput.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a CSVInputOutput object instance.
    """
    return isinstance(v, CSVInputOutput)

def check_is_CSVInputOutput(v: Any, vname: str) -> None:
    """
    Check if object is of type CSVInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CSVInputOutput):
        raise TypeError("variable '{0}' must be a CSVInputOutput".format(vname))
