"""
Module containing the abstract class GenotypingProtocol and its service functions.
"""

class GenotypingProtocol:
    """docstring for GenotypingProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class GenotypingProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(GenotypingProtocol, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def genotype(self, pgmat, miscout, **kwargs):
        """
        Genotype a genome.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A phased genotype matrix representing the whole simulated genome.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : GenotypeMatrix
            A GenotypeMatrix of genotyped individuals.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenotypingProtocol(v):
    """
    Determine whether an object is a GenotypingProtocol.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GenotypingProtocol object instance.
    """
    return isinstance(v, GenotypingProtocol)

def check_is_GenotypingProtocol(v, varname):
    """
    Check if object is of type GenotypingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GenotypingProtocol):
        raise TypeError("'%s' must be a GenotypingProtocol." % varname)

def cond_check_is_GenotypingProtocol(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type GenotypingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a GenotypingProtocol.
    """
    if cond(v):
        check_is_GenotypingProtocol(v, varname)
