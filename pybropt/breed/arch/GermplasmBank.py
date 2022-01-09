from pybropt.breed.arch.BreedingNode import BreedingNode

class GermplasmBank(BreedingNode):
    """docstring for GermplasmBank."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class GermplasmBank.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(GermplasmBank, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_GermplasmBank(v):
    """
    Determine whether an object is a GermplasmBank.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GermplasmBank object instance.
    """
    return isinstance(v, GermplasmBank)

def check_is_GermplasmBank(v, varname):
    """
    Check if object is of type GermplasmBank. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GermplasmBank):
        raise TypeError("'%s' must be a GermplasmBank." % varname)

def cond_check_is_GermplasmBank(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type GermplasmBank. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a GermplasmBank.
    """
    if cond(v):
        check_is_GermplasmBank(v, varname)
