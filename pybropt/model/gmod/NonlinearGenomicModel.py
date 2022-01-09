from . import GenomicModel

class NonlinearGenomicModel(GenomicModel):
    """docstring for NonlinearGenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for NonlinearGenomicModel class.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(NonlinearGenomicModel, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_NonlinearGenomicModel(v):
    """
    Determine whether an object is a NonlinearGenomicModel.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a NonlinearGenomicModel object instance.
    """
    return isinstance(v, NonlinearGenomicModel)

def check_is_NonlinearGenomicModel(v, vname):
    """
    Check if object is of type NonlinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, NonlinearGenomicModel):
        raise TypeError("variable '{0}' must be a NonlinearGenomicModel".format(vname))

def cond_check_is_NonlinearGenomicModel(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type NonlinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a NonlinearGenomicModel.
    """
    if cond(v):
        check_is_NonlinearGenomicModel(v, vname)
