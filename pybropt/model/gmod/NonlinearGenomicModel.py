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
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(NonlinearGenomicModel, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_NonlinearGenomicModel(v):
    return isinstance(v, NonlinearGenomicModel)

def check_is_NonlinearGenomicModel(v, vname):
    if not isinstance(v, NonlinearGenomicModel):
        raise TypeError("variable '{0}' must be a NonlinearGenomicModel".format(vname))

def cond_check_is_NonlinearGenomicModel(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_NonlinearGenomicModel(v, vname)
