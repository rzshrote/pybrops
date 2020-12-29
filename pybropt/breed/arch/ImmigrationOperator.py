from . import BreedingEdge

class ImmigrationOperator(BreedingEdge):
    """docstring for ImmigrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(ImmigrationOperator, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def immigrate(self, bnode, **kwargs):
        """
        Immigrate individuals from a BreedingNode.

        Parameters
        ----------
        bnode : BreedingNode
            A BreedingNode object from which to pull individuals.

        Returns
        -------
        kwindiv : dict
            Dictionary of data from individuals selected.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_ImmigrationOperator(v):
    return isinstance(v, ImmigrationOperator)

def check_is_ImmigrationOperator(v, varname):
    if not isinstance(v, ImmigrationOperator):
        raise TypeError("'%s' must be a ImmigrationOperator." % varname)

def cond_check_is_ImmigrationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_ImmigrationOperator(v, varname)
