from . import BreedingEdge

class EmigrationOperator(BreedingEdge):
    """docstring for EmigrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(EmigrationOperator, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def emigrate(self, bnode, **kwargs):
        """
        Emigrate individuals to a BreedingNode.

        Parameters
        ----------
        bnode : BreedingNode
            A BreedingNode object to add individuals.
        kwargs : dict
            Dictionary of data for individuals to be added.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_EmigrationOperator(v):
    return isinstance(v, EmigrationOperator)

def check_is_EmigrationOperator(v, varname):
    if not isinstance(v, EmigrationOperator):
        raise TypeError("'%s' must be a EmigrationOperator." % varname)

def cond_check_is_EmigrationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_EmigrationOperator(v, varname)
