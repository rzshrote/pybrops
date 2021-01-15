class GroupableMatrixInterface:
    """docstring for GroupableMatrixInterface."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GroupableMatrixInterface, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Grouping Methods ###################
    def group(self, axis):
        """
        Sort genetic map, then populate grouping indices.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        If not sorted, raise RuntimeError.
        """
        raise NotImplementedError("method is abstract")

    def is_grouped(self, axis):
        """
        Determine whether the GeneticMap has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GroupableMatrixInterface(v):
    return isinstance(v, GroupableMatrixInterface)

def check_is_GroupableMatrixInterface(v, varname):
    if not isinstance(v, GroupableMatrixInterface):
        raise TypeError("'%s' must be a GroupableMatrixInterface." % varname)

def cond_check_is_GroupableMatrixInterface(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GroupableMatrixInterface(v, varname)
