from . import SortableMatrix

class GroupableMatrix(SortableMatrix):
    """docstring for GroupableMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GroupableMatrix, self).__init__()

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
def is_GroupableMatrix(v):
    return isinstance(v, GroupableMatrix)

def check_is_GroupableMatrix(v, varname):
    if not isinstance(v, GroupableMatrix):
        raise TypeError("'%s' must be a GroupableMatrix." % varname)

def cond_check_is_GroupableMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GroupableMatrix(v, varname)
