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
    def group(self, axis, **kwargs):
        """
        Sort the Matrix, then populate grouping indices.

        Parameters
        ----------
        axis : int
            The axis along which values are grouped.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def is_grouped(self, axis, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped.

        Parameters
        ----------
        axis: int
            Axis along which to determine whether elements have been sorted and
            grouped.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
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
