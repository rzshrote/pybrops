"""
Module containing matrix utilities.
"""

__all__ = [
    "get_axis",
]

def get_axis(
        axis: int, 
        ndim: int
    ) -> int:
    """
    Return an index (unsigned) from a provided axis integer (signed)

    Parameters
    ----------
    axis : int
        Integer representation of the axis. Can be in range (-ndim,ndim).
        If outside this range, will raise an AxisError.
    ndim : int
        Number of dimensions available to index along.

    Returns
    -------
    index : int
        Index representation of the axis. In range [0,ndim).
    """
    # handle axis argument
    if (axis >= ndim) or (axis < -ndim):
        raise IndexError("axis {0} is out of bounds for array of dimension {1}".format(axis, ndim))

    # modulo the axis number to get the axis (in the case of negative axis)
    axis %= ndim

    return axis
