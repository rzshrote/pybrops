import numpy

__all__ = ["is_ndarray"]

# TODO: consider removal of this checking function. It is rather trivial
def is_ndarray(obj):
    """
    Determine if an object is a numpy.ndarray.

    Parameters
    ----------
    obj : object
        Any Python object.

    Returns
    -------
    out : bool
        ``True`` if ``obj`` is a numpy.ndarray, otherwise ``False``.
    """
    return isinstance(obj, numpy.ndarray)
