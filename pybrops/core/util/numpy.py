import numpy

__all__ = ["is_ndarray"]

def is_ndarray(obj):
    return isinstance(obj, numpy.ndarray)
