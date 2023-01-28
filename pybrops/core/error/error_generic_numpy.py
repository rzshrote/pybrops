import numpy

################################################################################
########################### ndarray check functions ############################
################################################################################
def generic_check_ndarray_eq(v, vname, w, wname):
    """
    Generic check ndarray value function.

    Parameters
    ----------
    v : numpy.ndarray
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    w : object
        Any object or primitive to be tested for equality.
    wname : str
        Name associated with 'w'.
    """
    if not numpy.all(v == w):
        raise ValueError("variable '{0}' must have values equal to {1}".format(vname, wname))

def generic_check_ndarray_ndim(v, vname, vndim):
    """
    Generic check ndarray dimension function.

    Parameters
    ----------
    v : numpy.ndarray
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vndim : int
        Number of dimensions ndarray must have.
    """
    if v.ndim != vndim:
        raise ValueError("variable '{0}' must have dimension equal to {1}".format(vname, vndim))

def generic_check_ndarray_ndim_gteq(v, vname, vndim):
    """
    Generic check ndarray dimension number greater than or equal to function.

    Parameters
    ----------
    v : numpy.ndarray
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vndim : int
        Minimum number of dimensions ndarray must have.
    """
    if v.ndim < vndim:
        raise ValueError("variable '{0}' must have dimension greater than or equal to {1}".format(vname, vndim))

def generic_check_ndarray_size(v, vname, vsize):
    """
    Generic check ndarray size function.

    Parameters
    ----------
    v : numpy.ndarray
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vsize : int
        Size ndarray must have.
    """
    if v.size != vsize:
        raise ValueError("variable '{0}' must have size equal to {1}".format(vname, vsize))

def generic_check_ndarray_sum(v, vname, vsum, vaxis = None):
    """
    Generic check ndarray sum function.

    Parameters
    ----------
    v : numpy.ndarray
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vsum : number, numpy.ndarray, list-like
        Sum ndarray must have.
    vaxis : int
        Axis or axes along which to sum.
    """
    if numpy.any(v.sum(vaxis) != vsum):
        raise ValueError("variable '{0}' must have sum equal to {1} along axis {2}".format(vname, vsum, vaxis))

def generic_check_ndarray_shape(v, vname, vshape, vaxis = None):
    """
    Generic check ndarray shape function.

    Parameters
    ----------
    v : numpy.ndarray
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vshape : int, tuple
        Shape 'v' must have.
    vaxis : None, int, tuple
        Axis or axes along which to check shape.
    """
    error = True
    if vaxis is None:
        error = (v.shape != vshape)
    elif isinstance(vaxis, int):
        error = (v.shape[vaxis] != vshape)
    elif isinstance(vaxis, tuple):
        error = any(v.shape[e] != vshape[i] for e,i in enumerate(vaxis))
    else:
        raise TypeError("vaxis must be of type None, int, or tuple")
    if error:
        raise ValueError("variable '{0}' must have shape equal to {1} along axis {2}".format(vname, vshape, vaxis))

def generic_check_ndarray_is_square(v: Any, vname: str) -> None:
    """
    Generic check that an ndarray is square.

    Parameters
    ----------
    v : numpy.ndarray
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    """
    s = v.shape # get shape
    if any(s[0] != e for e in s):
        raise ValueError("variable '{0}' must have equal lengths along all axes".format(vname))
