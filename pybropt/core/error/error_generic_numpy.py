import numpy

from . import generic_default_cond

def generic_check_dtype_issubdtype(v, vname, vdtype):
    """
    Generic check ndarray subdtype function.

    Parameters
    ----------
    v : dtype
        Reference to dtype variable.
    vname : str
        Name associated with the dtype variable.
    vdtype : dtype, str, tuple
        type : numpy.dtype or dtype string associated with the object variable.
        tuple : Tuple of numpy.dtype or dtype strings (logical or) for the object variable.
    """
    istuple = isinstance(vdtype, tuple)
    logic = any(numpy.issubdtype(v, e) for e in vdtype) if istuple else numpy.issubdtype(v, vdtype)
    if not logic:
        mname = ""
        if istuple:
            l = len(vdtype)
            for i in range(l):
                mname += str(vdtype[i])
                if i < l - 2:
                    mname += ", "
                elif i < l - 1:
                    mname += ", or "
        else:
            mname = str(vdtype)
        raise TypeError("variable '{0}' must be of dtype {1}".format(vname, mname))

def generic_cond_check_dtype_issubdtype(v, vname, vdtype, cond = generic_default_cond):
    """
    Generic check ndarray subdtype function.

    Parameters
    ----------
    v : dtype
        Reference to dtype variable.
    vname : str
        Name associated with the dtype variable.
    vdtype : dtype, str, tuple
        type : numpy.dtype or dtype string associated with the object variable.
        tuple : Tuple of numpy.dtype or dtype strings (logical or) for the object variable.
    """
    if cond(v):
        generic_check_dtype_issubdtype(v, vname, vdtype)

################################################################################
########################### ndarray check functions ############################
################################################################################
def generic_check_ndarray_dtype_issubdtype(v, vname, vdtype):
    """
    Generic check ndarray subdtype function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vdtype : dtype, str, tuple
        type : numpy.dtype or dtype string associated with the object variable.
        tuple : Tuple of numpy.dtype or dtype strings (logical or) for the object variable.
    """
    istuple = isinstance(vdtype, tuple)
    logic = any(numpy.issubdtype(v.dtype, e) for e in vdtype) if istuple else numpy.issubdtype(v.dtype, vdtype)
    if not logic:
        mname = ""
        if istuple:
            l = len(vdtype)
            for i in range(l):
                mname += str(vdtype[i])
                if i < l - 2:
                    mname += ", "
                elif i < l - 1:
                    mname += ", or "
        else:
            mname = str(vdtype)
        raise TypeError("variable '{0}' must be of dtype {1}".format(vname, mname))

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

################################################################################
##################### conditional ndarray check functions ######################
################################################################################
def generic_cond_check_ndarray_dtype_issubdtype(v, vname, vdtype, cond = generic_default_cond):
    """
    Generic conditional check ndarray subdtype function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vdtype : dtype, str, tuple
        type : numpy.dtype or dtype string associated with the object variable.
        tuple : Tuple of numpy.dtype or dtype strings (logical or) for the object variable.
    cond : callable
    """
    if cond(v):
        generic_check_ndarray_dtype_issubdtype(v, vname, vdtype)

def generic_cond_check_ndarray_eq(v, vname, w, wname, cond = generic_default_cond):
    """
    Generic check ndarray value function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    w : object, tuple
        object : any object or primitive to be tested for equality.
        tuple : tuple of objects or primitives (logical or) to be tested for equality.
    wname : str, tuple
        Name associated with 'w'.
    cond : callable
    """
    if cond(v):
        generic_check_ndarray_eq(v, vname, w, wname)

def generic_cond_check_ndarray_ndim(v, vname, vndim, cond = generic_default_cond):
    """
    Generic check ndarray dimension function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vndim : int, tuple
        int : number of dimensions ndarray must have.
        tuple : tuple of number of dimensions (logical or) ndarray must have.
    cond : callable
    """
    if cond(v):
        generic_check_ndarray_ndim(v, vname, vndim)

def generic_cond_check_ndarray_size(v, vname, vsize, cond = generic_default_cond):
    """
    Generic check ndarray size function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vsize : int, tuple
        int : size ndarray must have.
        tuple : tuple of sizes (logical or) ndarray must have.
    cond : callable
    """
    if cond(v):
        generic_check_ndarray_size(v, vname, vsize)

def generic_cond_check_ndarray_sum(v, vname, vsum, vaxis = None, cond = generic_default_cond):
    """
    Generic check ndarray sum function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vsum : number, tuple
        number : sum ndarray must have.
        tuple : tuple of sums (logical or) ndarray must have.
    vaxis : int
        Axis or axes along which to sum.
    """
    if cond(v):
        generic_cond_check_ndarray_sum(v, vname, vsum, vaxis)

def generic_cond_check_ndarray_shape(v, vname, vshape, vaxis = None, cond = generic_default_cond):
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
    if cond(v):
        generic_check_ndarray_shape(v, vname, vshape, vaxis)
