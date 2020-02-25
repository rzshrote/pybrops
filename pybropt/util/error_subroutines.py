import numpy

def check_is_matrix(mat, varname):
    if not isinstance(mat, numpy.ndarray):
        raise TypeError("'%s' must be a numpy.ndarray." % varname)

def check_is_numeric_or_bool_dtype(dtype, varname):
    if ((not numpy.issubdtype(dtype, numpy.numeric)) and
        (not numpy.issubdtype(dtype, numpy.bool_))):
        raise TypeError("'%s' must be a numeric or bool dtype." % varname)

def check_is_numeric_dtype(dtype, varname):
    if not numpy.issubdtype(dtype, numpy.numeric):
        raise ValueError("'%s' must be a numeric dtype." % varname)

def check_is_numeric(n, varname):
    if not numpy.issubdtype(type(n), numpy.number):
        raise TypeError("'%s' must be a numeric type." % varname)

def check_is_integer_dtype(dtype, varname):
    if not numpy.issubdtype(dtype, numpy.integer):
        raise ValueError("'%s' must be an integer dtype." % varname)

def check_is_integer(i, varname):
    if not numpy.issubdtype(type(i), numpy.integer):
        raise TypeError("'%s' must be an integer type." % varname)

def check_is_floating_dtype(dtype, varname):
    if not numpy.issubdtype(dtype, numpy.floating):
        raise TypeError("'%s' must be a floating dtype." % varname)

def check_is_floating(f, varname):
    if not numpy.issubdtype(type(f), numpy.floating):
        raise TypeError("'%s' must be a floating type." % varname)

def check_matrix_dtype_is_numeric_or_bool(mat, varname):
    if ((not numpy.issubdtype(mat.dtype, numpy.numeric)) and
        (not numpy.issubdtype(mat.dtype, numpy.bool_))):
        raise TypeError("'%s' must be a numeric or bool matrix." % varname)

def check_matrix_dtype_is_integer(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.integer):
        raise ValueError("'%s' must be an integer matrix." % varname)

def check_matrix_dtype_is_floating(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.floating):
        raise ValueError("'%s' must be an floating point matrix." % varname)

def check_matrix_dtype_is_numeric(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.number):
        raise ValueError("'%s' must be a numeric matrix." % varname)

def check_matrix_dtype(mat, varname, dtype):
    if not mat.dtype == dtype:
        raise ValueError("'%s' must have a dtype of %s.", (varname, dtype))

def check_matrix_ndim(mat, varname, ndim):
    if mat.ndim != ndim:
        raise ValueError("'%s' must have %d dimensions." % (varname, ndim))

def check_matrix_is_binary(mat, varname):
    if (mat.max() > 1) or (mat.min() < 0):
        raise ValueError("'%s' must be a binary matrix." % varname)

def check_matrix_size(mat, varname, matsize):
    if mat.size != matsize:
        raise ValueError("'%s' must have a size of %d." % (varname, matsize))

def check_matrix_axis_len(mat, varname, axis, axislen):
    if mat.shape[axis] != axislen:
        raise ValueError("'%s' must have length %d alonge axis %d." % (varname, axislen, axis))

def check_matrix_all_value(mat, varname, value):
    if not numpy.all(mat == value):
        raise ValueError("All elments in '%s' must be %s." % (varname, value))

def error_readonly(varname):
    raise AttributeError("'%s' is read-only." % varname)

def check_divisibility(a, aname, b, bname):
    if a % b != 0:
        raise ValueError(
            "'%s' (%d) is not divisible by '%s' (%d)." % (aname, a, bname, b)
        )
