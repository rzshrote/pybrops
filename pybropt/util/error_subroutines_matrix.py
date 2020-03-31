### Error subroutines specifically for numpy.ndarray objects

import numpy

def check_is_matrix(mat, varname):
    if not isinstance(mat, numpy.ndarray):
        raise TypeError("'%s' must be a numpy.ndarray." % varname)

def cond_check_is_matrix(mat, varname, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_is_matrix(mat, varname)

def check_matrix_all_value(mat, varname, value):
    if not numpy.all(mat == value):
        raise ValueError("All elments in '%s' must be %s." % (varname, value))

def check_matrix_axis_len(mat, varname, axis, axislen):
    if mat.shape[axis] != axislen:
        raise ValueError("'%s' must have length %d alonge axis %d." % (varname, axislen, axis))

def check_matrix_is_binary(mat, varname):
    if (mat.max() > 1) or (mat.min() < 0):
        raise ValueError("'%s' must be a binary matrix." % varname)

def check_matrix_ndim(mat, varname, ndim):
    if mat.ndim != ndim:
        raise ValueError("'%s' must have %d dimensions." % (varname, ndim))

def cond_check_matrix_ndim(mat, varname, ndim, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_ndim(mat, varname, ndim)

def check_matrix_size(mat, varname, matsize):
    if mat.size != matsize:
        raise ValueError("'%s' must have a size of %d." % (varname, matsize))

def cond_check_matrix_size(mat, varname, matsize, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_size(mat, varname, matsize)

def check_matrix_sum(mat, varname, sum):
    if mat.sum() != sum:
        raise ValueError("'%s' must have a sum of %s" % (varname, sum))


################################################################################
# dtype checking of numpy.ndarray

def check_matrix_dtype(mat, varname, dtype):
    if not mat.dtype == dtype:
        raise ValueError("'%s' must have a dtype of %s.", (varname, dtype))

def check_matrix_dtype_is_floating(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.floating):
        raise ValueError("'%s' must be an floating point matrix." % varname)

def check_matrix_dtype_is_integer(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.integer):
        raise ValueError("'%s' must be an integer matrix." % varname)

def check_matrix_dtype_is_numeric(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.number):
        raise ValueError("'%s' must be a numeric matrix." % varname)

def check_matrix_dtype_is_numeric_or_bool(mat, varname):
    if ((not numpy.issubdtype(mat.dtype, numpy.number)) and
        (not numpy.issubdtype(mat.dtype, numpy.bool_))):
        raise TypeError("'%s' must be a numeric or bool matrix." % varname)

def check_matrix_dtype_is_object_(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.object_):
        raise ValueError("'%s' must be an object_ matrix." % varname)

def cond_check_matrix_dtype_is_object_(mat, varname, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_dtype_is_object_(mat, varname)

def check_matrix_dtype_is_string_(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.string_):
        raise ValueError("'%s' must be a string_ matrix." % varname)

def cond_check_matrix_dtype_is_string_(mat, varname, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_dtype_is_string_(mat, varname)
