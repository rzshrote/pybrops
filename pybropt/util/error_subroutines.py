# import 3rd party modules
import numpy

def check_is_matrix(mat, varname):
    if not isinstance(mat, numpy.ndarray):
        raise TypeError("'%s' must be a numpy.ndarray." % varname)

def check_is_string(s, varname):
    if not isinstance(s, str):
        raise TypeError("'%s' must be a string." % varname)

def check_is_list(l, varname):
    if not isinstance(l, list):
        raise TypeError("'%s' must be a list." % varname)

def check_is_tuple(t, varname):
    if not isinstance(t, tuple):
        raise TypeError("'%s' must be a tuple." % varname)

def check_is_bool(b, varname):
    if not isinstance(b, bool):
        raise TypeError("'%s' must be a bool." % varname)

def check_is_iterable(i, varname):
    if not hasattr(i, "__iter__"):
        raise TypeError("'%s' must be iterable." % varname)

def check_is_tuple_or_list(v, varname):
    if not isinstance(v, (tuple,list)):
        raise TypeError("'%s' must be a tuple or list." % varname)

def check_is_integer_or_floating_dtype(dtype, varname):
    if ((not numpy.issubdtype(dtype, numpy.integer)) and
        (not numpy.issubdtype(dtype, numpy.floating))):
        raise TypeError("'%s' must be an integer or floating dtype." % varname)

def check_is_integer_or_floating(var, varname):
    dtype_var = type(var)
    if ((not numpy.issubdtype(dtype_var, numpy.integer)) and
        (not numpy.issubdtype(dtype_var, numpy.floating))):
        raise TypeError("'%s' must be an integer or floating dtype." % varname)

def check_is_integer_or_inf(var, varname):
    if ((not numpy.issubdtype(type(var), numpy.integer)) and (var != numpy.inf)):
        raise TypeError("'%s' must be an integer dtype or infinity." % varname)

def check_is_string_or_object_dtype(dtype, varname):
    if ((not numpy.issubdtype(dtype, numpy.str_)) and
        (not numpy.issubdtype(dtype, numpy.object_))):
        raise TypeError("'%s' must be a string or object dtype." % varname)

def check_is_numeric_or_bool_dtype(dtype, varname):
    if ((not numpy.issubdtype(dtype, numpy.number)) and
        (not numpy.issubdtype(dtype, numpy.bool_))):
        raise TypeError("'%s' must be a numeric or bool dtype." % varname)

def check_is_numeric_dtype(dtype, varname):
    if not numpy.issubdtype(dtype, numpy.number):
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
    if ((not numpy.issubdtype(mat.dtype, numpy.number)) and
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

def check_matrix_dtype_is_object_(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.object_):
        raise ValueError("'%s' must be an object_ matrix." % varname)

def check_matrix_dtype_is_string_(mat, varname):
    if not numpy.issubdtype(mat.dtype, numpy.string_):
        raise ValueError("'%s' must be a string_ matrix." % varname)

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

def check_divisibility(a, aname, b, bname):
    if a % b != 0:
        raise ValueError(
            "'%s' (%d) is not divisible by '%s' (%d)." % (aname, a, bname, b)
        )

def check_is_dict(d, varname):
    if not isinstance(d, dict):
        raise TypeError("'%s' must be a dictionary." % varname)

def check_keys_in_dict(d, varname, *argv):
    key_absent = [arg not in d for arg in argv]
    if any(key_absent):
        # build error string
        err_str = "'%s' does not have the required fields.\n" % varname
        for i,absent in enumerate(key_absent):
            if absent:
                err_str += "    %s is missing.\n" % argv[i]
        raise ValueError(err_str)

def check_matrix_sum(mat, varname, sum):
    if mat.sum() != sum:
        raise ValueError("'%s' must have a sum of %s" % (varname, sum))

def check_is_not_none(var, varname):
    if var is None:
        raise ValueError("'%s' cannot be None." % varname)

def check_is_MarkerSet(mkrset, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.MarkerSet
    if not isinstance(mkrset, pybropt.popgen.MarkerSet):
        raise TypeError("'%s' must be a MarkerSet object." % varname)

def check_is_GeneticMap(gmap, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.GeneticMap
    if not isinstance(gmap, pybropt.popgen.GeneticMap):
        raise TypeError("'%s' must be a GeneticMap object." % varname)

def check_is_GenomicModel(gmod, varname):
    # HACK: elminates circular dependency issues
    import pybropt.model.GenomicModel
    if not isinstance(gmod, pybropt.model.GenomicModel):
        raise TypeError("'%s' must be a GenomicModel object." % varname)

def check_is_SearchSpace(sspace, varname):
    if not isinstance(sspace, SearchSpace):
        raise TypeError("'%s' must be a SearchSpace object." % varname)

def check_is_CategoricalSearchSpace(sspace, varname):
    if not isinstance(sspace, CategoricalSearchSpace):
        raise TypeError("'%s' must be a CategoricalSearchSpace object." % varname)

def check_is_SetSearchSpace(sspace, varname):
    if not isinstance(sspace, SetSearchSpace):
        raise TypeError("'%s' must be a SetSearchSpace object." % varname)

def check_is_Population(pop, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.Population
    if not isinstance(pop, pybropt.popgen.Population):
        raise TypeError("'%s' must be a Population object." % varname)

def check_is_ParametricGenomicModel(gmod, varname):
    # HACK: elminates circular dependency issues
    import pybropt.model.ParametricGenomicModel
    if not isinstance(gmod, pybropt.model.ParametricGenomicModel):
        raise TypeError("'%s' must be a ParametricGenomicModel object." % varname)

def check_is_NonparametricGenomicModel(gmod, varname):
    # HACK: elminates circular dependency issues
    import pybropt.model.NonparametricGenomicModel
    if not isinstance(gmod, pybropt.model.NonparametricGenomicModel):
        raise TypeError("'%s' must be a NonparametricGenomicModel object." % varname)

def check_is_Cross(cross, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.Cross
    if not isinstance(cross,pybropt.popgen.Cross):
        raise TypeError("'%s' must be a Cross object." % varname)

################################################################################

def error_readonly(varname):
    raise AttributeError("'%s' is read-only." % varname)

################################################################################

def cond_check_is_matrix(mat, varname, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_is_matrix(mat, varname)

def cond_check_is_string(s, varname, cond=(lambda s: s is not None)):
    if cond(s):
        check_is_string(s, varname)

def cond_check_is_MarkerSet(mkrset, varname, cond=(lambda mkrset: mkrset is not None)):
    if cond(mkrset):
        check_is_MarkerSet(mkrset, varname)

def cond_check_is_GeneticMap(gmap, varname, cond=(lambda gmap: gmap is not None)):
    if cond(gmap):
        check_is_GeneticMap(gmap, varname)

def cond_check_is_GenomicModel(gmod, varname, cond=(lambda gmod: gmod is not None)):
    if cond(gmod):
        check_is_GenomicModel(gmod, varname)

def cond_check_matrix_dtype_is_object_(mat, varname, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_dtype_is_object_(mat, varname)

def cond_check_matrix_ndim(mat, varname, ndim, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_ndim(mat, varname, ndim)

def cond_check_matrix_size(mat, varname, matsize, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_size(mat, varname, matsize)

def cond_check_matrix_dtype_is_string_(mat, varname, cond=(lambda mat: mat is not None)):
    if cond(mat):
        check_matrix_dtype_is_string_(mat, varname)
