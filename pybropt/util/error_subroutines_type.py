### Error subroutines specifically for data type checking

import numpy

################################################################################
# reliant on language itself
def check_is_bool(b, varname):
    if not isinstance(b, bool):
        raise TypeError("'%s' must be a bool." % varname)

def check_is_dict(d, varname):
    if not isinstance(d, dict):
        raise TypeError("'%s' must be a dictionary." % varname)

def check_is_iterable(i, varname):
    if not hasattr(i, "__iter__"):
        raise TypeError("'%s' must be iterable." % varname)

def check_is_list(l, varname):
    if not isinstance(l, list):
        raise TypeError("'%s' must be a list." % varname)

def check_is_not_none(var, varname):
    if var is None:
        raise ValueError("'%s' cannot be None." % varname)

def check_is_string(s, varname):
    if not isinstance(s, str):
        raise TypeError("'%s' must be a string." % varname)

def cond_check_is_string(s, varname, cond=(lambda s: s is not None)):
    if cond(s):
        check_is_string(s, varname)

def check_is_tuple(t, varname):
    if not isinstance(t, tuple):
        raise TypeError("'%s' must be a tuple." % varname)

def check_is_tuple_or_list(v, varname):
    if not isinstance(v, (tuple,list)):
        raise TypeError("'%s' must be a tuple or list." % varname)

################################################################################
# class type check functions

############################################################
# pybropt.algo check functions
def check_is_Algorithm(algo, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.Algorithm
    if not isinstance(algo, pybropt.algo.Algorithm):
        raise TypeError("'%s' must be an Algorithm object." % varname)

def check_is_CategoricalSearchSpace(sspace, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.CategoricalSearchSpace
    if not isinstance(sspace, pybropt.algo.CategoricalSearchSpace):
        raise TypeError("'%s' must be a CategoricalSearchSpace object." % varname)

def check_is_ContinuousSearchSpace(sspace, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.ContinuousSearchSpace
    if not isinstance(sspace, pybropt.algo.ContinuousSearchSpace):
        raise TypeError("'%s' must be a ContinuousSearchSpace object." % varname)

def check_is_HillClimber(algo, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.HillClimber
    if not isinstance(algo, pybropt.algo.HillClimber):
        raise TypeError("'%s' must be a HillClimber object." % varname)

def check_is_ICPSO(algo, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.ICPSO
    if not isinstance(algo, pybropt.algo.ICPSO):
        raise TypeError("'%s' must be an ICPSO object." % varname)

def check_is_ParticleSwarmOptimization(algo, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.ParticleSwarmOptimization
    if not isinstance(algo, pybropt.algo.ParticleSwarmOptimization):
        raise TypeError("'%s' must be a ParticleSwarmOptimization object." % varname)

def check_is_SearchSpace(sspace, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.SearchSpace
    if not isinstance(sspace, pybropt.algo.SearchSpace):
        raise TypeError("'%s' must be a SearchSpace object." % varname)

def check_is_SetHC(algo, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.SetHC
    if not isinstance(algo, pybropt.algo.SetHC):
        raise TypeError("'%s' must be a SetHC object." % varname)

def check_is_SetSearchSpace(sspace, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.SetSearchSpace
    if not isinstance(sspace, pybropt.algo.SetSearchSpace):
        raise TypeError("'%s' must be a SetSearchSpace object." % varname)

def check_is_StateHC(algo, varname):
    # HACK: elminates circular dependency issues
    import pybropt.algo.StateHC
    if not isinstance(algo, pybropt.algo.StateHC):
        raise TypeError("'%s' must be a StateHC object." % varname)

############################################################
# pybropt.breed check functions

############################################################
# pybropt.model check functions
def check_is_GenomicModel(gmod, varname):
    # HACK: elminates circular dependency issues
    import pybropt.model.GenomicModel
    if not isinstance(gmod, pybropt.model.GenomicModel):
        raise TypeError("'%s' must be a GenomicModel object." % varname)

def cond_check_is_GenomicModel(gmod, varname, cond=(lambda gmod: gmod is not None)):
    if cond(gmod):
        check_is_GenomicModel(gmod, varname)

def check_is_NonparametricGenomicModel(gmod, varname):
    # HACK: elminates circular dependency issues
    import pybropt.model.NonparametricGenomicModel
    if not isinstance(gmod, pybropt.model.NonparametricGenomicModel):
        raise TypeError("'%s' must be a NonparametricGenomicModel object." % varname)

def check_is_ParametricGenomicModel(gmod, varname):
    # HACK: elminates circular dependency issues
    import pybropt.model.ParametricGenomicModel
    if not isinstance(gmod, pybropt.model.ParametricGenomicModel):
        raise TypeError("'%s' must be a ParametricGenomicModel object." % varname)

############################################################
# pybropt.popgen check functions
def check_is_Cross(cross, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.Cross
    if not isinstance(cross,pybropt.popgen.Cross):
        raise TypeError("'%s' must be a Cross object." % varname)

def check_is_GeneticMap(gmap, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.GeneticMap
    if not isinstance(gmap, pybropt.popgen.GeneticMap):
        raise TypeError("'%s' must be a GeneticMap object." % varname)

def cond_check_is_GeneticMap(gmap, varname, cond=(lambda gmap: gmap is not None)):
    if cond(gmap):
        check_is_GeneticMap(gmap, varname)

def check_is_MarkerSet(mkrset, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.MarkerSet
    if not isinstance(mkrset, pybropt.popgen.MarkerSet):
        raise TypeError("'%s' must be a MarkerSet object." % varname)

def cond_check_is_MarkerSet(mkrset, varname, cond=(lambda mkrset: mkrset is not None)):
    if cond(mkrset):
        check_is_MarkerSet(mkrset, varname)

def check_is_Population(pop, varname):
    # HACK: elminates circular dependency issues
    import pybropt.popgen.Population
    if not isinstance(pop, pybropt.popgen.Population):
        raise TypeError("'%s' must be a Population object." % varname)

################################################################################
# reliant on numpy

############################################################
# direct dtype check functions
def check_is_integer_or_floating_dtype(dtype, varname):
    if ((not numpy.issubdtype(dtype, numpy.integer)) and
        (not numpy.issubdtype(dtype, numpy.floating))):
        raise TypeError("'%s' must be an integer or floating dtype." % varname)

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

def check_is_integer_dtype(dtype, varname):
    if not numpy.issubdtype(dtype, numpy.integer):
        raise ValueError("'%s' must be an integer dtype." % varname)

def check_is_floating_dtype(dtype, varname):
    if not numpy.issubdtype(dtype, numpy.floating):
        raise TypeError("'%s' must be a floating dtype." % varname)

############################################################
# indirect dtype check functions
def check_is_integer_or_inf(var, varname):
    if ((not numpy.issubdtype(type(var), numpy.integer)) and (var != numpy.inf)):
        raise TypeError("'%s' must be an integer dtype or infinity." % varname)

def check_is_integer_or_floating(var, varname):
    dtype_var = type(var)
    if ((not numpy.issubdtype(dtype_var, numpy.integer)) and
        (not numpy.issubdtype(dtype_var, numpy.floating))):
        raise TypeError("'%s' must be an integer or floating dtype." % varname)

def check_is_numeric(n, varname):
    if not numpy.issubdtype(type(n), numpy.number):
        raise TypeError("'%s' must be a numeric type." % varname)

def check_is_integer(i, varname):
    if not numpy.issubdtype(type(i), numpy.integer):
        raise TypeError("'%s' must be an integer type." % varname)

def check_is_floating(f, varname):
    if not numpy.issubdtype(type(f), numpy.floating):
        raise TypeError("'%s' must be a floating type." % varname)
