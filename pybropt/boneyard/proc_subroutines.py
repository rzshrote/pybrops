import numpy
from .error_subroutines import *

def proc_hc_sa_set_varg(algorithm_varg, geno):
    # if algorithm_varg is None
    if algorithm_varg is None:
        # make a default dictionary
        algorithm_varg = {
            "states": numpy.arange(geno.shape[1])
        }
    # at this point, algorithm_varg must be a dict
    check_is_dict(algorithm_varg, "algorithm_varg")

    # if algorithm_varg does not have a 'states' key
    if "states" not in algorithm_varg.keys():
        # add a default 'states' key
        algorithm_varg["states"] = numpy.arange(geno.shape[1])

    return algorithm_varg

def proc_hc_sa_state_varg(algorithm_varg, geno, sel_size):
    if algorithm_varg is None:
        # make a default dictionary
        algorithm_varg = {
            "states": numpy.tile(numpy.arange(geno.shape[1]), (sel_size,)),
            "dim_sizes": numpy.repeat(geno.shape[1], sel_size)
        }
    # at this point, algorithm_varg must be a dict
    check_is_dict(algorithm_varg, "algorithm_varg")
    check_keys_in_dict(algorithm_varg, "algorithm_varg", "states", "dim_sizes")

    return algorithm_varg

def proc_icpso_varg(algorithm_varg, ):
    check_is_dict(algorithm_varg, "algorithm_varg")
    check_keys_in_dict(algorithm_varg, "algorithm_varg",
        "n", "states", "dim_sizes", "inertia_wt", "accel_coeff_pbest",
        "accel_coeff_gbest", "scale_factor", "stpfn"
    )
    return algorithm_varg

def proc_sel_size_pop_size(sel_size, pop_size, geno):

    type_sel_size = type(sel_size)# get type
    check_is_integer_or_floating_dtype(type_sel_size, "sel_size")

    if numpy.issubdtype(type(sel_size), numpy.floating):
        sel_size = numpy.int(numpy.around(sel_size))
    if sel_size > geno.shape[1]:
        raise ValueError("'sel_size' is greater than geno.shape[1].")

    return sel_size, pop_size
