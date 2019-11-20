# import 3rd party libraries
import numpy
import time
import pandas
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
import objfn
import algo

def opv(geno,
        coeff,
        k,
        algo,
        algo_varg = None,
        seed = None,
        nthreads = 1,
        zwidth = 3,
        verbose = True):
    """
    Optimal Population Value (OPV) selection

    Parameters
    ==========
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.
    coeff : numpy.ndarray
        An array of coefficients for allele effects. The dtype of 'coeff'
        should be either 'float32' or 'float64'. This array should be single
        dimensional (column,) = (N,) where 'N' represents number of markers.
    k : int, numpy.integer, float, numpy.float
        Number of individuals to select OR proportion of individuals to select.
        If a proportion is given, round to the nearest individual.
    algo : {'hc_sa_set', 'hc_sa_state', 'icpso'}
        Specification for algorithm to use for selecting individuals.
        Algorithm overview:
            hc_sa_set       Greedy steepest ascent hillclimb for sets
            hc_sa_state     Greedy steepest ascent hillclimb for states
            icpso           Non-greedy integer/categorical PSO
        Algorithm speeds:
            hc_sa_set       Fast
            hc_sa_state     Medium
            icpso           Variable, but typically slow
        Algorithm solution quality:
            hc_sa_set       Okay, at least locally optimal
            hc_sa_state     Okay, at least locally optimal
            icpso           Good, at least locally optimal
    algo_varg : dictionary, None
        Dictionary of algorithm variable arguments. If None, will provide
        default options for hillclimb algorithms. Will throw error for icpso.
    seed : int
        Random number seed.
    nthreads : int
        Number of threads.
    zwidth : int
        Number of zeros to preceed in positions column names in the output
        pandas dataframes.
        Example:
            zwidth = 2:
                positions = ["x00", "x01", "x02", ...]
            zwidth = 3:
                positions = ["x000", "x001", "x002", ...]
    verbose : boolean
        Show verbose messages


    """
    ########################################
    # Step 0) process arguments
    ########################################

    # check for correct genotype shape
    if len(geno.shape) != 3:
        raise TypeError(
            "'geno' must have 3 dimensions: (phase, indiv, locus)\n"\
            "    geno.shape = %s" %
            (geno.shape,)
        )

    # check for correct coeff shape
    if len(coeff.shape) != 1:
        raise TypeError("'coeff' must have 1 dimension: (locus,)")

    # check if loci lengths are the same
    if geno.shape[2] != coeff.shape[0]:
        raise ValueError(
            "Length of 'coeff' does not match 'geno':\n"\
            "    geno.shape = (%s, %s, locus=%s)\n"\
            "    coeff.shape = (locus=%s,)\n"\
            "'locus' dimensions need to match." %
            (geno.shape[0], geno.shape[1], geno.shape[2], coeff.shape[0])
        )

    # process 'k' and make sure it's the right data type
    if isinstance(k, (int, numpy.integer)):
        k = k if k < geno.shape[1] else geno.shape[1]
    elif isinstance(k, (float, numpy.float)):
        k = numpy.int32(numpy.around(k))
        k = k if k < geno.shape[1] else geno.shape[1]
    else:
        raise TypeError(
            "'k' must be int, numpy.integer, float, or numpy.float\n"\
            "    type(k) = %s" %
            (type(k),)
        )

    # tests for algorithm choice
    if algo == 'hc_sa_set':
        # if algo_varg is None
        if isinstance(algo_varg, type(None)):
            # make a default dictionary
            algo_varg = {
                "states": numpy.arange(geno.shape[1])
            }
        # make sure we've received a dictionary
        elif not isinstance(algo_varg, dict):
            raise TypeError(
                "Expected 'algo_varg' to be a dict\n"\
                "    type(algo_varg) = %s" %
                (type(algo_varg))
            )
        # if algo_varg does not have a 'states' key
        elif "states" not in algo_varg.keys():
            # add a default 'states' key
            algo_varg["states"] = numpy.arange(geno.shape[1])
    elif algo == 'hc_sa_state':
        if isinstance(algo_varg, type(None)):
            # make a default dictionary
            algo_varg = {
                "states": numpy.tile(numpy.arange(geno.shape[1]), (k,1)),
                "dim_sizes": numpy.repeat(geno.shape[1], k)
            }
        # make sure we've received a dictionary
        elif not isinstance(algo_varg, dict):
            raise TypeError(
                "Expected 'algo_varg' to be a dict\n"\
                "    type(algo_varg) = %s" %
                (type(algo_varg))
            )
        elif "states" or "dim_sizes" not in algo_varg.keys():
            raise ValueError(
                "Expected 'algo_varg' dict to have fields "\
                "'states' and 'dim_sizes'"
            )
    elif algo == 'icpso':
        # make sure we've received a dictionary
        if not isinstance(algo_varg, dict):
            raise TypeError(
                "Expected 'algo_varg' to be a dict\n"\
                "    type(algo_varg) = %s" %
                (type(algo_varg))
            )
        # make sure we have all of our necessary arguments.
        elif ("n"                 not in algo_varg.keys() or
              "states"            not in algo_varg.keys() or
              "dim_sizes"         not in algo_varg.keys() or
              "inertia_wt"        not in algo_varg.keys() or
              "accel_coeff_pbest" not in algo_varg.keys() or
              "accel_coeff_gbest" not in algo_varg.keys() or
              "scale_factor"      not in algo_varg.keys() or
              "stpfn"             not in algo_varg.keys()):
            raise ValueError(
                "Expected 'algo_varg' to have all of the following fields:\n"\
                "    Field              Received\n"\
                "    n                  %s\n"\
                "    states             %s\n"\
                "    dim_sizes          %s\n"\
                "    inertia_wt         %s\n"\
                "    accel_coeff_pbest  %s\n"\
                "    accel_coeff_gbest  %s\n"\
                "    scale_factor       %s\n"\
                "    stpfn              %s\n" %
                ("n" in algo_varg.keys(),
                 "states" in algo_varg.keys(),
                 "dim_sizes" in algo_varg.keys(),
                 "inertia_wt" in algo_varg.keys(),
                 "accel_coeff_pbest" in algo_varg.keys(),
                 "accel_coeff_gbest" in algo_varg.keys(),
                 "scale_factor" in algo_varg.keys(),
                 "stpfn" in algo_varg.keys())
            )
    else:
        raise ValueError(
            "Expected 'algo' to be one of the following\n"\
            "    'hc_sa_set'\n"\
            "    'hc_sa_state'\n"\
            "    'icpso'"
        )

    ########################################
    # Step 1) establish engine states
    ########################################

    # if no seed is provided, use time as random seed
    if seed == None:
        seed = numpy.uint32(time.time())

        # print engine state before beginning.
        if verbose:
            print("OPV engine state:")
            print("Miscellaneous", "=============", sep='\n')
            print("numpy.random seed =", seed)
            print("Number of threads =", nthreads)

    ########################################
    # Step 2) begin computing
    ########################################

    # compute hcoeff
    opv_varg = {"hcoeff": geno * coeff}

    # define a variable for output by the algorithm
    algo_out = None

    # use the correct algorithm
    if algo == "hc_sa_set":
        algo_out = algo.hc_sa_set(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            k = k,
            **algo_varg,        # variable args: contains 'states'
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    elif algo == "hc_sa_state":
        algo_out = algo.hc_sa_state(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            **algo_varg,        # variable args: contains 'states', 'dim_sizes'
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    elif algo == "icpso":
        algo_out = algo.icpso(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            **algo_varg,        # variable args
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    
