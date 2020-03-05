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
import util
from util.error_subroutines import *
from util.proc_subroutines import *
from popgen.simul.breeding_program import breeding_program
from popgen.simul import mate

def opv(geno, coeff, sel_size, algorithm,
    algorithm_varg = None, seed = None, nthreads = 1, zwidth = 3,
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
    sel_size : int, numpy.integer, float, numpy.float
        Number of individuals to select OR proportion of individuals to select.
        If a proportion is given, round to the nearest individual.
    algorithm : {'hc_sa_set', 'hc_sa_state', 'icpso'}
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
    algorithm_varg : dictionary, None
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

    check_matrix_ndim(geno, "geno", 3)# check for correct genotype shape
    check_matrix_ndim(coeff, "coeff", 1)# check for correct coeff shape
    check_matrix_axis_len(geno, "geno", 2, len(coeff))# check if loci lengths are the same
    sel_size, pop_size = proc_sel_size_pop_size(sel_size, geno.shape[1], geno)

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
    # Step 2) compute results
    ########################################

    # compute hcoeff
    opv_varg = {"hcoeff": geno * coeff}

    # define a variable for output by the algorithm
    algo_out = None

    # tests for algorithm choice
    if algorithm == 'hc_sa_set':
        algorithm_varg = proc_hc_sa_set_varg(algorithm_varg, geno)
        algo_out = algo.hc_sa_set(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            k = sel_size,
            **algorithm_varg,        # variable args: contains 'states'
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    elif algorithm == 'hc_sa_state':
        algorithm_varg = proc_hc_sa_state_varg(algorithm_varg, geno, sel_size)
        algo_out = algo.hc_sa_state(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            **algorithm_varg,        # variable args: contains 'states', 'dim_sizes'
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    elif algorithm == 'icpso':
        algorithm_varg = proc_icpso_varg(algorithm_varg)
        algo_out = algo.icpso(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            **algorithm_varg,        # variable args
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    else:
        raise ValueError(
            "Expected 'algorithm' to be one of the following\n"\
            "    'hc_sa_set'\n"\
            "    'hc_sa_state'\n"\
            "    'icpso'"
        )

    ########################################
    # Step 3) create a results dataframe
    ########################################

    return algo_out



# simulation function
def opv_sim(geno, coeff, pop_size, sel_size, cycles, gmap, lgroup, algorithm,
    algorithm_varg = None, interference = None, seed = None, nthreads = 1,
    zwidth = 3, verbose = True):
    """
    Optimal Population Value (OPV) selection simulation

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
    pop_size : int
        Number of individuals to constitute a population to run selection on.
    sel_size : int, numpy.integer, float, numpy.float
        Number of individuals to select OR proportion of individuals to select.
        If a proportion is given, round to the nearest individual.
    cycles : int
        Number of breeding cycles to simulate.
    d : numpy.ndarray
        Genetic map in Morgan units.
    lgroup : numpy.ndarray
        Array of linkage group sizes. The sum of the elements should equal the
        length of 'd'.
    algorithm : {'hc_sa_set', 'hc_sa_state', 'icpso'}
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
    algorithm_varg : dictionary, None
        Dictionary of algorithm variable arguments. If None, will provide
        default options for hillclimb algorithms. Will throw error for icpso.
    interference : int
        Meiosis interference.
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
    # type checks
    check_is_matrix(geno, "geno")
    check_is_matrix(coeff, "coeff")
    check_is_integer(pop_size, "pop_size")
    check_is_integer_or_floating(sel_size, "sel_size")
    check_is_integer(cycles, "cycles")
    check_is_matrix(gmap, "gmap")
    check_is_matrix(lgroup, "lgroup")
    check_is_string(algorithm, "algorithm")

    # shape checks
    check_matrix_ndim(geno, "geno", 3)
    check_matrix_ndim(coeff, "coeff", 1)
    check_matrix_ndim(gmap, "gmap", 1)
    check_matrix_ndim(lgroup, "lgroup", 1)
    check_matrix_axis_len(geno, "geno", 2, len(coeff))
    check_matrix_sum(lgroup, "lgroup", len(gmap))

    # variable processing
    sel_size, pop_size = proc_sel_size_pop_size(sel_size, pop_size, geno)

    # tests for algorithm choice
    if algorithm == 'hc_sa_set':
        algorithm_varg = proc_hc_sa_set_varg(algorithm_varg, geno)
    elif algorithm == 'hc_sa_state':
        algorithm_varg = proc_hc_sa_state_varg(algorithm_varg, geno, sel_size)
    elif algorithm == 'icpso':
        algorithm_varg = proc_icpso_varg(algorithm_varg)
    else:
        raise ValueError(
            "Expected 'algorithm' to be one of the following\n"\
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

    # seed the rng
    numpy.random.seed(seed)

    # print engine state before beginning.
    if verbose:
        print("OPV siumlation engine state:")
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    # generate random seeds for algorithm seeding at each cycle
    algorithm_seed = numpy.random.randint(
        low = numpy.iinfo(numpy.uint32).min,
        high = numpy.iinfo(numpy.uint32).max,
        size = cycles,
        dtype = numpy.uint32
    )

    # generate random seeds for meiosis simulations
    meiosis_seed = numpy.random.randint(
        low = numpy.iinfo(numpy.uint32).min,
        high = numpy.iinfo(numpy.uint32).max,
        size = cycles,
        dtype = numpy.uint32
    )

    history = breeding_program(
        dtype_method = numpy.dtype("object_"),
        dtype_algorithm = numpy.dtype("object_"),
        dtype_seed = numpy.dtype("int64"),
        dtype_cycle = numpy.dtype("int64"),
        dtype_score = numpy.dtype("float64"),
        dtype_gebv = numpy.dtype("float64"),
        dtype_phase = numpy.dtype("float64")
    )

    # declare a working copy for genotype data
    pop_geno = geno.copy()
    pop_gebv = numpy.dot(pop_geno, coeff).sum(0)
    pop_score = objfn.opv(slice(None), pop_geno * coeff)

    history.history_add_population(
        method = ["opv"],
        algorithm = [algorithm],
        seed = [seed],
        cycle = [0],
        score = numpy.repeat(pop_score, len(pop_gebv)),
        gebv = pop_gebv,
        phase1 = pop_geno[0,:,:],
        phase2 = pop_geno[1,:,:]
    )

    for cycle in range(cycles):
        # declare opv_varg
        opv_varg = {"hcoeff": pop_geno * coeff}

        # declare optimization output variable
        optout = None

        # do the algorithm
        if algorithm == "hc_sa_set":
            optout = algo.hc_sa_set(
                objfn = objfn.opv,
                objfn_varg = opv_varg,
                k = sel_size,
                **algorithm_varg,
                seed = algorithm_seed[cycle],
                nthreads = nthreads,
                verbose = False
            )
        elif algorithm == "hc_sa_state":
            optout = algo.hc_sa_state(
                objfn = objfn.opv,
                objfn_varg = opv_varg,
                **algorithm_varg,
                seed = algorithm_seed[cycle],
                nthreads = nthreads,
                verbose = False
            )
        elif algorithm == "icpso":
            optout = algo.icpso(
                objfn = objfn.opv,
                objfn_varg = opv_varg,
                **algorithm_varg,
                seed = algorithm_seed[cycle],
                nthreads = nthreads,
                verbose = False
            )

        best = optout.gbest()

        selection_indices = best[4] if algorithm == "icpso" else best[2]
        selection_score = best[1]

        print(selection_indices)

        history.history_add_selection(
            method = ["opv"],
            algorithm = [algorithm],
            seed = [algorithm_seed[cycle]],
            cycle = [cycle+1],
            score = numpy.repeat(selection_score, len(selection_indices)),
            gebv = pop_gebv[selection_indices],
            phase1 = pop_geno[0,selection_indices,:],
            phase2 = pop_geno[1,selection_indices,:]
        )

        # mate and recalculate gebv, score
        pop_geno = mate.erand_dh(
            male = selection_indices,
            female = selection_indices,
            geno = pop_geno,
            gmap = gmap,
            lgroup = lgroup,
            n = pop_size // sel_size
        )
        pop_gebv = numpy.dot(pop_geno, coeff).sum(0)
        pop_score = objfn.opv(slice(None), pop_geno * coeff)

        history.history_add_population(
            method = ["opv"],
            algorithm = [algorithm],
            seed = [seed],
            cycle = [cycle+1],
            score = numpy.repeat(pop_score, len(pop_gebv)),
            gebv = pop_gebv,
            phase1 = pop_geno[0,:,:],
            phase2 = pop_geno[1,:,:]
        )

        # print progress
        if verbose:
            print("Cycle:", cycle+1)

    ########################################
    # Step 3) return history object
    ########################################

    return history
