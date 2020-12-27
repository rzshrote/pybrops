# import 3rd party libraries
import numpy
import time
import pandas
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
from util.human2bytes import human2bytes
from popgen.simul.meiosis import meiosis
from util.arg_proc import proc_pop_sel_sizes
from util.arg_proc import proc_geno_coeff
from util.arg_proc import proc_wcoeff
from util.arg_proc import proc_tfreq
from util.arg_proc import proc_gmap_lgroup
from util.arg_proc import proc_dcoeff
from util.arg_proc import proc_dcoeff_t
from util.arg_proc import proc_bcycles
from util.arg_proc import proc_mem
from util.arg_proc import proc_hc_sa_set_varg
from util.arg_proc import proc_hc_sa_state_varg
from util.arg_proc import proc_icpso_varg
from popgen.popvar import dh
from popgen.simul import cross
import objfn
import algo

# define a lookup dictionary for mapping function strings to function pointers
MOGS_ALGORITHM_DICT = {
    'hc_sa_set': algo.hc_sa_set,
    'hc_sa_state': algo.hc_sa_state,
    'icpso': algo.icpso
}

def mogm(sel_size, geno, coeff, tfreq, gmap, lgroup, dcoeff, algorithm,
         varAfn = '2-way', mapfn = 'haldane', mem = None,
         dtype = numpy.dtype("float64"),
         algorithm_varg = None, seed = None, nthreads = 1, zwidth = 3,
         verbose = True, verbose_algo = False):
    """
    Multi-Objective Genomic Mating (MOGM)

    Parameters
    ==========
    sel_size : int, numpy.integer, float, numpy.float
        Number of individuals to select OR proportion of individuals to select.
        If a proportion is given, round to the nearest individual.
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
    coeff : numpy.ndarray
        An array of coefficients for allele effects. The dtype of 'coeff'
        should be a floating point. This array should be single
        dimensional (column,) = (N,) where 'N' represents number of markers.
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
    verbose_algo : boolean
        Show algorithm verbose messages
    """
    ############################################################################
    # Step 0) process arguments: check for data entry errors, data type, etc.
    ############################################################################
    # process 'sel_size'
    pop_size, sel_size = proc_pop_sel_sizes(geno.shape[1], sel_size, geno)

    geno, coeff = proc_geno_coeff(geno, coeff)      # process 'geno' && 'coeff'
    gmap, lgroup = proc_gmap_lgroup(gmap, lgroup)   # process 'gmap' && 'lgroup'
    wcoeff = proc_wcoeff(None, coeff)               # calculate 'wcoeff'
    tfreq = proc_tfreq(tfreq)                       # process 'tfreq'
    dcoeff = proc_dcoeff(dcoeff, 3)                 # process 'dcoeff'
    mem = proc_mem(mem, gmap)                       # process mem

    # calculate all pairwise cross variances
    varA = dh.varA_2way(geno, coeff, gmap, lgroup, mem=mem)

    # test for algorithm choice and process 'algorithm_varg' accordingly
    if algorithm == 'hc_sa_set':
        algorithm_varg = proc_hc_sa_set_varg(algorithm_varg, sel_size, geno)
    elif algorithm == 'hc_sa_state':
        algorithm_varg = proc_hc_sa_state_varg(algorithm_varg, sel_size, geno)
    elif algorithm == 'icpso':
        algorithm_varg = proc_icpso_varg(algorithm_varg)
    else:
        raise ValueError(
            "Expected 'algorithm' to be one of the following\n"\
            "    'hc_sa_set'\n"\
            "    'hc_sa_state'\n"\
            "    'icpso'\n"\
            "Received: %s" % algorithm
        )

    ############################################################################
    # Step 1) establish engine states
    ############################################################################

    # if no seed is provided, use time as random seed
    if seed is None:
        seed = numpy.uint32(time.time())

    # print engine state before beginning.
    if verbose:
        print("MOGM engine state:")
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ########################################
    # Step 2) compute results
    ########################################

    # assemble function variable arguments
    mogm_varg = {
        "geno": geno,
        "wcoeff": wcoeff,
        "tfreq": tfreq,
        "varA": varA,
        "varAfn": varAfn,
        "dcoeff": dcoeff,
        "dtype": dtype
    }

    # lookup the algorithm in the dict, and execute the function with arguments
    algo_out = MOGS_ALGORITHM_DICT[algorithm](
        objfn = objfn.mogm_max,
        objfn_varg = mogm_varg,
        **algorithm_varg,
        seed = seed,
        nthreads = nthreads,
        verbose = verbose_algo
    )

    ########################################
    # Step 3) create a results dataframe
    ########################################

    # declare several dataframe variables
    info_df = None
    scores_df = None
    positions_df = None
    results_df = None
    history_df = None

    # create an info dataframe for storing algorithm details, etc.
    if algorithm == "hc_sa_set" or algorithm == "hc_sa_state":
        info_df = pandas.DataFrame(
            data = [["mogm",algorithm,seed]],
            columns = ["method","algorithm","seed"]
        )
        scores_df = pandas.DataFrame(
            data = [algo_out["X_gbest_scr"]],
            columns = ["score"]
        )
        positions_df = pandas.DataFrame(
            data = algo_out["X_gbest_pos"][None,:],
            columns = ["x"+str(i).zfill(zwidth)
                       for i in range(len(algo_out["X_gbest_pos"]))]
        )
    elif algorithm == "icpso":
        info_df = pandas.DataFrame(
            data = [["mogm",algorithm,seed]
                    for _ in range(algo_out["X_pbest_smpl"].shape[0])],
            columns = ["method","algorithm","seed"]
        )
        scores_df = pandas.DataFrame(
            data = algo_out["X_pbest_scr"],
            columns = ["score"]
        )
        positions_df = pandas.DataFrame(
            data = algo_out["X_pbest_smpl"],
            columns = ["x"+str(i).zfill(zwidth)
                       for i in range(algo_out["X_pbest_smpl"].shape[1])]
        )

    # concatenate results together
    results_df = pandas.concat(
        objs = [info_df, scores_df, positions_df],
        axis = 1
    )

    # handle the history
    history_df = pandas.DataFrame(
        data = algo_out["history"],
        columns = ["iter", "x_pos", "score"]
    )
    tmp = pandas.DataFrame(
        data = history_df["x_pos"].tolist(),
        columns = ["x"+str(i).zfill(zwidth)
                   for i in range(len(history_df["x_pos"][0]))]
    )
    history_df = pandas.concat(
        objs = [history_df, tmp],
        axis = 1
    )
    # drop the x_pos column
    history_df = history_df.drop(columns=["x_pos"])

    # return results and history
    return results_df, history_df



# simulation function
def mogm_sim(pop_size, sel_size, geno, coeff, tfreq, gmap, lgroup, dcoeff_t,
             bcycles, algorithm,
             varAfn = '2-way', mapfn = 'haldane', mem = None,
             dtype = numpy.dtype("float64"),
             matefn = None, algorithm_varg = None, interference = None,
             seed = None, nthreads = 1, zwidth = 3, verbose = True,
             verbose_algo = False):
    """
    Population Architect (PA) selection simulation

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
    gmap : numpy.ndarray
        Genetic map in Morgan units.
    lgroup : numpy.ndarray
        Array of linkage group sizes. The sum of the elements should equal the
        length of 'gmap'.
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
    ############################################################################
    # Step 0) process arguments
    ############################################################################
    pop_size, sel_size = proc_pop_sel_sizes(pop_size, sel_size, geno)
    geno, coeff = proc_geno_coeff(geno, coeff)      # process 'geno' && 'coeff'
    gmap, lgroup = proc_gmap_lgroup(gmap, lgroup)   # process 'gmap' && 'lgroup'
    wcoeff = proc_wcoeff(None, coeff)               # process 'wcoeff'
    tfreq = proc_tfreq(tfreq)                       # process 'tfreq'
    dcoeff_t = proc_dcoeff_t(dcoeff_t, bcycles, 3)  # process 'dcoeff_t'
    bcycles = proc_bcycles(bcycles)                 # process 'bcycles'
    mem = proc_mem(mem, gmap)                       # process mem

    # process 'algorithm'
    if algorithm == 'hc_sa_set':
        algorithm_varg = proc_hc_sa_set_varg(algorithm_varg, sel_size, geno)
    elif algorithm == 'hc_sa_state':
        algorithm_varg = proc_hc_sa_state_varg(algorithm_varg, sel_size, geno)
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
    if seed is None:
        seed = numpy.uint32(time.time())

    # seed the rng
    numpy.random.seed(seed)

    # print engine state before beginning.
    if verbose:
        print("MOGM simulation engine state:")
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random.seed =", seed)
        print("Number of threads =", nthreads)

    # generate random seeds for bcycle, shuffle, meiosis, for each bcycle
    seeds = numpy.random.randint(
        low = numpy.iinfo(numpy.uint32).min,
        high = numpy.iinfo(numpy.uint32).max,
        size = (2,bcycles),
        dtype = numpy.uint32
    )

    ########################################
    # Step 2) Initialize data collection lists
    ########################################
    # declare a working copy for genotype data, this will be modified by meiosis
    mogm_geno = geno.copy()

    # build method string
    method = "mogm"

    # calculate population GEBVs, put into list
    pop_gebvs = numpy.dot(mogm_geno, coeff).sum(0)

    # score the population using objfn.opv
    pop_opv_score = objfn.opv(slice(None), mogm_geno * coeff)

    # score the population using objfn.mogs
    # pop_mogm_score = objfn.mogm(
    #     numpy.arange(mogm_geno.shape[1]),
    #     **mogm_varg
    # )

    # declare lists to store results
    selections_list = list()
    population_list = [[
        method,         # method string
        "NA",           # f_PAU priority to generate this population
        "NA",           # f_PAFD priority to generate this population
        "NA",           # f_stdA to generate this population
        algorithm,      # algorithm string
        seed,           # seed for mogm_sim
        0,              # breeding cycle number (0 for start)
        pop_opv_score,  # upper selection limit score
        "NA",           # MOGM score for everything
        pop_gebvs.mean(),
        pop_gebvs.std(),
        pop_gebvs.min(),
        numpy.median(pop_gebvs),
        pop_gebvs.max()
    ] + pop_gebvs.tolist()]



    ########################################
    # Step 3) run simulations
    ########################################

    ############################################################################
    for bcycle in range(bcycles):
        # calculate all pairwise cross variances
        varA = dh.varA_2way(mogm_geno, coeff, gmap, lgroup, mem=mem)

        # build MOGM variable arguments
        mogm_varg = {
            "geno": mogm_geno,
            "wcoeff": wcoeff,
            "tfreq": tfreq,
            "varA" : varA,
            "varAfn": varAfn,
            "dcoeff": dcoeff_t[bcycle,:], # get row for current breeding cycle
            "dtype": dtype
        }

        # use algorithm selected, get output
        optout = MOGS_ALGORITHM_DICT[algorithm](
            objfn = objfn.mogm_max,
            objfn_varg = mogm_varg,
            **algorithm_varg,
            seed = seeds[0,bcycle],
            nthreads = nthreads,
            verbose = verbose_algo
        )

        sel_ix = None
        sel_score = None

        # process output
        if algorithm == "hc_sa_set" or algorithm == "hc_sa_state":
            sel_ix = optout["X_gbest_pos"]      # get sel_ix
            sel_score = optout["X_gbest_scr"]   # get sel_score
        elif algorithm == "icpso":
            sel_ix = optout["X_gbest_smpl"]     # get sel_ix
            sel_score = optout["X_gbest_scr"]   # get sel_score

        # calculate GEBVs, put into list
        sel_gebvs = numpy.dot(mogm_geno[:,sel_ix,:], coeff).sum(0)

        # calculate selection upper selection limit
        sel_opv_score = objfn.opv(sel_ix, mogm_geno * coeff)

        # store everything into selections_list
        selections_list.append([
            method,             # method string
            dcoeff_t[bcycle,0], # f_PAU priority to generate this selection
            dcoeff_t[bcycle,1], # f_PAFD priority to generate this selection
            dcoeff_t[bcycle,2], # f_stdA to generate this selection
            algorithm,          # algorithm string
            seed,               # seed for mogm_sim
            bcycle+1,           # breeding cycle number (0 for start)
            sel_opv_score,      # upper selection limit score
            sel_score,          # MOGM score for selections
            sel_gebvs.mean(),
            sel_gebvs.std(),
            sel_gebvs.min(),
            numpy.median(sel_gebvs),
            sel_gebvs.max()
        ] + sel_gebvs.tolist())

        # mate and make progeny
        mogm_geno = cross.f1_dh(
            male = sel_ix[0::2],
            female = sel_ix[1::2],
            geno = mogm_geno,
            gmap = gmap,
            lgroup = lgroup,
            n = 2 * pop_size // sel_size,
            seed = seeds[1,bcycle]
        )

        # calculate population GEBVs, put into list
        pop_gebvs = numpy.dot(mogm_geno, coeff).sum(0)

        # score the population using objfn.opv
        pop_opv_score = objfn.opv(slice(None), mogm_geno * coeff)

        # score the population using objfn.mogs
        # pop_mogm_score = objfn.mogm(
        #     numpy.arange(mogm_geno.shape[1]),
        #     **mogm_varg
        # )

        # append population stats to list
        population_list.append([
            method,             # method string
            dcoeff_t[bcycle,0], # f_PAU priority to generate this population
            dcoeff_t[bcycle,1], # f_PAFD priority to generate this population
            dcoeff_t[bcycle,2], # f_stdA to generate this population
            algorithm,          # algorithm string
            seed,               # seed for mogm_sim
            bcycle+1,           # breeding cycle number (0 for start)
            pop_opv_score,      # upper selection limit score
            "NA",               # MOGM score for everything
            pop_gebvs.mean(),
            pop_gebvs.std(),
            pop_gebvs.min(),
            numpy.median(pop_gebvs),
            pop_gebvs.max()
        ] + pop_gebvs.tolist())

        # print progress
        if verbose:
            print("Breeding cycle:", bcycle+1)
    ############################################################################

    ########################################
    # Step 3) create a results dataframes
    ########################################

    # make some common column header names
    common_cols = [
        "method",
        "f_PAU",
        "f_PAFD",
        "f_stdA",
        "algorithm",
        "seed",
        "bcycle",
        "opv_score",
        "mogm_score",
        "mean_gebv",
        "std_gebv",
        "min_gebv",
        "median_gebv",
        "max_gebv"
    ]

    # make population dataframe
    population_df = pandas.DataFrame(
        population_list,
        columns = (
            common_cols + ["gebv"+str(i).zfill(zwidth) for i in range(pop_size)]
        )
    )

    # make selections dataframe
    selections_df = pandas.DataFrame(
        selections_list,
        columns = (
            common_cols + ["gebv"+str(i).zfill(zwidth) for i in range(sel_size)]
        )
    )

    return population_df, selections_df
