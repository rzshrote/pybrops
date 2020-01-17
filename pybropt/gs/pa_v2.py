# import 3rd party libraries
import numpy
import time
import pandas
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
from util.human2bytes import human2bytes
from util.meiosis import meiosis
from util.arg_proc import proc_pop_sel_sizes
from util.arg_proc import proc_geno_coeff
from util.arg_proc import proc_wcoeff
from util.arg_proc import proc_tfreq_tld
from util.arg_proc import proc_gmap_lgroup
from util.arg_proc import proc_dcoeff
from util.arg_proc import proc_mem
from util.arg_proc import proc_hc_sa_set_varg
from util.arg_proc import proc_hc_sa_state_varg
from util.arg_proc import proc_icpso_varg
import objfn
import algo

# define a lookup dictionary for mapping function strings to function pointers
PA_ALGORITHM_DICT = {
    'hc_sa_set': algo.hc_sa_set,
    'hc_sa_state': algo.hc_sa_state,
    'icpso': algo.icpso
}


def pa_v2(sel_size, geno, coeff, tfreq, tld, gmap, lgroup, generations, dcoeff,
          algorithm,
          mapfn = 'haldane', ldfn = 'r_sq', mtype = 'tril', mem = None,
          dtype = numpy.dtype("float64"),
          algorithm_varg = None, seed = None, nthreads = 1, zwidth = 3,
          verbose = True, verbose_algo = False):
    """
    Population Architect (PA) selection

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
    wcoeff = proc_wcoeff(None, coeff)               # calculate 'wcoeff'
    tfreq, tld = proc_tfreq_tld(tfreq, tld)         # process 'tfreq' && 'tld'
    gmap, lgroup = proc_gmap_lgroup(gmap, lgroup)   # process 'gmap' && 'lgroup'
    dcoeff = proc_dcoeff(dcoeff, 3)                 # process 'dcoeff'
    mem = proc_mem(mem, gmap)                       # process 'mem'

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
        print("PA_v2 engine state:")
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ########################################
    # Step 2) compute results
    ########################################

    # assemble function variable arguments
    pa_v2_varg = {
        "geno": geno,
        "wcoeff": wcoeff,
        "tfreq": tfreq,
        "tld": tld,
        "gmap": gmap,
        "lgroup": lgroup,
        "generations": generations,
        "dcoeff": dcoeff,
        "mapfn": mapfn,
        "ldfn": ldfn,
        "mtype": mtype,
        "mem": mem,
        "dtype": dtype
    }

    # lookup the algorithm in the dict, and execute the function with arguments
    algo_out = PA_ALGORITHM_DICT[algorithm](
        objfn = objfn.pa_v2,
        objfn_varg = pa_v2_varg,
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
            data = [["opv",algorithm,seed]],
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
            data = [["opv",algorithm,seed]
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
def pa_v2_sim(pop_size, sel_size, geno, coeff, tfreq, tld, gmap, lgroup,
              dcoeff_t, bcycles, algorithm,
              mapfn = 'haldane', ldfn = 'r_sq', mtype = 'tril', mem = None,
              dtype = numpy.dtype("float64"),
              cgenerations = 1, matefn = None, algorithm_varg = None,
              interference = None,
              seed = None, nthreads = 1, zwidth = 3, verbose = True,
              averbose = False):
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
    wcoeff = proc_wcoeff(None, coeff)               # process 'wcoeff'
    tfreq, tld = proc_tfreq_tld(tfreq, tld)         # process 'tfreq' && 'tld'
    gmap, lgroup = proc_gmap_lgroup(gmap, lgroup)   # process 'gmap' && 'lgroup'
    dcoeff_t = proc_dcoeff_t(dcoeff_t, bcycles, 3)  # process 'dcoeff_t'
    bcycles = proc_bcycles(bcycles)                 # process 'bcycles'
    mem = proc_mem(mem, gmap)                       # process 'mem'

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
        print("PA simulation engine state:")
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    # generate random seeds for cycle, shuffle, meiosis, for each cycle
    seeds = numpy.random.randint(
        low = numpy.iinfo(numpy.uint32).min,
        high = numpy.iinfo(numpy.uint32).max,
        size = (3,bcycles),
        dtype = numpy.uint32
    )

    # declare lists to store results into
    selections_list = list()
    population_list = list()

    # declare a working copy for genotype data, this will be modified by meiosis
    pa_geno = geno.copy()

    for bcycle in range(bcycles):
        pa_varg = {
            "geno": pa_geno,
            "wcoeff": wcoeff,
            "tfreq": tfreq,
            "tld": tld,
            "gmap": gmap,
            "lgroup": lgroup,
            # generations between deadline and current breeding cycle
            "generations": (bcycles-bcycle) * cgenerations,
            "dcoeff": dcoeff_t[bcycle,:], # get row for current breeding cycle
            "mapfn": mapfn,
            "ldfn": ldfn,
            "mtype": mtype,
            "mem": mem,
            "dtype": dtype
        }

        # use algorithm selected, get output
        optout = PA_ALGORITHM_DICT[algorithm](
            objfn = objfn.pa,
            objfn_varg = pa_varg,
            **algorithm_varg,
            seed = seeds[0,bcycle],
            nthreads = nthreads,
            verbose = verbose
        )

        selection_indices = None
        selection_score = None

        # process output
        if algorithm == "hc_sa_set" or algorithm == "hc_sa_state":
            selection_indices = optout["X_gbest_pos"]   # get selection_indices
            selection_score = optout["X_gbest_scr"]     # get selection_score
        elif algorithm == "icpso":
            selection_indices = optout["X_gbest_smpl"]  # get selection_indices
            selection_score = optout["X_gbest_scr"]     # get selection_score

        # calculate GEBVs, put into list
        gebvs = numpy.dot(
            pa_geno[:,selection_indices,:],
            coeff
        ).sum(0).tolist()

        # build method
        method = "pa"+"_"+str(mapfn)+"_"+str(ldfn)

        # store everything into selections_list
        selections_list.append(
            [method, algorithm, seed, cycle+1, selection_score] + gebvs
        )

        # seed the rng for shuffling
        numpy.random.seed(seeds[1,bcycle])

        # TODO: implement other mating methods

        # make parents
        females = selection_indices.copy()
        males = selection_indices.copy()

        # shuffle parents
        numpy.random.shuffle(females)
        numpy.random.shuffle(males)

        # make hybrids
        hybrids = numpy.empty(
            (2, sel_size, pa_geno.shape[2]),
            dtype=pa_geno.dtype
        )
        hybrids[0,:,:] = pa_geno[0,females,:]
        hybrids[1,:,:] = pa_geno[1,males,:]

        # make gametes production assignments

        # calculate how much to repeat, how much to fill
        reps = pop_size // sel_size
        rands = pop_size % sel_size

        # make empty gamete sources array
        sources = numpy.empty(pop_size, dtype=numpy.int)
        sources[0:reps*sel_size] = numpy.repeat(
            numpy.arange(sel_size),
            reps
        )
        if rands > 0:
            sources[-rands:] = numpy.random.choice(
                numpy.arange(sel_size),
                rands,
                replace = False
            )

        # generate gametes
        gout = meiosis(
            hybrids,
            gmap,
            lgroup,
            sources,
            interference = interference,
            seed = seeds[2,bcycle],
            verbose = False
        )

        # double these gamets to make DH; replace prev. data in 'pa_geno'
        pa_geno = numpy.array([gout, gout])

        # calculate population GEBVs, put into list
        gebvs = numpy.dot(
            pa_geno,
            coeff
        ).sum(0).tolist()

        # score the population using objfn.opv
        pop_opv_score = objfn.opv(
            slice(None),
            pa_geno * coeff
        )

        # score the population using objfn.pa
        pop_pa_score = objfn.pa_v2(
            slice(None),
            **pa_varg
        )

        # append population stats to list
        population_list.append(
            [method, algorithm, seed, cycle+1, pop_opv_score, pop_pa_score] + gebvs
        )

        # print progress
        if verbose:
            print("Cycle:", cycle+1)

    ########################################
    # Step 3) create a results dataframes
    ########################################

    # make population dataframe
    population_df = pandas.DataFrame(
        population_list,
        columns = (["method", "algorithm", "seed", "cycle", "opv_score", "pa_score"] +
            ["gebv"+str(i).zfill(zwidth) for i in range(pop_size)])
    )

    # make selections dataframe
    selections_df = pandas.DataFrame(
        selections_list,
        columns = (["method", "algorithm", "seed", "cycle", "score"] +
            ["gebv"+str(i).zfill(zwidth) for i in range(sel_size)])
    )

    return population_df, selections_df
