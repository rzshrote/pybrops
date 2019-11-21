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

def opv(geno,
        coeff,
        sel_size,
        algorithm,
        algorithm_varg = None,
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

    # process 'sel_size' and make sure it's the right data type
    if isinstance(sel_size, (int, numpy.integer)):
        sel_size = sel_size if sel_size < geno.shape[1] else geno.shape[1]
    elif isinstance(sel_size, (float, numpy.float)):
        sel_size = numpy.int32(numpy.around(sel_size))
        sel_size = sel_size if sel_size < geno.shape[1] else geno.shape[1]
    else:
        raise TypeError(
            "'sel_size' must be int, numpy.integer, float, or numpy.float\n"\
            "    type(sel_size) = %s" %
            (type(sel_size),)
        )

    # tests for algorithm choice
    if algorithm == 'hc_sa_set':
        # if algorithm_varg is None
        if isinstance(algorithm_varg, type(None)):
            # make a default dictionary
            algorithm_varg = {
                "states": numpy.arange(geno.shape[1])
            }
        # make sure we've received a dictionary
        elif not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg),)
            )
        # if algorithm_varg does not have a 'states' key
        elif "states" not in algorithm_varg.keys():
            # add a default 'states' key
            algorithm_varg["states"] = numpy.arange(geno.shape[1])
    elif algorithm == 'hc_sa_state':
        if isinstance(algorithm_varg, type(None)):
            # make a default dictionary
            algorithm_varg = {
                "states": numpy.tile(numpy.arange(geno.shape[1]), (sel_size,1)),
                "dim_sizes": numpy.repeat(geno.shape[1], sel_size)
            }
        # make sure we've received a dictionary
        elif not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg))
            )
        elif "states" or "dim_sizes" not in algorithm_varg.keys():
            raise ValueError(
                "Expected 'algorithm_varg' dict to have fields "\
                "'states' and 'dim_sizes'"
            )
    elif algorithm == 'icpso':
        # make sure we've received a dictionary
        if not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg))
            )
        # make sure we have all of our necessary arguments.
        elif ("n"                 not in algorithm_varg.keys() or
              "states"            not in algorithm_varg.keys() or
              "dim_sizes"         not in algorithm_varg.keys() or
              "inertia_wt"        not in algorithm_varg.keys() or
              "accel_coeff_pbest" not in algorithm_varg.keys() or
              "accel_coeff_gbest" not in algorithm_varg.keys() or
              "scale_factor"      not in algorithm_varg.keys() or
              "stpfn"             not in algorithm_varg.keys()):
            raise ValueError(
                "Expected 'algorithm_varg' to have all of the following fields:\n"\
                "    Field              Received\n"\
                "    n                  %s\n"\
                "    states             %s\n"\
                "    dim_sizes          %s\n"\
                "    inertia_wt         %s\n"\
                "    accel_coeff_pbest  %s\n"\
                "    accel_coeff_gbest  %s\n"\
                "    scale_factor       %s\n"\
                "    stpfn              %s\n" %
                ("n" in algorithm_varg.keys(),
                 "states" in algorithm_varg.keys(),
                 "dim_sizes" in algorithm_varg.keys(),
                 "inertia_wt" in algorithm_varg.keys(),
                 "accel_coeff_pbest" in algorithm_varg.keys(),
                 "accel_coeff_gbest" in algorithm_varg.keys(),
                 "scale_factor" in algorithm_varg.keys(),
                 "stpfn" in algorithm_varg.keys())
            )
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

    # use the correct algorithm
    if algorithm == "hc_sa_set":
        algo_out = algo.hc_sa_set(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            k = sel_size,
            **algorithm_varg,        # variable args: contains 'states'
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    elif algorithm == "hc_sa_state":
        algo_out = algo.hc_sa_state(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            **algorithm_varg,        # variable args: contains 'states', 'dim_sizes'
            seed = seed,
            nthreads = nthreads,
            verbose = False
        )
    elif algorithm == "icpso":
        algo_out = algo.icpso(
            objfn = objfn.opv,
            objfn_varg = opv_varg,
            **algorithm_varg,        # variable args
            seed = seed,
            nthreads = nthreads,
            verbose = False
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
            data = ["opv",algorithm,seed],
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
            data = algo_out["X_pbest_smpl"].T,
            columns = ["x"+str(i).zfill(zwidth)
                       for i in range(algo_out["X_gbest_pos"].shape[1])]
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
                   for i in range(history_df["x_pos"][0])]
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
def opv_sim(geno,
            coeff,
            pop_size,
            sel_size,
            d,
            lgroup_size,
            cycles,
            algorithm,
            algorithm_varg = None,
            interference = None,
            seed = None,
            nthreads = 1,
            zwidth = 3,
            verbose = True):
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
    d : numpy.ndarray
        Genetic map in Morgan units.
    lgroup_size : numpy.ndarray
        Array of linkage group sizes. The sum of the elements should equal the
        length of 'd'.
    cycles : int
        Number of breeding cycles to simulate.
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

    # process 'pop_size' and make sure it's the right data type
    if not isinstance(pop_size, (int, numpy.integer)):
        raise TypeError(
            "Expected 'pop_size' to be an int or numpy.integer\n"\
            "    type(pop_size) = %s" %
            (type(pop_size),)
        )

    # process 'sel_size' and make sure it's the right data type
    if isinstance(sel_size, (int, numpy.integer)):
        sel_size = sel_size if sel_size < geno.shape[1] else geno.shape[1]
    elif isinstance(sel_size, (float, numpy.float)):
        sel_size = numpy.int32(numpy.around(sel_size))
        sel_size = sel_size if sel_size < geno.shape[1] else geno.shape[1]
    else:
        raise TypeError(
            "'sel_size' must be int, numpy.integer, float, or numpy.float\n"\
            "    type(sel_size) = %s" %
            (type(sel_size),)
        )

    # process d and lgroup_size
    if not isinstance(d, numpy.ndarray):
        raise TypeError(
            "'d' must be of type numpy.ndarray\n"\
            "    type(d) = %s" %
            (type(d),)
        )
    if not isinstance(lgroup_size, numpy.ndarray):
        raise TypeError(
            "'lgroup_size' must be of type numpy.ndarray\n"\
            "    type(lgroup_size) = %s" %
            (type(lgroup_size),)
        )
    if len(d) != lgroup_size.sum():
        raise ValueError(
            "Length of 'd' must equal the sum of elements in 'lgroup_size'\n"\
            "    len(d) = %s\n"\
            "    sum(lgroup_size) = %s" %
            (len(d), lgroup_size.sum())
        )

    # test to make sure 'cycles' is good
    if not isinstance(cycles, (int, numpy.integer)):
        raise TypeError(
            "Expected 'cycles' to be an int or numpy.integer\n"\
            "    type(cycles) = %s" %
            (type(cycles),)
        )

    # tests for algorithm choice
    if algorithm == 'hc_sa_set':
        # if algorithm_varg is None
        if isinstance(algorithm_varg, type(None)):
            # make a default dictionary
            algorithm_varg = {
                "states": numpy.arange(geno.shape[1])
            }
        # make sure we've received a dictionary
        elif not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg),)
            )
        # if algorithm_varg does not have a 'states' key
        elif "states" not in algorithm_varg.keys():
            # add a default 'states' key
            algorithm_varg["states"] = numpy.arange(geno.shape[1])
    elif algorithm == 'hc_sa_state':
        if isinstance(algorithm_varg, type(None)):
            # make a default dictionary
            algorithm_varg = {
                "states": numpy.tile(numpy.arange(geno.shape[1]), (sel_size,1)),
                "dim_sizes": numpy.repeat(geno.shape[1], sel_size)
            }
        # make sure we've received a dictionary
        elif not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg))
            )
        elif "states" or "dim_sizes" not in algorithm_varg.keys():
            raise ValueError(
                "Expected 'algorithm_varg' dict to have fields "\
                "'states' and 'dim_sizes'"
            )
    elif algorithm == 'icpso':
        # make sure we've received a dictionary
        if not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg))
            )
        # make sure we have all of our necessary arguments.
        elif ("n"                 not in algorithm_varg.keys() or
              "inertia_wt"        not in algorithm_varg.keys() or
              "accel_coeff_pbest" not in algorithm_varg.keys() or
              "accel_coeff_gbest" not in algorithm_varg.keys() or
              "scale_factor"      not in algorithm_varg.keys() or
              "stpfn"             not in algorithm_varg.keys()):
            raise ValueError(
                "Expected 'algorithm_varg' to have all of the following fields:\n"\
                "    Field              Received\n"\
                "    n                  %s\n"\
                "    inertia_wt         %s\n"\
                "    accel_coeff_pbest  %s\n"\
                "    accel_coeff_gbest  %s\n"\
                "    scale_factor       %s\n"\
                "    stpfn              %s\n" %
                ("n" in algorithm_varg.keys(),
                 "inertia_wt" in algorithm_varg.keys(),
                 "accel_coeff_pbest" in algorithm_varg.keys(),
                 "accel_coeff_gbest" in algorithm_varg.keys(),
                 "scale_factor" in algorithm_varg.keys(),
                 "stpfn" in algorithm_varg.keys())
            )
        # if the user has not specified the states, dim_sizes, default to all
        if ("states"    not in algorithm_varg.keys() and
            "dim_sizes" not in algorithm_varg.keys()):
            algorithm_varg["states"] = numpy.tile(
                numpy.arange(geno.shape[1]),
                (sel_size,1)
            )
            algorithm_varg["dim_sizes"] = numpy.repeat(geno.shape[1], sel_size)
        # elif 'states' and 'dim_sizes' have both been provided
        elif ("states"    in algorithm_varg.keys() and
              "dim_sizes" in algorithm_varg.keys()):
            pass
        else:
            raise ValueError(
                "Inconsistent 'state' and 'dim_sizes' in 'algorithm_varg'\n"\
                "    Must provide both or none"
            )
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
        print("OPV engine state:")
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    # generate random seeds for algorithm seeding at each cycle
    cycle_seed = numpy.random.randint(
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

    # generate random seeds for shuffle seeding at each cycle
    shuffle_seed = numpy.random.randint(
        low = numpy.iinfo(numpy.uint32).min,
        high = numpy.iinfo(numpy.uint32).max,
        size = cycles,
        dtype = numpy.uint32
    )

    # declare lists to store results into
    selections_list = list()
    population_list = list()

    # declare a working copy for genotype data
    opv_geno = geno.copy()

    for cycle in range(cycles):
        # declare opv_varg
        opv_varg = {"hcoeff": opv_geno * coeff}

        # declare optimization output variable
        optout = None

        # do the algorithm
        if algorithm == "hc_sa_set":
            optout = algo.hc_sa_set(
                objfn = objfn.opv,
                objfn_varg = opv_varg,
                k = sel_size,
                **algorithm_varg,
                seed = cycle_seed[cycle],
                nthreads = nthreads,
                verbose = False
            )
        elif algorithm == "hc_sa_state":
            optout = algo.hc_sa_state(
                objfn = objfn.opv,
                objfn_varg = opv_varg,
                **algorithm_varg,
                seed = cycle_seed[cycle],
                nthreads = nthreads,
                verbose = False
            )
        elif algorithm == "icpso":
            optout = algo.icpso(
                objfn = objfn.opv,
                objfn_varg = opv_varg,
                **algorithm_varg,
                seed = cycle_seed[cycle],
                nthreads = nthreads,
                verbose = False
            )

        selection_indices = None
        selection_score = None

        # process output
        if algorithm == "hc_sa_set" or algorithm == "hc_sa_state":
            # get selection_indices
            selection_indices = optout["X_gbest_pos"]
            # get selection_score
            selection_score = optout["X_gbest_scr"]
        elif algorithm == "icpso":
            # get selection_indices
            selection_indices = optout["X_gbest_smpl"]
            # get selection_score
            selection_score = optout["X_gbest_scr"]

        # calculate GEBVs, put into list
        gebvs = numpy.dot(
            opv_geno[:,selection_indices,:],
            coeff
        ).sum(0).tolist()

        # store everything into selections_list
        selections_list.append(
            ["opv", algorithm, seed, cycle+1, selection_score] + gebvs
        )

        # seed the rng for shuffling
        numpy.random.seed(shuffle_seed[cycle])

        # TODO: implement other mating methods

        # make parents
        females = selection_indices.copy()
        males = selection_indices.copy()

        # shuffle parents
        numpy.random.shuffle(females)
        numpy.random.shuffle(males)

        # make hybrids
        hybrids = numpy.empty(
            (2, sel_size, opv_geno.shape[2]),
            dtype=opv_geno.dtype
        )
        hybrids[0,:,:] = opv_geno[0,females,:]
        hybrids[1,:,:] = opv_geno[1,males,:]

        # make gametes production assignments

        # calculate how much to repeat, how much to
        reps = pop_size // sel_size
        rands = pop_size % sel_size

        # make empty gamete sources array
        sources = numpy.empty(pop_size, dtype=numpy.int)
        sources[0:reps*sel_size] = numpy.repeat(
            numpy.arange(sel_size),
            reps
        )
        sources[-rands:] = numpy.random.choice(
            numpy.arange(sel_size),
            rands,
            replace = False
        )

        # generate gametes
        gout = util.meiosis(
            hybrids,
            d,
            lgroup_size,
            sources,
            verbose=False
        )

        # double these gamets to make DH; replace prev. data in 'opv_geno'
        opv_geno = numpy.array([gout, gout])

        # calculate population GEBVs, put into list
        gebvs = numpy.dot(
            opv_geno,
            coeff
        ).sum(0).tolist()

        # score the population using objfn.opv
        population_score = objfn.opv(
            numpy.arange(pop_size),
            opv_geno * hcoeff
        )

        # append population stats to list
        population_list.append(
            ["opv", algorithm, seed, cycle+1, population_score] + gebvs
        )

        # recalculate haplotype coefficients
        varg = {"hcoeff": phases * effects}

        # print progress
        if verbose:
            print("Cycle:", cycle+1)

    ########################################
    # Step 3) create a results dataframes
    ########################################

    # make population dataframe
    population_df = pandas.DataFrame(
        population_list,
        columns = (["method", "algorithm", "seed", "cycle", "score"] +
            ["gebv"+str(i).zfill(zwidth) for i in range(pop_size)])
    )

    # make selections dataframe
    selections_df = pandas.DataFrame(
        selections_list,
        columns = (["method", "algorithm", "seed", "cycle", "score"] +
            ["gebv"+str(i).zfill(zwidth) for i in range(sel_size)])
    )

    return population_df, selections_df
