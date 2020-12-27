# import 3rd party libraries
import numpy
import time
import pandas
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
from util.human2bytes import human2bytes
import objfn
import algo
import util

def pa(sel_size, geno, coeff, tfreq, tld, gmap, lgroup, cycles, algorithm,
       efreq = None, eld = None, mapfn = None, ldfn = None, mem = None,
       mtype = None, dtype = numpy.dtype("float64"),
       algorithm_varg = None, seed = None, nthreads = 1, zwidth = 3,
       verbose = True):
    """
    Population Architect (PA) Selection

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


    """
    ########################################
    # Step 0) process arguments
    ########################################

    ############################ Parse 'sel_size' #############################
    # process 'sel_size' and make sure it's the right data type
    if isinstance(sel_size, (int, numpy.integer)):
        pass
    elif isinstance(sel_size, (float, numpy.float)):
        sel_size = numpy.int32(numpy.around(sel_size * geno.shape[1]))
    else:
        raise TypeError(
            "'sel_size' must be int, numpy.integer, float, or numpy.float\n"\
            "    type(sel_size) = %s" %
            (type(sel_size),)
        )

    ######################## Parse 'geno' and 'coeff' ########################
    # check for correct genotype shape
    if len(geno.shape) != 3:
        raise TypeError(
            "'geno' must have 3 dimensions: (phase, indiv, locus)\n"\
            "    geno.shape = %s" % (geno.shape,)
        )
    # check for correct coeff shape
    if len(coeff.shape) != 1:
        raise TypeError(
            "'coeff' must have 1 dimension: (locus,)\n"\
            "    coeff.shape = %s" % (coeff.shape,)
        )
    # check if loci lengths are the same
    if geno.shape[2] != coeff.shape[0]:
        raise ValueError(
            "Length of 'coeff' does not match 'geno':\n"\
            "    geno.shape = (%s, %s, locus=%s)\n"\
            "    coeff.shape = (locus=%s,)\n"\
            "'locus' dimensions need to match." %
            (geno.shape[0], geno.shape[1], geno.shape[2], coeff.shape[0])
        )
    # throw an error if the selection size is >= to the number of individuals
    if sel_size >= geno.shape[1]:
        raise ValueError(
            "'sel_size' must be less than the number of indiv:\n"\
            "    sel_size = %s\n"\
            "    geno.shape = (%s, indiv=%s, %s)\n" %
            (sel_size, geno.shape[0], geno.shape[1], geno.shape[2])
        )

    ############################ Calculate 'wcoeff' ############################
    # here we calculate 'wcoeff', which will be input into objfn.pa
    wcoeff = numpy.absolute(coeff)              # take absolute values
    divisor = wcoeff.sum(dtype=numpy.float64)   # take sum of wcoeff array
    divisor *= divisor                          # square the sum
    if mtype in ["tril", "triu"]:               # if we want a triangle matrix
        divisor += (wcoeff*wcoeff).sum(         # add sum of squared elements
            dtype=numpy.float64                 # force double-precision
        )
        divisor /= 2.0                          # divide by two
    wcoeff /= numpy.sqrt(divisor)               # scale so sums equal one

    ######################### Parse 'tfreq' and 'tld' ##########################
    # test if 'tfreq' is an acceptable type
    if not isinstance(tfreq, (float, numpy.float, numpy.float32, numpy.float64,
    numpy.ndarray)):
        raise TypeError(
            "'tfreq' is not a floating point or numpy.ndarray.\n"\
            "    type(tfreq) = %s" % (type(tfreq),)
        )
    if not isinstance(tld, (float, numpy.float, numpy.float32, numpy.float64,
    numpy.ndarray)):
        raise TypeError(
            "'tld' is not a floating point or numpy.ndarray.\n"\
            "    type(tld) = %s" % (type(tfreq),)
        )

    ######################## Parse 'gmap' and 'lgroup' #########################
    # check if 'gmap' and 'lgroup' align
    if len(gmap) != lgroup.sum():
        raise ValueError(
            "The 'gmap' and 'lgroup' do not align:\n"\
            "    len(gmap)    = %s\n"\
            "    lgroup.sum() = %s\n"
            % (len(gmap), lgroup.sum())
        )

    ################# Parse 'algorithm' and 'algorithm_varg' ##################
    # tests for algorithm choice
    if algorithm == 'hc_sa_set':
        # if algorithm_varg is None
        if algorithm_varg is None:
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
        if algorithm_varg is None:
            # make a default dictionary
            algorithm_varg = {
                "states": numpy.tile(numpy.arange(geno.shape[1]), (sel_size,)),
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
                "Expected 'algorithm_varg' dict to have fields:\n"\
                "    states\n"\
                "    dim_sizes"
            )
    elif algorithm == 'icpso':
        # make sure we've received a dictionary
        if not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg),)
            )
        # make sure we have all of our necessary arguments.
        elif any(e not in algorithm_varg.keys() for e in ["n", "states",
                 "dim_sizes", "inertia_wt", "accel_coeff_pbest",
                 "accel_coeff_gbest", "scale_factor", "stpfn"]):
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
                tuple(e not in algorithm_varg.keys() for e in
                        ["n", "states", "dim_sizes", "inertia_wt",
                         "accel_coeff_pbest", "accel_coeff_gbest",
                         "scale_factor", "stpfn"])
            )
    else:
        raise ValueError(
            "Expected 'algorithm' to be one of the following\n"\
            "    'hc_sa_set'\n"\
            "    'hc_sa_state'\n"\
            "    'icpso'"
        )

    ######################### Parse 'efreq' and 'eld' ##########################
    # create 'efreq' if it is None.
    if efreq is None:
        # if 'tfreq' is a numpy.ndarray, generate target matrix
        if isinstance(tfreq, numpy.ndarray):
            efreq = numpy.asarray(
                tfreq < 0.5,        # if tfreq < 0.5 set 1.0; else set 0.0
                dtype=tfreq.dtype   # make identical data types
            )
        # else 'tfreq' should be a floating point since it passed above
        else:
            efreq = 1.0 if tfreq < 0.5 else 0.0
    # create 'eld' if it is None.
    if eld is None:
        # if 'eld' is a numpy.ndarray, generate target matrix
        if isinstance(tld, numpy.ndarray):
            eld = numpy.asarray(
                tld < 0.5,          # if tld < 0.5 set 1.0; else set 0.0
                dtype=tld.dtype     # make identical data types
            )
        # else 'tld' should be a floating point since it passed above
        else:
            eld = 1.0 if tld < 0.5 else 0.0

    ############################### Parse 'mem' ###############################
    # if 'mem' is None, set to compute whole matrix
    if mem is None:
        mem = len(gmap)
    # if 'mem' is an integer type, do nothing
    elif isinstance(mem, (int, numpy.int)):
        pass
    # if 'mem' is a string, parse it and calculate chunk size
    elif isinstance(mem, str):
        mem = int(                              # cast as integer
            math.floor(                         # get floor
                math.sqrt(                      # take square root
                    human2bytes(mem) /          # parse number of required bytes
                    numpy.dtype(dtype).itemsize # divide by dtype size in bytes
                )
            )
        )
    else:
        raise TypeError(
            "Incorrect 'mem' data type.\n"\
            "    Options: None, int, numpy.int, str\n"\
            "    type(mem) = %s" % (type(mem),)
        )
    # if mem is <= 0, raise an error
    if mem <= 0:
        raise ValueError(
            "'mem' must be greater than zero:\n"\
            "    Received: %s" % mem
        )


    ########################################
    # Step 1) establish engine states
    ########################################

    # if no seed is provided, use time as random seed
    if seed is None:
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

    # assemble function variable arguments
    pa_varg = {
        "geno": geno,
        "wcoeff": wcoeff,
        "tfreq": tfreq,
        "tld": tld,
        "gmap": gmap,
        "lgroup": lgroup,
        "cycles": cycles,
        "efreq": efreq,
        "eld": eld,
        "mapfn": mapfn,
        "ldfn": ldfn,
        "mem": mem,
        "mtype": mtype,
        "dtype": dtype
    }

    # define a variable for output by the algorithm
    algo_out = None

    # use the correct algorithm
    if algorithm == "hc_sa_set":
        algo_out = algo.hc_sa_set(
            objfn = objfn.pa,
            objfn_varg = pa_varg,
            k = sel_size,
            **algorithm_varg,        # variable args: contains 'states'
            seed = seed,
            nthreads = nthreads,
            verbose = verbose
        )
    elif algorithm == "hc_sa_state":
        algo_out = algo.hc_sa_state(
            objfn = objfn.pa,
            objfn_varg = pa_varg,
            **algorithm_varg,        # variable args: contains 'states', 'dim_sizes'
            seed = seed,
            nthreads = nthreads,
            verbose = verbose
        )
    elif algorithm == "icpso":
        algo_out = algo.icpso(
            objfn = objfn.pa,
            objfn_varg = pa_varg,
            **algorithm_varg,        # variable args
            seed = seed,
            nthreads = nthreads,
            verbose = verbose
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
def pa_sim(pop_size, sel_size, geno, coeff, tfreq, tld, gmap, lgroup, cycles,
           algorithm, efreq = None, eld = None, mapfn = None, ldfn = None,
           mem = None, mtype = None, dtype = numpy.dtype("float64"),
           algorithm_varg = None, matefn = None, interference = None,
           seed = None, nthreads = 1, zwidth = 3, verbose = True):
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
    ########################################
    # Step 0) process arguments
    ########################################

    ######################## Parse 'geno' and 'coeff' #########################
    # check for correct genotype shape
    if len(geno.shape) != 3:
        raise TypeError(
            "'geno' must have 3 dimensions: (phase, indiv, locus)\n"\
            "    geno.shape = %s" %
            (geno.shape,)
        )
    # check for correct coeff shape
    if len(coeff.shape) != 1:
        raise TypeError(
            "'coeff' must have 1 dimension: (locus,)\n"\
            "    coeff.shape = %s" % (coeff.shape,)
        )
    # check if loci lengths are the same
    if geno.shape[2] != coeff.shape[0]:
        raise ValueError(
            "Length of 'coeff' does not match 'geno':\n"\
            "    geno.shape = (%s, %s, locus=%s)\n"\
            "    coeff.shape = (locus=%s,)\n"\
            "'locus' dimensions need to match." %
            (geno.shape[0], geno.shape[1], geno.shape[2], coeff.shape[0])
        )

    ############################ Calculate 'wcoeff' ############################
    # here we calculate 'wcoeff', which will be input into objfn.pa
    wcoeff = numpy.absolute(coeff)              # take absolute values
    divisor = wcoeff.sum(dtype=numpy.float64)   # take sum of wcoeff array
    divisor *= divisor                          # square the sum
    if mtype in ["tril", "triu"]:               # if we want a triangle matrix
        divisor += (wcoeff*wcoeff).sum(         # add sum of squared elements
            dtype=numpy.float64                 # force double-precision
        )
        divisor /= 2.0                          # divide by two
    wcoeff /= numpy.sqrt(divisor)               # scale so sums equal one

    #################### Parse 'pop_size' and 'sel_size' ######################
    # process 'pop_size' and make sure it's the right data type
    if not numpy.issubdtype(type(pop_size), numpy.integer):
        raise TypeError(
            "Expected 'pop_size' to be a integer type.\n"\
            "    type(pop_size) = %s" % (type(pop_size),)
        )
    # process 'sel_size' and make sure it's the right data type
    if numpy.issubdtype(type(sel_size), numpy.integer):
        pass
    elif numpy.issubdtype(type(sel_size), numpy.floating):
        sel_size = numpy.int32(numpy.around(sel_size * geno.shape[1]))
    else:
        raise TypeError(
            "'sel_size' must be an integer or floating point type.\n"\
            "    type(sel_size) = %s" % (type(sel_size),)
        )
    # test if pop_size and sel_size are sane
    if pop_size < sel_size:
        raise ValueError(
            "'pop_size' must be greater than 'sel_size'.\n"\
            "    pop_size = %s\n"\
            "    sel_size = %s\n" % (pop_size, sel_size)
        )

    # throw an error if the selection size is > to the number of individuals
    if sel_size > geno.shape[1]:
        raise ValueError(
            "'sel_size' must be less than the number of indiv:\n"\
            "    sel_size = %s\n"\
            "    geno.shape = (%s, indiv=%s, %s)\n" %
            (sel_size, geno.shape[0], geno.shape[1], geno.shape[2])
        )

    ######################## Parse 'tfreq' and 'tld' ##########################
    # check if 'tfreq' type is: !ndarray && !floating == !(ndarray || floating)
    if not (isinstance(tfreq, numpy.ndarray) or
            numpy.issubdtype(type(tfreq), numpy.floating)):
        raise TypeError(
            "'tfreq' must be a numpy.ndarray or floating point.\n"\
            "    type(tfreq) = %s" % (type(tfreq),)
        )
    # check if 'tld' type is: !ndarray && !floating == !(ndarray || floating)
    if not (isinstance(tld, numpy.ndarray) or
            numpy.issubdtype(type(tld), numpy.floating)):
        raise TypeError(
            "'tld' must be a numpy.ndarray or floating point.\n"\
            "    type(tld) = %s" % (type(tld),)
        )

    ######################## Parse 'efreq' and 'eld' ##########################
    # check if efreq is None, if it is make a default.
    if efreq is None:
        if isinstance(tfreq, numpy.ndarray):  # optimized for most common
            efreq = numpy.where(tfreq < 0.5, 1.0, 0.0)
        elif numpy.issubdtype(type(tfreq), numpy.floating):
            efreq = 1.0 if tfreq < 0.5 else 0.0
    # check if eld is None, if it is make a default.
    if eld is None:
        if numpy.issubdtype(type(tfreq), numpy.floating): # optimized for most common
            efreq = 1.0 if tfreq < 0.5 else 0.0
        elif isinstance(tld, numpy.ndarray):
            efreq = numpy.where(tfreq < 0.5, 1.0, 0.0)

    ####################### Parse 'gmap' and 'lgroup' #########################
    # make sure gmap is a matrix
    if not isinstance(gmap, numpy.ndarray):
        raise TypeError(
            "'gmap' must be of type numpy.ndarray\n"\
            "    type(gmap) = %s" % (type(gmap),)
        )
    # make sure gmap is a 1D matrix
    if len(gmap.shape) != 1:
        raise ValueError(
            "'gmap' must be of shape (1,):\n"\
            "    gmap.shape = %s" % (gmap.shape,)
        )
    # make sure lgroup is a matrix
    if not isinstance(lgroup, numpy.ndarray):
        raise TypeError(
            "'lgroup' must be of type numpy.ndarray\n"\
            "    type(lgroup) = %s" % (type(lgroup),)
        )
    # make sure lgroup is a 1D matrix
    if len(lgroup.shape) != 1:
        raise ValueError(
            "'lgroup' must be of shape (1,):\n"\
            "    lgroup.shape = %s" % (lgroup.shape,)
        )
    # make sure gmap and lgroup align correctly
    if len(gmap) != lgroup.sum():
        raise ValueError(
            "Length of 'gmap' must equal the sum of elements in 'lgroup'\n"\
            "    len(gmap) = %s\n"\
            "    sum(lgroup) = %s" %
            (len(gmap), lgroup.sum())
        )

    ############################ Parse 'cycles' ##############################
    # test to make sure 'cycles' is good
    if not numpy.issubdtype(type(cycles), numpy.integer):
        raise TypeError(
            "Expected 'cycles' to be an integer type.\n"\
            "    type(cycles) = %s" % (type(cycles),)
        )

    ########################### Parse 'algorithm' #############################
    # process in the case of hc_sa_set algorithm
    if algorithm == 'hc_sa_set':
        # if algorithm_varg is None, make a default dictionary
        if algorithm_varg is None:
            algorithm_varg = {"states": numpy.arange(geno.shape[1])}

        # make sure we've received a dictionary
        elif not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg),)
            )

        # if algorithm_varg does not have a 'states' key, add a default key
        elif "states" not in algorithm_varg.keys():
            algorithm_varg["states"] = numpy.arange(geno.shape[1])
    ########################################################
    # process in the case of hc_sa_state algorithm
    elif algorithm == 'hc_sa_state':
        # if algorithm_varg is None, make a default dictionary
        if algorithm_varg is None:
            algorithm_varg = {
                "states": numpy.tile(numpy.arange(geno.shape[1]), (sel_size,1)),
                "dim_sizes": numpy.repeat(geno.shape[1], sel_size)
            }

        # else if: make sure we've received a dictionary
        elif not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg))
            )

        # else if the 'states' or 'dim_sizes' keys are missing
        elif any(e not in algorithm_varg.keys() for e in ["states","dim_sizes"]):
                raise ValueError(
                    "Expected 'algorithm_varg' dict to have fields:\n"\
                    "    states\n"\
                    "    dim_sizes"
                )
    ########################################################
    # process arguments for the icpso algorithm
    elif algorithm == 'icpso':
        # make sure we've received a dictionary
        if not isinstance(algorithm_varg, dict):
            raise TypeError(
                "Expected 'algorithm_varg' to be a dict\n"\
                "    type(algorithm_varg) = %s" %
                (type(algorithm_varg),)
            )
        # make sure we have all of our necessary arguments.
        elif any(e not in algorithm_varg.keys() for e in ["n", "inertia_wt",
                 "accel_coeff_pbest", "accel_coeff_gbest", "scale_factor",
                 "stpfn"]):
            raise ValueError(
                "Expected 'algorithm_varg' to have all of the following fields:\n"\
                "    Field              Received\n"\
                "    n                  %s\n"\
                "    inertia_wt         %s\n"\
                "    accel_coeff_pbest  %s\n"\
                "    accel_coeff_gbest  %s\n"\
                "    scale_factor       %s\n"\
                "    stpfn              %s\n" %
                tuple(e not in algorithm_varg.keys() for e in ["n",
                      "inertia_wt", "accel_coeff_pbest", "accel_coeff_gbest",
                      "scale_factor", "stpfn"])
            )
        # if the user has not specified the states, dim_sizes, default to all
        if all(e not in algorithm_varg.keys() for e in ["states", "dim_sizes"]):
            algorithm_varg["states"] = numpy.tile(
                numpy.arange(geno.shape[1]),
                (sel_size,)
            )
            algorithm_varg["dim_sizes"] = numpy.repeat(geno.shape[1], sel_size)
        # elif 'states' and 'dim_sizes' have both been provided
        elif all(e in algorithm_varg.keys() for e in ["states", "dim_sizes"]):
            pass
        else:
            raise ValueError(
                "Inconsistent 'state' and 'dim_sizes' in 'algorithm_varg'\n"\
                "    Must provide both or none"
            )
    ########################################################
    # else, the algorithm hasn't been specified
    else:
        raise ValueError(
            "Expected 'algorithm' to be one of the following\n"\
            "    'hc_sa_set'\n"\
            "    'hc_sa_state'\n"\
            "    'icpso'"
        )

    ############################### Parse 'mem' ###############################
    # if 'mem' is None, set to compute whole matrix
    if mem is None:
        mem = len(gmap)
    # if 'mem' is an integer type, do nothing
    elif isinstance(mem, (int, numpy.int)):
        pass
    # if 'mem' is a string, parse it and calculate chunk size
    elif isinstance(mem, str):
        mem = int(                              # cast as integer
            math.floor(                         # get floor
                math.sqrt(                      # take square root
                    human2bytes(mem) /          # parse number of required bytes
                    numpy.dtype(dtype).itemsize # divide by dtype size in bytes
                )
            )
        )
    else:
        raise TypeError(
            "Incorrect 'mem' data type.\n"\
            "    Options: None, int, numpy.int, str\n"\
            "    type(mem) = %s" % (type(mem),)
        )
    # if mem is <= 0, raise an error
    if mem <= 0:
        raise ValueError(
            "'mem' must be greater than zero:\n"\
            "    Received: %s" % mem
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

    # declare a working copy for genotype data, this will be modified by meiosis
    pa_geno = geno.copy()

    for cycle in range(cycles):
        pa_varg = {
            "geno": pa_geno,
            "wcoeff": wcoeff,
            "tfreq": tfreq,
            "tld": tld,
            "gmap": gmap,
            "lgroup": lgroup,
            "cycles": cycles-cycle, # time between deadline and current cycle
            "efreq": efreq,
            "eld": eld,
            "mapfn": mapfn,
            "ldfn": ldfn,
            "mem": mem,
            "mtype": mtype,
            "dtype": dtype
        }

        # declare optimization output variable
        optout = None

        # do the algorithm
        if algorithm == "hc_sa_set":
            optout = algo.hc_sa_set(
                objfn = objfn.pa,
                objfn_varg = pa_varg,
                k = sel_size,
                **algorithm_varg,
                seed = cycle_seed[cycle],
                nthreads = nthreads,
                verbose = verbose
            )
        elif algorithm == "hc_sa_state":
            optout = algo.hc_sa_state(
                objfn = objfn.pa,
                objfn_varg = pa_varg,
                **algorithm_varg,
                seed = cycle_seed[cycle],
                nthreads = nthreads,
                verbose = verbose
            )
        elif algorithm == "icpso":
            optout = algo.icpso(
                objfn = objfn.pa,
                objfn_varg = pa_varg,
                **algorithm_varg,
                seed = cycle_seed[cycle],
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
            (2, sel_size, pa_geno.shape[2]),
            dtype=pa_geno.dtype
        )
        hybrids[0,:,:] = pa_geno[0,females,:]
        hybrids[1,:,:] = pa_geno[1,males,:]

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
        if rands > 0:
            sources[-rands:] = numpy.random.choice(
                numpy.arange(sel_size),
                rands,
                replace = False
            )

        # generate gametes
        gout = util.meiosis(
            hybrids,
            gmap,
            lgroup,
            sources,
            verbose=False
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
        pop_pa_score = objfn.pa(
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
