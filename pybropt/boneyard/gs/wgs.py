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

def wgs(geno,
        coeff,
        sel_size,
        wgtfn = None,
        algorithm = None,
        zwidth = 3,
        verbose = True):
    """
    Weighted Genomic Selection (CGS)

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
    algorithm : None, {'quicksort', 'mergesort', 'heapsort', 'stable'}
        Optional specification for algorithm to use for sorting GEBVs. Default
        is 'quicksort'. Only worry about this for an extremely large number of
        individuals.

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
    if algorithm not in [None, 'quicksort', 'mergesort', 'heapsort', 'stable']:
        raise ValueError(
            "Expected 'algorithm' to be one of the following\n"\
            "    None\n"\
            "    'quicksort'\n"\
            "    'mergesort'\n"\
            "    'heapsort'\n"\
            "    'stable'"
        )

    ########################################
    # Step 1) calculate wGEBVs
    ########################################

    # calculate weights
    weight = None
    if wgtfn is None:
        pfreq = geno.sum((0,1)) / geno.shape[1]
        weight = numpy.where(pfreq == 0.0, 1.0, 1.0/numpy.sqrt(pfreq))
    else:
        weight = wgtfn(rslice, geno)

    # calculate wGEBVs
    wgebv = objfn.wgs(
        rslice = slice(None), # select everything
        geno = geno,
        coeff = coeff,
        weight = weight
    )

    ########################################
    # Step 2) calculate indices for sorted array
    ########################################
    wgebv_argsort = numpy.argsort(wgebv, kind=algorithm)

    ########################################
    # Step 3) create a results dataframe
    ########################################

    selection_indices = wgebv_argsort[-sel_size:]

    # make results dataframe
    results_df = pandas.DataFrame(
        data = [
            ["wgs", algorithm, "NA", wgebv[selection_indices].sum()] +
            selection_indices.tolist()
        ],
        columns = (
            ["method", "algorithm", "seed", "score"] +
            ["x" + str(i).zfill(zwidth) for i in range(sel_size)]
        )
    )

    # make history dataframe (this algorithm does not generate a history)
    history_df = None

    # return the 'sel_size' highest scoring indices
    return results_df, history_df

def wgs_sim(geno,
            coeff,
            pop_size,
            sel_size,
            cycles,
            d,
            lgroup_size,
            wgtfn = None,
            algorithm = None,
            interference = None,
            seed = None,
            nthreads = 1,
            zwidth = 3,
            verbose = True):
    """
    Conventional Genomic Selection (CGS)

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
    lgroup_size : numpy.ndarray
        Array of linkage group sizes. The sum of the elements should equal the
        length of 'd'.
    algorithm : None, {'quicksort', 'mergesort', 'heapsort', 'stable'}
        Optional specification for algorithm to use for sorting GEBVs. Default
        is 'quicksort'. Only worry about this for an extremely large number of
        individuals.
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

    # make sure d is a numpy.ndarray
    if not isinstance(d, numpy.ndarray):
        raise TypeError(
            "'d' must be of type numpy.ndarray\n"\
            "    type(d) = %s" %
            (type(d),)
        )

    # make sure lgroup_size is a numpy.ndarray
    if not isinstance(lgroup_size, numpy.ndarray):
        raise TypeError(
            "'lgroup_size' must be of type numpy.ndarray\n"\
            "    type(lgroup_size) = %s" %
            (type(lgroup_size),)
        )

    # make sure that the linkage map and the bins line up
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
    if algorithm not in [None, 'quicksort', 'mergesort', 'heapsort', 'stable']:
        raise ValueError(
            "Expected 'algorithm' to be one of the following\n"\
            "    None\n"\
            "    'quicksort'\n"\
            "    'mergesort'\n"\
            "    'heapsort'\n"\
            "    'stable'"
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
        print("CGS simulation engine state:")
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
    wgs_geno = geno.copy()

    # declare a working copy for marker weights
    wgs_weight = None

    # simulate the breeding cycles
    for cycle in range(cycles):
        # calculate weights
        if wgtfn is None:
            pfreq = geno.sum((0,1)) / geno.shape[1]
            wgs_weight = numpy.where(pfreq == 0.0, 1.0, 1.0/numpy.sqrt(pfreq))
        else:
            wgs_weight = wgtfn(rslice, geno)

        # calculate wGEBVs
        wgebv = objfn.wgs(
            rslice = slice(None), # select everything
            geno = wgs_geno,
            coeff = coeff,
            weight = wgs_weight
        )

        # sort based on GEBV
        wgebv_argsort = numpy.argsort(
            wgebv,
            kind = algorithm
        )

        # get the top 'sel_size' individuals
        selection_indices = wgebv_argsort[-sel_size:]

        # get the sum of scores to rate the population
        selection_score = wgebv[selection_indices].sum()

        # calculate GEBVs, put into list
        gebvs = numpy.dot(
            wgs_geno[:,selection_indices,:],
            coeff
        ).sum(0).tolist()

        # store everything into selections_list
        selections_list.append(
            ["cgs", algorithm, seed, cycle+1, selection_score] + gebvs
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
            (2, sel_size, wgs_geno.shape[2]),
            dtype=wgs_geno.dtype
        )
        hybrids[0,:,:] = wgs_geno[0,females,:]
        hybrids[1,:,:] = wgs_geno[1,males,:]

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
            d,
            lgroup_size,
            sources,
            verbose=False
        )

        # double these gamets to make DH; replace prev. data in 'wgs_geno'
        wgs_geno = numpy.array([gout, gout])

        # calculate population GEBVs, put into list
        gebvs = numpy.dot(
            wgs_geno,
            coeff
        ).sum(0).tolist()

        # score the population using objfn.opv
        population_score = objfn.wgs(
            rslice = slice(None),
            geno = wgs_geno,
            coeff = coeff,
            weight = wgs_weight
        ).sum()

        # append population stats to list
        population_list.append(
            ["cgs", algorithm, seed, cycle+1, population_score] + gebvs
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
