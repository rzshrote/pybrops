import numpy
import time

def bqpso(objfn,
          objfn_varg,
          n,
          d,
          alpha,
          cofreq,
          stpfn,
          initprob = 0.5,
          seed = None,
          nthreads = 1,
          verbose = True):
    """
    Perform binary quantum particle swarm optimization (BQPSO)

    Parameters
    ==========
    objfn : function
        Objective function to optimize.
    objfn_varg : dictionary
        Objective function variable arguments.
    n : int
        Number of particles in the swarm.
    d : int
        Number of dimensions in the swarm.
    alpha : double
        Expansion-contraction coefficient.
    cofreq : double
        Crossover frequency.
    stpfn : function
        A stop function. Should accept iteration number, ... some other things.
        Returns False to stop the iterations.
    initprob : double
        Probability of '1' bit being placed in the initial swarm positions.
    seed : int
        Random number seed.
    nthreads : int
        Number of threads.
    verbose : boolean
        Show verbose messages.
    """

    ###################################
    #### Step 0) Process arguments ####
    ###################################

    # if no seed is provided, use time as random seed
    if seed == None:
        seed = numpy.uint32(time.time())

    # seed random number generator
    numpy.random.seed(seed)

    # print engine state before beginning.
    if verbose:
        print("Binary Quantum PSO engine state:")
        print("Objective function", "==================", sep='\n')
        print("Objective function =", objfn.__name__)
        print("Objective function variable arguments =", objfn_varg)
        print("Swarm configuration", "===================", sep='\n')
        print("Number of particles =", n)
        print("Number of dimensions =", d)
        print("Expansion-contraction coefficient =", alpha)
        print("Crossover frequency =", cofreq)
        print("Stopping function =", stpfn.__name__)
        print("Miscellaneous", "=============", sep='\n')
        print("Initialization probability =", initprob)
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ################################
    #### Step 1) Initialization ####
    ################################

    ######################################################
    # Step 1a) allocate memory, initialize select arrays #
    ######################################################

    # create list to store history of the swarm.
    history = list()

    # create an iteration counter
    iter = 0

    # create X_pos matrix
    X_pos = numpy.random.binomial(1, initprob, (n,d)).astype(numpy.uint8)

    # create X_pbest_pos matrix
    X_pbest_pos = X_pos.copy()

    # create X_pbest_scr matrix, fill with lowest value possible
    # need to fill to avoid memory bugs
    X_pbest_scr = numpy.empty(n, dtype=numpy.float64)
    X_pbest_scr.fill(numpy.finfo(numpy.float64).min)

    # create X_mbest_pos matrix
    X_mbest_pos = numpy.empty(d, dtype=numpy.uint8)

    # create X_gbest_pos matrix
    X_gbest_pos = numpy.empty(d, dtype=numpy.uint8)

    # create localAttractors matrix
    localAttractors = numpy.empty((n,d), dtype=numpy.uint8)

    # create quantumWell_probs matrix
    quantumWell_probs = numpy.empty((n,d), dtype=numpy.float64)

    while stpfn(iter):
        ########################################
        # step 1) calculate mean best positions
        ########################################

        # calculate mean positions (1D array)
        cmean = numpy.mean(X_pbest_pos)

        # determine random flips for all positions
        rflip = numpy.random.binomial(
            1,
            0.5,
            cmean.shape
        ).astype(numpy.uint8)

        # if the mean == 0.5, then use a random flip, else set the bit to most
        # frequent bit; this should be an array of uint8
        X_mbest_pos = numpy.where(
            cmean == 0.5,
            rflip,
            cmean > 0.5
        )

        ########################################
        # step 2) calculate score positions
        ########################################

        # score our positions, update the personal positions
        if nthreads > 1:
            print("Greater than 1 threads not supported.")
        else:
            # for each particle
            for i in range(n):
                # score the position
                tmpscore = objfn(X_pos[i,:], **objfn_varg)
                # if the particle scores better
                if tmpscore > X_pbest_scr:
                    X_pbest_scr[i] = tmpscore       # copy the score
                    X_pbest_pos[i,:] = X_pos[i,:]   # copy the position

        # calculate the global best index
        X_gbest_index = numpy.argmax(X_pbest_scr)

        # set the global best position
        X_gbest_pos = X_pbest_pos[X_gbest_index,:]

        # calculate quantum well probabilities
        quantumWell_probs = \
            alpha * \
            (X_pos ^ X_mbest_pos) * \
            numpy.log(1.0 / (1.0 - numpy.random.uniform(0.0, 1.0, (n,d))))

        # transform X_pos
        X_pos = numpy.where(
            numpy.random.uniform(0.0, 1.0, (n,d)) < quantumWell_probs,
            X_pos ^ 1,
            X_pos
        )

    return X_gbest_pos
