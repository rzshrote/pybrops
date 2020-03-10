import numpy
import time

def ga_fps(objfn,
           objfn_varg,
           n,
           states,
           dim_sizes,
           p_rec,
           p_mut,
           stpfn,
           fpsfn = None,
           seed = None,
           nthreads = 1,
           verbose = True):
    """
    Perform Integer/Categorical Particle Swarm Optimization.

    Parameters
    ----------
    objfn : function
        Objective function to optimize.
    objfn_varg : dictionary
        Objective function variable arguments.
    n : int
        Number of particles in the swarm
    states : numpy.ndarray
        1D array of possible discrete positional states for the objective
        function. The boundaries between different dimensions are determined by
        the 'dim_sizes' array.
        Example:
            states = numpy.array([0,1,2,3,1,2,3,0,2,3])
            Dimension blocks are: [0,1,2,3], [1,2,3], [0,2,3]
    dim_sizes : numpy.ndarray
        1D array of integer dimension sizes for the states vector. This is used
        to determine boundaries for possible states in the 'states' array. The
        sum of elements in this array must equal the length of 'states'.
        The dtype of dim_sizes should be numpy.uint64.
        NOTE: Assumes all state sizes are greater than 1.
        Example:
            states = numpy.array([0,1,2,3,1,2,3,0,2,3])
            dim_sizes = numpy.array([4,3,3])
    p_rec : numpy.float64
        Probability of recombination between two individuals.
    p_mut : numpy.float64
        Probability of locus mutation.
    stpfn : function
        A stop function. Should accept iteration number, ... some other things.
        Returns False to stop the iterations.
    fpsfn : function, None
        A function for generating fitness proportion fitness values (FPS). This
        function accepts a single argument. The argument is an array of scores
        from the population as scored by 'objfn'. The function should return an
        array of values. The array should sum to unity.
        If None, then a default FPS is applied. The default FPS scales
        population scores such that they all sum to unity. NOTE: This assumes
        fitness values are never negative! If negative values are possible, set
        a fpsfn to be something that accounts for this!
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
        print("Genetic algorithm with FPS engine state:")
        print("Objective function", "==================", sep='\n')
        print("Objective function =", objfn.__name__)
        print("Objective function variable arguments =", objfn_varg)
        print("Population configuration", "===================", sep='\n')
        print("Population size =", n)
        print("Chromosome length =", len(dim_sizes))
        print("Possible chromosome alleles =", states)
        print("Number of alleles per chromosome locus =", dim_sizes)
        print("Stopping function =", stpfn.__name__)
        print("Fitness proportion selection function = ", fpsfn.__name__)
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ################################
    #### Step 1) Initialization ####
    ################################

    ######################################################
    # Step 1a) allocate memory, initialize select arrays #
    ######################################################

    # create an array of dimension state stop indices
    # 1D array for dimension stop indices
    dsp = numpy.cumsum(dim_sizes)

    # create an array of dimension state start indices
    # 1D arrray for dimension start indices
    # take 'dsp' and subtract number at first index
    dst = dsp - dim_sizes[0]

    # Matrix of chromosomes. Each row of the matrix contains all states, but
    # only the indices indicated by 'dst' are used.
    # Example:
    # Chromosome:    {1} 2 3 | {2} 4 3 1 | {2} 1 3 5 4
    # This chromosome represents the solution {1,2,2}
    # Tile the input states to get 'n' individuals
    X_pos = numpy.tile(states, (n,1))

    # randomly shuffle the chromosome states to get a random population.
    for i in range(n):
        for st,sp in zip(dst,dsp):
            numpy.random.shuffle(X_pos[i,st:sp])

    # allocate a vector to contain our individual scores
    X_scr = numpy.empty(n, dtype=numpy.float64)

    # score individuals
    if nthreads > 1:
        print("Greater than 1 thread not supported yet.")
        return
    else:
        # for each individual
        for i in range(n):
            # score the individual and store in score vector
            X_scr[i] = objfn(X_pos[i,dst], **objfn_varg)

    # calculate the fitness proportionate selection vector
    # NOTE: this assumes all fitness values are positive!
    fps = X_scr / X_scr.sum() if fpsfn == None else fpsfn(X_scr)

    # find the index of the highest scoring individual
    X_scr_argmax = numpy.argmax(X_scr)

    # define a variable to track iteration number; initialization counts as 0.
    iter = 0

    # make a history list
    history = list()

    # store the history
    history.append([iter,
                    X_pos[X_scr_argmax,dst].copy(),
                    X_scr[X_scr_argmax]])

    #####################################
    #### Step 2) Main iteration loop ####
    #####################################

    # while stpfn is True, we haven't been told to stop
    while stpfn(iter):
        # increment iteration counter
        iter += 1

        # determine random female mates (mates[0]) and male mates (mates[1])
        # females decide where to start chromosome transmission
        # males determine new introgressions
        # rows is male/female; columns is individual
        mates = numpy.empty((2,n), dtype=numpy.int64)

        # calculate females
        mates[0,:] = numpy.searchsorted(
            numpy.cumsum(fps),
            numpy.random.uniform(0.0, 1.0, n),
            side='left' # cannot be 'right'
        )

        # calculate males
        mates[1,:] = numpy.searchsorted(
            numpy.cumsum(fps),
            numpy.random.uniform(0.0, 1.0, n),
            side='left' # cannot be 'right'
        )

        # make an empty offspring array to hold offspring
        offspring = numpy.empty(X_pos.shape, dtype=X_pos.dtype)

        # generate a binary vector to indicate whether to introgress male
        # chromosome locus; use '>=' because 1.0 is not inclusive
        introgress = (
            numpy.random.uniform(
                0.0, 1.0, (n,len(dim_sizes))
            ) >= p_rec
        ).astype(numpy.uint8) # convert bools to ints for indexing

        # generate a binary vector to indicate whether to mutate a chromosome
        # locus; use '>=' because 1.0 is not inclusive
        mutate = (
            numpy.random.uniform(
                0.0, 1.0, (n,len(dim_sizes))
            ) >= p_rec
        )

        # peform crossing over and mutations in one loop
        for i in range(n):
            for j,(st,sp,k) in enumerate(zip(dst,dsp,introgress[i])):
                # assign chromosome locus
                offspring[i,st:sp] = X_pos[mates[k,i],st:sp]
                # mutate if need be
                if mutate[i,j]:
                    # decide a random index to swap.
                    s = numpy.random.randint(st+1,sp)
                    # perform the swap
                    offspring[i,st], offspring[i,0] = \
                        offspring[i,0], offspring[i,st]

        # assign offspring to X_pos
        X_pos = offspring

        # calculate the scores for each of the positions
        if nthreads > 1:
            print("Greater than 1 thread not supported yet.")
            return
        else:
            # for each individual
            for i in range(n):
                # score the individual and store in score vector
                X_scr[i] = objfn(X_pos[i,dst], **objfn_varg)

        # calculate the fitness proportionate selection vector for next iter
        # NOTE: this assumes all fitness values are positive!
        fps = X_scr / X_scr.sum() if fpsfn == None else fpsfn(X_scr)

        # find the index of the highest scoring solution
        X_scr_argmax = numpy.argmax(X_scr)

        # store the history
        history.append([iter,
                        X_pos[X_scr_argmax,dst].copy(),
                        X_scr[X_scr_argmax]])

        if verbose:
            print("Iteration", iter, "\tgBest =", X_scr[X_scr_argmax])

    # Step 3) return positions of particles in a dictionary

    return {"history": history,
            "X_gbest_pos": X_pos[X_scr_argmax,dst].copy(),
            "X_gbest_scr": X_scr[X_scr_argmax],
            "iter": iter,
            "gbest_index": X_scr_argmax}


# WARNING: under development
def ga_elitism(objfn,
               objfn_varg,
               n,
               states,
               dim_sizes,
               p_rec,
               p_mut,
               stpfn,
               elitism = 0.5,
               seed = None,
               nthreads = 1,
               verbose = True):
    """
    Perform Integer/Categorical Particle Swarm Optimization.

    Parameters
    ----------
    objfn : function
        Objective function to optimize.
    objfn_varg : dictionary
        Objective function variable arguments.
    n : int
        Number of particles in the swarm
    states : numpy.ndarray
        1D array of possible discrete positional states for the objective
        function. The boundaries between different dimensions are determined by
        the 'dim_sizes' array.
        Example:
            states = numpy.array([0,1,2,3,1,2,3,0,2,3])
            Dimension blocks are: [0,1,2,3], [1,2,3], [0,2,3]
    dim_sizes : numpy.ndarray
        1D array of integer dimension sizes for the states vector. This is used
        to determine boundaries for possible states in the 'states' array. The
        sum of elements in this array must equal the length of 'states'.
        The dtype of dim_sizes should be numpy.uint64.
        NOTE: Assumes all state sizes are greater than 1.
        Example:
            states = numpy.array([0,1,2,3,1,2,3,0,2,3])
            dim_sizes = numpy.array([4,3,3])
    p_rec : numpy.float64
        Probability of recombination between two individuals.
    p_mut : numpy.float64
        Probability of locus mutation.
    stpfn : function
        A stop function. Should accept iteration number, ... some other things.
        Returns False to stop the iterations.
    fpsfn : function, None
        A function for generating fitness proportion fitness values (FPS). This
        function accepts a single argument. The argument is an array of scores
        from the population as scored by 'objfn'. The function should return an
        array of values. The array should sum to unity.
        If None, then a default FPS is applied. The default FPS scales
        population scores such that they all sum to unity. NOTE: This assumes
        fitness values are never negative! If negative values are possible, set
        a fpsfn to be something that accounts for this!
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
        print("Genetic algorithm with FPS engine state:")
        print("Objective function", "==================", sep='\n')
        print("Objective function =", objfn.__name__)
        print("Objective function variable arguments =", objfn_varg)
        print("Population configuration", "===================", sep='\n')
        print("Population size =", n)
        print("Chromosome length =", len(dim_sizes))
        print("Possible chromosome alleles =", states)
        print("Number of alleles per chromosome locus =", dim_sizes)
        print("Stopping function =", stpfn.__name__)
        print("Fitness proportion selection function = ", fpsfn.__name__)
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ################################
    #### Step 1) Initialization ####
    ################################

    ######################################################
    # Step 1a) allocate memory, initialize select arrays #
    ######################################################

    # create an array of dimension state stop indices
    # 1D array for dimension stop indices
    dsp = numpy.cumsum(dim_sizes)

    # create an array of dimension state start indices
    # 1D arrray for dimension start indices
    # take 'dsp' and subtract number at first index
    dst = dsp - dim_sizes[0]

    # Matrix of chromosomes. Each row of the matrix contains all states, but
    # only the indices indicated by 'dst' are used.
    # Example:
    # Chromosome:    {1} 2 3 | {2} 4 3 1 | {2} 1 3 5 4
    # This chromosome represents the solution {1,2,2}
    # Tile the input states to get 'n' individuals
    X_pos = numpy.tile(states, (n,1))

    # randomly shuffle the chromosome states to get a random population.
    for i in range(n):
        for st,sp in zip(dst,dsp):
            numpy.random.shuffle(X_pos[i,st:sp])

    # allocate a vector to contain our individual scores
    X_scr = numpy.empty(n, dtype=numpy.float64)

    # score individuals
    if nthreads > 1:
        print("Greater than 1 thread not supported yet.")
        return
    else:
        # for each individual
        for i in range(n):
            # score the individual and store in score vector
            X_scr[i] = objfn(X_pos[i,dst], **objfn_varg)

    # calculate the fitness proportionate selection vector
    # NOTE: this assumes all fitness values are positive!
    fps = X_scr / X_scr.sum() if fpsfn == None else fpsfn(X_scr)

    # find the index of the highest scoring individual
    X_scr_argmax = numpy.argmax(X_scr)

    # define a variable to track iteration number; initialization counts as 0.
    iter = 0

    # make a history list
    history = list()

    # store the history
    history.append([iter,
                    X_pos[X_scr_argmax,dst].copy(),
                    X_scr[X_scr_argmax]])

    #####################################
    #### Step 2) Main iteration loop ####
    #####################################

    # while stpfn is True, we haven't been told to stop
    while stpfn(iter):
        # increment iteration counter
        iter += 1

        # determine random female mates (mates[0]) and male mates (mates[1])
        # females decide where to start chromosome transmission
        # males determine new introgressions
        # rows is male/female; columns is individual
        mates = numpy.empty((2,n), dtype=numpy.int64)

        # calculate females
        mates[0,:] = numpy.searchsorted(
            numpy.cumsum(fps),
            numpy.random.uniform(0.0, 1.0, n),
            side='left' # cannot be 'right'
        )

        # calculate males
        mates[1,:] = numpy.searchsorted(
            numpy.cumsum(fps),
            numpy.random.uniform(0.0, 1.0, n),
            side='left' # cannot be 'right'
        )

        # make an empty offspring array to hold offspring
        offspring = numpy.empty(X_pos.shape, dtype=X_pos.dtype)

        # generate a binary vector to indicate whether to introgress male
        # chromosome locus; use '>=' because 1.0 is not inclusive
        introgress = (
            numpy.random.uniform(
                0.0, 1.0, (n,len(dim_sizes))
            ) >= p_rec
        ).astype(numpy.uint8) # convert bools to ints for indexing

        # generate a binary vector to indicate whether to mutate a chromosome
        # locus; use '>=' because 1.0 is not inclusive
        mutate = (
            numpy.random.uniform(
                0.0, 1.0, (n,len(dim_sizes))
            ) >= p_rec
        )

        # peform crossing over and mutations in one loop
        for i in range(n):
            for j,(st,sp,k) in enumerate(zip(dst,dsp,introgress[i])):
                # assign chromosome locus
                offspring[i,st:sp] = X_pos[mates[k,i],st:sp]
                # mutate if need be
                if mutate[i,j]:
                    # decide a random index to swap.
                    s = numpy.random.randint(st+1,sp)
                    # perform the swap
                    offspring[i,st], offspring[i,0] = \
                        offspring[i,0], offspring[i,st]

        # assign offspring to X_pos
        X_pos = offspring

        # calculate the scores for each of the positions
        if nthreads > 1:
            print("Greater than 1 thread not supported yet.")
            return
        else:
            # for each individual
            for i in range(n):
                # score the individual and store in score vector
                X_scr[i] = objfn(X_pos[i,dst], **objfn_varg)

        # calculate the fitness proportionate selection vector for next iter
        # NOTE: this assumes all fitness values are positive!
        fps = X_scr / X_scr.sum() if fpsfn == None else fpsfn(X_scr)

        # find the index of the highest scoring solution
        X_scr_argmax = numpy.argmax(X_scr)

        # store the history
        history.append([iter,
                        X_pos[X_scr_argmax,dst].copy(),
                        X_scr[X_scr_argmax]])

        if verbose:
            print("Iteration", iter, "\tgBest =", X_scr[X_scr_argmax])

    # Step 3) return positions of particles in a dictionary

    return {"history": history,
            "X_gbest_pos": X_pos[X_scr_argmax,dst].copy(),
            "X_gbest_scr": X_scr[X_scr_argmax],
            "iter": iter,
            "gbest_index": X_scr_argmax}
