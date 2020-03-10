import numpy
import time

from .icpso_opt import icpso_opt

def icpso(objfn, objfn_varg,
          n, states, dim_sizes, inertia_wt, accel_coeff_pbest,
          accel_coeff_gbest, scale_factor,
          stpfn, seed = None, nthreads = 1, verbose = True):
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
        Example:
            states = numpy.array([0,1,2,3,1,2,3,0,2,3])
            dim_sizes = numpy.array([4,3,3])
    inertia_wt : float
        Inertial weight for position update. The previous velocity of the
        particle is multiplied by this coefficient.
    accel_coeff_pbest : float
        Acceleration coefficient 1 for position update. The difference between
        the personal best of the particle minus the current position is
        multiplied by a sample from a random uniform distribution
        [0, accel_coeff1).
    accel_coeff_gbest : float
        Acceleration coefficient 2 for position update. The difference between
        the global best of the particle minus the current position is multiplied
        by a sample from a random uniform distribution [0, accel_coeff1).
    scale_factor : float
        Scaling factor for determining magnitude of shifts in distributions
        when a personal or global best position is updated. The range of this
        coefficient is [0,1).
    stpfn : function
        A stop function. Should accept iteration number, ... some other things.
        Returns False to stop the iterations.
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
        print("Integer/Categorical PSO engine state:")
        print("Objective function", "==================", sep='\n')
        print("Objective function =", objfn.__name__)
        print("Objective function variable arguments =", objfn_varg)
        print("Swarm configuration", "===================", sep='\n')
        print("Number of particles =", n)
        print("Number of dimensions =", len(dim_sizes))
        print("Possible position states =", states)
        print("Dimension number of states =", dim_sizes)
        print("Velocity inertial weight =", inertia_wt)
        print("Personal best acceleration coefficient =", accel_coeff1)
        print("Global best acceleration coefficient =", accel_coeff2)
        print("Stopping function =", stpfn.__name__)
        print("Miscellaneous", "=============", sep='\n')
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ################################
    #### Step 1) Initialization ####
    ################################

    ######################################################
    # Step 1a) allocate memory, initialize select arrays #
    ######################################################

    # create arrays of dimension state start and stop indices
    dst = dim_sizes.cumsum() - dim_sizes    # start indices
    dsp = dim_sizes.cumsum()                # stop indices

    #################################################
    # generate position probability array
    #################################################

    # initialize random X probability positions with n rows, len(states) columns
    pos = numpy.random.uniform(     # use uniform distribution (dtype=float64)
        low = 0.0,                  # min probability (inclusive)
        high = 1.0,                 # max probability (exclusive)
        size = (n, len(states))     # n rows by len(states) columns
    )

    #################################################
    # generate velocity vectors
    #################################################
    # initialize random X velocity with matrix shape 'pshape'
    vel = numpy.random.uniform( # use uniform distribution (dtype=float64)
        low = 0.0,             # min trajectory strength (inclusive)
        high = 1.0,             # max trajectory strength (exclusive)
        size = pos.shape        # use same shape as pos
    )

    # for each group adjust values to sum to 1
    for st,sp in zip(dst, dsp):
        pos[:,st:sp] /= pos[:,st:sp].sum(1)[:,None]
        vel[:,st:sp] /= vel[:,st:sp].sum(1)[:,None] # not needed???

    #################################################
    # generate sample array
    #################################################

    # allocate memory for an empty sample array
    smpl = numpy.empty(                 # holds samples from pos vector
        shape = (n, len(dim_sizes)),    # n rows by len(dim_sizes) columns
        dtype = states.dtype            # use same dtype as states
    )

    # take samples for each row, for each grouping
    for i in range(n):
        # for each index, start, stop (for each group)
        for j,(st,sp) in enumerate(zip(dst, dsp)):
            smpl[i,j] = numpy.random.choice(    # set sample to random choice
                states[st:sp],                  # with these states
                size = 1,                       # only sample 1
                p = pos[i,st:sp]                # with these probabilities
            )

    #################################################
    # declare score matrix variable
    #################################################

    score = None

    #################################################
    # initialize personal and global best positions and score variables
    #################################################

    # personal best variables
    pbest_pos = None    # personal best positions (probability matrix)
    pbest_score = None  # personal best scores (vector; dtype depends on objfn)

    # global best variables
    gbest_pos = None    # global best position (probability vector); a matrix view
    gbest_score = None  # global best score (dtype depends on objfn)

    #################################################
    # score initial positions
    #################################################

    # score samples that have been created (threading applicable)
    if nthreads > 1:
        # https://stackoverflow.com/questions/3033952/threading-pool-similar-to-the-multiprocessing-pool
        # https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python
        raise RuntimeError("Multithreading not supported yet...")
    else:
        # evaluate the objective function
        score = numpy.apply_along_axis(objfn, 1, smpl, **objfn_varg)


    #################################################
    # set personal best positions and scores
    #################################################

    # copy randomly initialized position matrix to an pbest_pos matrix
    pbest_pos = pos.copy()

    # copy initial scores
    pbest_score = score.copy()

    #################################################
    # identify global best
    #################################################

    # find the index of the highest scoring particle
    gbest_ix = pbest_score.argmax()

    # get matrix of global best position
    gbest_pos = pbest_pos[gbest_ix,:].copy()

    # set the gbest_score to pbest_score[argmax]
    gbest_score = pbest_score[gbest_ix]

    #################################################
    # define iteration counter global best
    #################################################

    # define a variable to track iteration number; initialization counts as 0.
    iter = 0

    #################################################
    # create, append history
    #################################################

    # create history object
    opt = icpso_opt(
        ssize = n,
        inertia_wt = inertia_wt,
        pbest_comp = accel_coeff_pbest,
        gbest_comp = accel_coeff_gbest,
        scale_factor = scale_factor,
        dtype_iter = numpy.dtype('int64'),
        dtype_score = score.dtype,
        dtype_pos = pos.dtype,
        dtype_vel = vel.dtype,
        dtype_smpl = smpl.dtype
    )

    # add history to history object
    opt.history_add(iter, score, pos, vel, smpl)

    #####################################
    #### Step 2) Main iteration loop ####
    #####################################

    # while stpfn is True, we haven't been told to stop
    while stpfn(iter):
        #####################################################
        # Step 2a) update particle velocities and positions #
        #####################################################

        # increment iteration counter
        iter += 1

        # generate random acceleration coefficients for personal and global
        acoef_pbest = numpy.random.uniform(
            low = 0.0,
            high = accel_coeff_pbest,
            size = n
        )
        acoef_gbest = numpy.random.uniform(
            low = 0.0,
            high = accel_coeff_gbest,
            size = n
        )

        # calculate the new velocities
        vel = (inertia_wt * vel) + \
                (acoef_pbest[:,None] * (pbest_pos - pos)) + \
                (acoef_gbest[:,None] * (gbest_pos - pos))

        # map to edge
        vel[vel > 1.0] = 1.0
        vel[vel < 0.0] = 0.0

        # scale velocities
        for st,sp in zip(dst, dsp):
            vel[:,st:sp] /= vel[:,st:sp].sum(1)[:,None] # do not scale???

        # add velocity vectors to X positions
        pos += vel

        # Map position probabilities to their edge
        # Logic: if pos > 1, set to 1; if pos < 0, set to 0
        pos[pos > 1.0] = 1.0
        pos[pos < 0.0] = 0.0

        # for each group adjust pos values to sum to 1
        for st,sp in zip(dst, dsp):
            pos[:,st:sp] /= pos[:,st:sp].sum(1)[:,None]

        ################################################################
        # Step 2b) generate samples from pos and score those samples #
        ################################################################
        # take samples for each row, for each grouping
        for i in range(n):
            # for each index, start, stop (for each group)
            for j,(st,sp) in enumerate(zip(dst, dsp)):
                smpl[i,j] = numpy.random.choice(    # set sample to random choice
                    states[st:sp],                  # with these states
                    size = 1,                       # only sample 1
                    p = pos[i,st:sp]                # with these probabilities
                )

        # score samples (threading applicable)
        if nthreads > 1:
            raise RuntimeError("Multithreading not supported yet...")
        else:
            # evaluate the objective function for each sample
            score = numpy.apply_along_axis(objfn, 1, smpl, **objfn_varg)

        ################################################################
        # add history
        ################################################################
        opt.history_add(iter, score, pos, vel, smpl)

        ################################################################
        # determine if scores improve
        ################################################################
        # create an improvement mask
        improve_mask = score > pbest_score

        for i in range(n):                  # for each score
            if improve_mask[i]:             # if score beats previous best
                pbest_score[i] = score[i]   # update best score

                # set new pbest_pos that is biased for values
                for j,(st,sp) in enumerate(zip(dst,dsp)):
                    # find where sample state is not equal to potential state
                    inv_mask = states[st:sp] != smpl[i,j]

                    # calculate positive bias for sample values
                    pos_bias = ((1-scale_factor) * pos[i,st:sp] * inv_mask).sum()

                    # bias and assign new personal best position
                    pbest_pos[i,st:sp] = numpy.where(
                        inv_mask,                     # if element not sample
                        scale_factor * pos[i,st:sp],  # add neg bias
                        pos[i,st:sp] + pos_bias       # else add positive bias
                    )

        ##########################################
        # Step 2c) update the global best values #
        ##########################################

        # find the index of the highest scoring
        gbest_ix = pbest_score.argmax()

        # set the X global best position to the view of the argmax row.
        gbest_pos = pbest_pos[gbest_ix,:].copy()

        # set the gbest_score to pbest_score[argmax]
        gbest_score = pbest_score[gbest_ix]

        # if verbose:
        print("Iteration", iter, "\tgBest =", gbest_score)

    # Step 3) return optimization history
    return opt
