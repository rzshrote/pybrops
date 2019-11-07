import numpy
import time

def icpso(objfn,
          objfn_varg,
          n,
          states,
          dim_sizes,
          inertia_wt,
          accel_coeff_pbest,
          accel_coeff_gbest,
          scale_factor,
          stpfn,
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
        seed = time.time()

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

    # create an array of dimension state stop indices
    # 1D array for dimension stop indices
    dsp = numpy.cumsum(dim_sizes)

    # create an array of dimension state start indices
    # 1D arrray for dimension start indices
    # take 'dsp' and subtract number at first index
    dst = dsp - dim_sizes[0]

    # particle population matrix shape
    pshape = (n, len(dim_sizes))

    # initialize random X probability positions with matrix shape 'pshape'
    # use uniform distribution; dtype of this array is numpy.float64
    X_pos = numpy.random.uniform(0.0, 1.0, pshape)

    # allocate memory for an empty sample array
    # holds samples from X_pos vectors
    X_pos_smpl = numpy.empty((n, len(dim_sizes)), dtype=states.dtype)

    # initialize random X velocity with matrix shape 'pshape'
    X_vel = numpy.random.uniform(0.0, 1.0, pshape)

    # set X_best_pos to None, it will be initialized in Step 1c
    # individual particle personal best values
    X_pbest_pos = None

    # personal best samples, not used by the algorithm per se, but updated and
    # returned
    X_pbest_smpl = None

    # set X_pbest_scr to None, it will be initialized in Step 1c
    # X personal best scores for the positions of the 'n' particles.
    X_pbest_scr = None

    # set X_gbest_pos to None, it will be initialized in Step 1c
    # global best position for the swarm
    # note: this is a view of X_pbest_pos
    X_gbest_pos = None

    # global best samples, not used by the algorithm per se, but updated and
    # returned
    X_gbest_smpl = None

    # global best score for the swarm
    X_gbest_scr = None

    ############################################################
    # Step 1b) restrict limits of values and scale to sum to 1 #
    ############################################################

    # Map probabilities >1 to their edge, 1
    # Logic: if value > 1 then set to 1
    X_pos[X_pos > 1] = 1
    X_vel[X_vel > 1] = 1

    # adjust values to sum to 1 in highly obfuscated manner
    # Logic: 1. zip dimension starts and stops into iterable object (single use)
    #        2. iterate through start, stop positions
    #        3. modify the matrix in blocks by dividing row-wise the sum of the
    #           row-wise block sum
    X_pos[:,st:sp] /= X_pos[:,st:sp].sum(1)[:,None] for st,sp in zip(dst, dsp)
    X_vel[:,st:sp] /= X_vel[:,st:sp].sum(1)[:,None] for st,sp in zip(dst, dsp)

    ############################################################################
    # Step 1c) initialize best positions, take samples, calculate start scores #
    ############################################################################

    # copy randomly initialized position matrix to an X_pbest_pos matrix
    X_pbest_pos = X_pos.copy()

    # generate sample from X personal best position
    for i in range(n):                          # for each row
        X_pos_smpl[i,j] = numpy.random.choice(  # set sample to random choice
                states[st:sp],                  # with these states
                size=1,                         # only sample 1
                p=X_pbest_pos[i,st:sp]          # with these probabilities
            ) for j,(st,sp) in enumerate(       # for each index, start, stop
                zip(dst, dsp)                   # zip up start, stop to iterate
            )

    # copy personal best samples
    X_pbest_smpl = X_pos_smpl.copy()

    # score samples that have been created (threading applicable)
    if nthreads > 1:
        # https://stackoverflow.com/questions/3033952/threading-pool-similar-to-the-multiprocessing-pool
        # https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python
        print("Multithreading not supported yet...")
        return
    else:
        # evaluate the objective function
        X_pbest_scr = numpy.apply_along_axis(objfn, 0, X_pos_smpl, **objfn_varg)

    # find the index of the highest scoring
    argmax_index = numpy.argmax(X_pbest_scr)

    # since up until this point, we're agnostic about what previous scores have
    # been, we set the X global best position to the view of the argmax row.
    X_gbest_pos = X_pbest_pos[argmax_index,:]

    # set the X_gbest_scr to X_pbest_scr[argmax]
    X_gbest_scr = X_pbest_scr[argmax_index]

    #####################################
    #### Step 2) Main iteration loop ####
    #####################################

    # define a variable to track iteration number; initialization counts as 0.
    iter = 0

    # while stpfn is True, we haven't been told to stop
    while stpfn(iter):
        #####################################################
        # Step 2a) update particle velocities and positions #
        #####################################################

        # generate random acceleration coefficients for personal and global
        acoef_pbest = numpy.random.uniform(low=0,high=accel_coeff_pbest,size=n)
        acoef_gbest = numpy.random.uniform(low=0,high=accel_coeff_gbest,size=n)

        # calculate the new velocities
        X_vel = (inertia_wt * X_vel) + \
                (acoef_pbest[:,None] * (X_pbest_pos - X_pos)) + \
                (acoef_gbest[:,None] * (X_gbest_pos - X_pos))

        # restrict limits of velocity values ...
        X_vel[X_vel > 1] = 1

        # ... and scale to sum to 1 for each dimension
        X_vel[:,st:sp] /= X_vel[:,st:sp].sum(1)[:,None] for st,sp in zip(dst, dsp)

        # add velocity vectors to X positions
        X_pos += X_vel

        # restrict limits of position values ...
        X_pos[X_pos > 1] = 1

        # ... and scale to sum to 1 for each dimension
        X_pos[:,st:sp] /= X_pos[:,st:sp].sum(1)[:,None] for st,sp in zip(dst, dsp)

        ################################################################
        # Step 2b) generate samples from X_pos and score those samples #
        ################################################################
        # TODO: maybe replace this with for one-liner?
        # generate sample from X position
        for i in range(n):                          # for each row
            X_pos_smpl[i,j] = numpy.random.choice(  # set sample to random choice
                    states[st:sp],                  # with these states
                    size=1,                         # only sample 1
                    p=X_pos[i,st:sp]                # with these probabilities
                ) for j,(st,sp) in enumerate(       # for each index, start, stop
                    zip(dst, dsp)                   # zip up start, stop to iterate
                )

        # score samples (threading applicable)
        if nthreads > 1:
            print("Multithreading not supported yet...")
            return
        else:
            # evaluate the objective function for each sample
            scores = numpy.apply_along_axis(objfn, 0, X_pos_smpl, **objfn_varg)

        # determine if scores improve
        for i in range(n):                  # for each score
            if scores[i] > X_pbest_scr[i]:  # if score beats previous best
                X_pbest_scr[i] = scores[i]  # update best score

                # update best samples (not used by algorithm, for output)
                X_pbest_smpl[i,:] = X_pos_smpl[i,:]

                # set new X_pbest_pos that is biased for values
                for j,(st,sp) in enumerate(zip(dst,dsp)):
                    # make mask for when potential state == sample state
                    mask = (states[st:sp] == X_pos_smpl[i,j])

                    # find the index of the state making new best score
                    smpl_index = numpy.flatnonzero(mask)[0]

                    # calculate positive bias for sample values
                    # multiply the sum of non-sample positions by (1 -
                    # scale_factor)
                    pos_bias = (1 - scale_factor) * \
                                numpy.delete(
                                    X_pos[i,st:sp],
                                    smpl_index
                                ).sum()

                    # bias and assign new personal best position
                    X_pbest_pos[i,st:sp] = numpy.where(
                        mask,                           # if in sample
                        X_pos[i,st:sp] + pos_bias,      # add positive bias
                        scale_factor * X_pos[i,st:sp]   # else add neg bias
                    )

        ##########################################
        # Step 2c) update the global best values #
        ##########################################

        # find the index of the highest scoring
        argmax_index = numpy.argmax(X_pbest_scr)

        # set the X global best position to the view of the argmax row.
        X_gbest_pos = X_pbest_pos[argmax_index,:]

        # set the X global best sample (not for use by algorithm, for output)
        X_gbest_smpl = X_pbest_smpl[argmax_index,:]

        # set the X_gbest_scr to X_pbest_scr[argmax]
        X_gbest_scr = X_pbest_scr[argmax_index]

        # increment iteration counter
        iter += 1

        if verbose:
            print("Iteration", iter, "\tgBest =", X_gbest_scr)

    # Step 3) return positions of particles in a dictionary

    return {"X_pos": X_pos,
            "X_vel": X_vel,
            "X_pbest_pos": X_pbest_pos,
            "X_pbest_smpl": X_pbest_smpl,
            "X_pbest_scr": X_pbest_scr,
            "X_gbest_pos": X_gbest_pos.copy(),
            "X_gbest_smpl": X_gbest_smpl.copy(),
            "X_gbest_scr": X_gbest_scr,
            "iter": iter,
            "gbest_index": numpy.argmax(X_pbest_scr)}
