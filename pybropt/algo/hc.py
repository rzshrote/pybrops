import numpy
import time

def hc_sa_set(objfn,
              objfn_varg,
              k,
              states,
              seed = None,
              nthreads = 1,
              verbose = True):
    """
    Perform a hillclimbing using steepest ascent strategy. A set of size
    'k' is used. Single states within the set are iteratively exchanged
    with one state not within the set.

    Parameters
    ==========
    objfn : function
        Objective function.
    objfn_varg : dictionary
        Dictionary containing variable arguments
    k : int
        Number of states to select and conduct hillclimb on.
    states : numpy.ndarray
        1D array of possible discrete positional states for the objective
        function.
        Example:
            states = numpy.array([0,1,2,3,4,5,6,7,8,9])
    seed : int
        Random number seed.
    nthreads : int
        Number of threads.
    verbose : boolean
        Show verbose messages
    """
    #####################
    # Process arguments #
    #####################

    # if no seed is provided, use time as random seed
    if seed == None:
        seed = time.time()
    # seed random number generator
    numpy.random.seed(seed)

    # print engine state before beginning.
    if verbose:
        print("Engine state:")
        print("Objective function =", objfn.__name__)
        print("Objective function variable arguments =")
        print(objfn_varg)
        print("Subset size =", k)
        print("Possible position states =")
        print(states)
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ####################
    # Computation part #
    ####################

    # create list to store results from the hillclimb.
    hc_result = list()

    # iteration counter
    iter = 0

    # make X position array
    X_pos = states.copy()

    # shuffle the X_pbest_pos for random initialization
    numpy.random.shuffle(X_pos)

    # score the initial configuration
    X_scr = objfn(X_pos[0:k], **objfn_varg)

    # stash the initialization result into the hc_result list.
    hc_result.append([iter, X_pos[0:k].copy(), X_scr])

    # make exit variable
    found_local_optima = False

    # determine if we're threading
    if nthreads > 1:
        print("Greater than 1 thread not supported yet.")
        return
    # else we're single threading
    else:
        # begin exchange process
        while not found_local_optima:
            # increment iteration
            iter += 1

            e_best = None       # index for chosen elements best swap
            s_best = None       # index for available elements best swap
            new_scr = X_scr     # variable to store better scores if found

            for e in range(k):                              # for each chosen
                for s in range(k, len(X_pos)):              # for each available
                    X_pos[e], X_pos[s] = X_pos[s], X_pos[e] # swap to new
                    scr = objfn(X_pos[0:k], **objfn_varg)   # score matrix
                    if scr > new_scr:                       # if better score
                        new_scr = scr                       # set new best score
                        e_best = e                          # record e_best
                        s_best = s                          # record s_best
                    X_pos[e], X_pos[s] = X_pos[s], X_pos[e] # swap to original

            # if we found something better (e_best got set)
            if e_best != None:
                # swap to new best position
                X_pos[e_best], X_pos[s_best] = X_pos[s_best], X_pos[e_best]
                # update X score
                X_scr = new_scr
                # append result to list
                hc_result.append([iter, X_pos[0:k].copy(), X_scr])
            # else we didn't and no further improvements can be made
            else:
                # we've found a local optima
                found_local_optima = True

            if verbose:
                print("Iteration", iter, "\tX_scr =", X_scr)

    # return a history list from the hillclimb
    return hc_result



def hc_sa_state(objfn,
                objfn_varg,
                states,
                dim_sizes,
                seed = None,
                nthreads = 1,
                verbose = True):
    """
    Perform a hillclimbing using steepest ascent strategy. Variables are
    represented by dimensions with possible states. States are iteratively
    tested one site at a time.

    Parameters
    ==========
    objfn : function
        Objective function.
    objfn_varg : dictionary
        Dictionary containing variable arguments
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
    seed : int
        Random number seed.
    nthreads : int
        Number of threads.
    verbose : boolean
        Show verbose messages
    """
    #####################
    # Process arguments #
    #####################

    # if no seed is provided, use time as random seed
    if seed == None:
        seed = time.time()
    # seed random number generator
    numpy.random.seed(seed)

    # print engine state before beginning.
    if verbose:
        print("Engine state:")
        print("Objective function =", objfn.__name__)
        print("Objective function variable arguments =")
        print(objfn_varg)
        print("Possible position states =")
        print(states)
        print("Dimension sizes =")
        print(dim_sizes)
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)

    ####################
    # Computation part #
    ####################

    # create list to store results from the hillclimb.
    hc_result = list()

    # iteration counter
    iter = 0

    # create array of dimension stop indices
    dsp = numpy.cumsum(dim_sizes)

    # create array of dimension start indices
    dst = dsp - dim_sizes[0]

    # copy states array
    X_pos = states.copy()

    # randomly shuffle from state options
    # the first index in X_pos[st:sp] represents the chosen state
    for st,sp in zip(dst,dsp):
        numpy.random.shuffle(X_pos[st:sp])

    # score the initial configuration
    # get an array view of only our first indices for each dimension.
    X_scr = objfn(X_pos[dst], **objfn_varg)

    # stash the initialization result into the hc_result list.
    hc_result.append([iter, X_pos[dst].copy(), X_scr])

    # make loop exit variable
    found_local_optima = False

    # determine if we're threading
    if nthreads > 1:
        print("Greater than 1 thread not supported yet.")
        return
    # else we're single threading
    else:
        # begin exchange process
        while not found_local_optima:
            # increment iteration
            iter += 1

            st_best = None
            sp_best = None
            s_best = None       # index for available elements best swap
            new_scr = X_scr     # variable to store better scores if found

            # for each dimension
            for st,sp in zip(dst,dsp):
                # for each alternative state within the dimension
                for s in range(1, len(X_pos[st:sp])):
                    # swap to new state
                    X_pos[st:sp][0], X_pos[st:sp][s] = \
                        X_pos[st:sp][s], X_pos[st:sp][0]

                    # score position
                    scr = objfn(X_pos[dst], **objfn_varg)

                    # if better a better score is generated
                    if scr > new_scr:
                        new_scr = scr   # set new best score
                        st_best = st    # record best start index
                        sp_best = sp    # record best stop index
                        s_best = s      # record s_best

                    # swap back to original state
                    X_pos[st:sp][0], X_pos[st:sp][s] = \
                        X_pos[st:sp][s], X_pos[st:sp][0]

            # if we found something better (e_best got set)
            if st_best != None:
                # swap to new best position
                X_pos[st_best:sp_best][0], X_pos[st_best:sp_best][s_best] = \
                    X_pos[st_best:sp_best][s_best], X_pos[st_best:sp_best][0]
                # update X score
                X_scr = new_scr
                # append result to list
                hc_result.append([iter, X_pos[dst].copy(), X_scr])
            # else we didn't and no further improvements can be made
            else:
                # we've found a local optima; set variable to exit loop
                found_local_optima = True

            if verbose:
                print("Iteration", iter, "\tX_scr =", X_scr)

    # return a history list from the hillclimb
    return hc_result
