import numpy
import time

from .algo_opt import algo_opt

def _set_exchange_matrix(vec, exch):
    out_vec = numpy.tile(vec, (len(vec)*len(exch),1))
    for i in range(len(vec)):
        out_vec[i*len(exch):(i+1)*len(exch),i] = exch
    out_exch = numpy.tile(exch, (len(vec)*len(exch),1))
    for i in range(len(exch)):
        out_exch[i::len(exch),i] = vec
    return out_vec, out_exch

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
        seed = numpy.uint32(time.time())
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

    # iteration counter
    iter = 0

    tmp = states.copy()         # copy state list
    numpy.random.shuffle(tmp)   # shuffle the copy
    gbest_pos = tmp[:k].copy()  # extract start position
    gbest_exch = tmp[k:].copy() # extract start exchange vector

    # score the initial configuration
    gbest_score = objfn(gbest_pos, **objfn_varg)

    # make and add history
    opt = algo_opt(
        dtype_iter = numpy.dtype(type(iter)),
        dtype_score = numpy.dtype(type(gbest_score)),
        dtype_pos = gbest_pos.dtype
    )
    opt.history_add(iter, [gbest_score], gbest_pos)

    # variables for position, exchange, score matrices
    pos = None
    exch = None
    score = None

    # print verbose statements
    if verbose:
        print("Iteration", iter, "\tscore =", score)

    # make exit variable
    found_local_optima = False

    # begin exchange process
    while not found_local_optima:
        # increment iteration
        iter += 1

        pos, exch = _set_exchange_matrix(gbest_pos, gbest_exch)

        # score positions that have been created (threading applicable)
        if nthreads > 1:
            # https://stackoverflow.com/questions/3033952/threading-pool-similar-to-the-multiprocessing-pool
            # https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python
            raise RuntimeError("Multithreading not supported yet...")
        else:
            # evaluate the objective function
            score = numpy.apply_along_axis(objfn, 1, pos, **objfn_varg)

        if numpy.any(score > gbest_score):
            gbest_ix = score.argmax()
            gbest_pos = pos[gbest_ix,:].copy()
            gbest_exch = exch[gbest_ix,:].copy()
            gbest_score = score[gbest_ix]
        else:
            found_local_optima = True

        # add history
        opt.history_add(iter, score, pos)

        print("Iteration", iter, "\tscore =", gbest_score)

    # return a history list from the hillclimb
    return opt



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
        seed = numpy.uint32(time.time())
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
    history = list()

    # iteration counter
    iter = 0

    # create array of dimension stop indices
    dsp = numpy.cumsum(dim_sizes)

    # create array of dimension start indices
    dst = dsp - dim_sizes[0]

    # copy states array
    pos = states.copy()

    # randomly shuffle from state options
    # the first index in pos[st:sp] represents the chosen state
    for st,sp in zip(dst,dsp):
        numpy.random.shuffle(pos[st:sp])

    # score the initial configuration
    # get an array view of only our first indices for each dimension.
    score = objfn(pos[dst], **objfn_varg)

    # stash the initialization result into the history list.
    history.append([iter, pos[dst].copy(), score])

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
            new_scr = score     # variable to store better scores if found

            # for each dimension
            for st,sp in zip(dst,dsp):
                # for each alternative state within the dimension
                for s in range(1, len(pos[st:sp])):
                    # swap to new state
                    pos[st:sp][0], pos[st:sp][s] = \
                        pos[st:sp][s], pos[st:sp][0]

                    # score position
                    scr = objfn(pos[dst], **objfn_varg)

                    # if better a better score is generated
                    if scr > new_scr:
                        new_scr = scr   # set new best score
                        st_best = st    # record best start index
                        sp_best = sp    # record best stop index
                        s_best = s      # record s_best

                    # swap back to original state
                    pos[st:sp][0], pos[st:sp][s] = \
                        pos[st:sp][s], pos[st:sp][0]

            # if we found something better (e_best got set)
            if st_best != None:
                # swap to new best position
                pos[st_best:sp_best][0], pos[st_best:sp_best][s_best] = \
                    pos[st_best:sp_best][s_best], pos[st_best:sp_best][0]
                # update X score
                score = new_scr
                # append result to list
                history.append([iter, pos[dst].copy(), score])
            # else we didn't and no further improvements can be made
            else:
                # we've found a local optima; set variable to exit loop
                found_local_optima = True

            if verbose:
                print("Iteration", iter, "\tscore =", score)

    # return a history list from the hillclimb
    return {"history": history,
            "X_gbest_pos": pos[dst].copy(),
            "X_gbest_scr": score}
