import numpy
import time

def hc(n,
       states,
       objfn,
       objfn_varg,
       seed = None,
       nthreads = 1,
       verbose = True):
    """
    n : int
        Number of hillclimbs.
    states : numpy.ndarray of numpy.ndarray objects
        Possible states for the hillclimb.
    objfn : function
        Objective function.
    objfn_varg : dictionary
        Dictionary containing variable arguments
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
        print("Hillclimb number =", n)
        print("Possible position states =")
        print(states)
        print("Objective function =", objfn.__name__)
        print("Objective function variable arguments =")
        print(objfn_varg)
        print("numpy.random seed =", seed)
        print("Number of threads =", nthreads)
        # printing verbose is redundant

    ####################
    # Computation part #
    ####################

    # create list to store results from each hillclimb run.
    hc_results = list()

    # determine if we're threading
    if nthreads > 1:
        print("Greater than 1 thread not supported yet.")
        return
    # else we're single threading
    else:
        # for each hillclimb run.
        for num_hc in range(n):
            # iteration counter
            iter = 0

            # make a list to store single hillclimb result
            hc_result = list()

            # initialize X_pos to zeros
            X_pos_best = numpy.zeros(len(states), dtype=states[0].dtype)

            # for each dimension, choose a random state.
            for dim in range(len(states)):
                X_pos_best[dim] = numpy.random.choice(states[dim])

            # score the random starting position
            X_score_best = objfn(X_pos_best, **objfn_varg)

            # stash the initialization result into the hc_result list.
            hc_result.append([iter, X_pos_best, X_score_best])

            # make exit variable
            found_local_optima = False

            # begin exchange process
            while !found_local_optima:
                # increment iteration
                iter += 1

                # make climbing copies of postion and score for comparisions
                X_pos_climb = X_pos_best.copy()
                X_score_climb = X_score_best

                # make climbing best copies of position and score to compare
                X_pos_climb_best = X_pos_best.copy()
                X_score_climb_best = X_score_best

                # for each dimension
                for dim in range(len(states)):
                    # for each possible state in dimension 'dim'
                    for state in range(len(states[dim])):
                        # if the state is same as the current best (what we're
                        # trying to improve by climbing higher), skip to avoid
                        # unnecessary computation
                        if states[dim][state] == X_pos_best[dim]:
                            continue

                        # set the working copy of position array to new state.
                        X_pos_climb[dim] = state[dim][state]

                        # score the working climb copy
                        X_score_climb = objfn(X_pos_climb, **objfn_varg)

                        # if we found a score that exceeds the best, store it
                        if X_score_climb > X_score_climb_best:
                            # copy array to prevent pointer issues
                            X_pos_climb_best = X_score_climb.copy()

                            # copy best climbing score
                            X_score_climb_best = X_score_climb

                    # end loop by resetting the working copy of the position
                    X_pos_climb[dim] = X_pos_best[dim]
                    #X_score_climb = X_score_best # unnecessary

                # if we found a better score in the climb, update best position,
                # update best score, append iteration result
                if X_score_climb_best > X_score_best:
                    # update best position
                    X_pos_best = X_pos_climb_best.copy()

                    # update best score
                    X_score_best = X_score_climb_best

                    # stash the iteration result into the hc_result list.
                    hc_result.append([iter, X_pos_best, X_score_best])

                # else we've reached a local maximum, so set found_local_optima
                else:
                    found_local_optima = True
            # END while

            # append hc_result list to hc_results (all hillclimb runs) list
            hc_results.append(hc_result)
        # END for

        # return results
        return hc_results
