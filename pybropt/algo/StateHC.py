class StateHC(HillClimber):
    """docstring for StateHC."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################

    def __init__(self, search_space):
        super(StateHC, self).__init__()

        check_is_CategoricalSearchSpace(search_space, "search_space")
        self._search_space = search_space

    ############################################################################
    ################################ Properties ################################
    ############################################################################

    def search_space():
        doc = "The search_space property."
        def fget(self):
            return self._search_space
        def fset(self, value):
            self._search_space = value
        def fdel(self):
            del self._search_space
        return locals()
    search_space = property(**search_space())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    def _state_exchange_matrix(vec, exch)

    def optimize(self, objfn, stpfn, seed = None, nthreads = None,
            verbose = False, *args, **kwargs):
        """
        Perform a hillclimbing using steepest ascent strategy. Variables are
        represented by dimensions with possible states. States are iteratively
        tested one site at a time.

        Parameters
        ==========
        objfn : function
            Objective function to optimize.
        stpfn : function
            A stop function. Should accept iteration number.
            Returns False to stop the iterations.
        seed : int, default = None.
            Random number seed.
        nthreads : int, default = None
            Number of threads. If nthreads is not None, this instantiates a
            multiprocessing pool with the given number of processes.
        verbose : boolean
            Show verbose messages.
        *args
            Additional arguments for objfn.
        **kwargs
            Additional keyword arguments for objfn.

        Returns
        -------
        None
        """
        #####################
        # Process arguments #
        #####################

        # if no seed is provided, use time as random seed
        if seed is None:
            seed = numpy.uint32(time.time())
        # seed random number generator
        numpy.random.seed(seed)

        # print engine state before beginning.
        if verbose:
            print("Steepest Ascent State HillClimber engine state:")
            print()
            print("Objective function")
            print("==================================================")
            print("Objective function =", objfn.__name__)
            print("Objective function arguments:", str(args))
            print("Objective function keyword arguments:", str(kwargs))
            print()
            print("HillClimber configuration")
            print("==================================================")
            print("Number of dimensions =", self._search_space.ndim)
            print("Possible position states =", self._search_space)
            print()
            print("Miscellaneous")
            print("==================================================")
            print("numpy.random seed =", seed)
            print("Number of threads =", nthreads)

        ####################
        # Computation part #
        ####################

        # extract states array and dimension sizes from search_space
        states = self._search_space.state_flat
        dim_sizes = self._search_space.dim_size

        # iteration counter
        iter = 0

        # create start, stop indices
        dsp = dim_sizes.cumsum()
        dst = dsp - dim_sizes

        # copy states array
        pos = states.copy()

        # randomly shuffle from state options
        # the first index in pos[st:sp] represents the chosen state
        for st,sp in zip(dst,dsp):
            numpy.random.shuffle(pos[st:sp])

        # score the initial configuration
        # get an array view of only our first indices for each dimension.
        score = objfn(pos[dst], *args, **kwargs)

        # stash the initialization result into the history list.
        self.history_add(iter, [score], pos[dst])

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
                        scr = objfn(pos[dst], *args, *kwargs)

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
                    self.history_add(iter, score, pos[dst])
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
