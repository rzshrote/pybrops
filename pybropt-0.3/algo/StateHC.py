# 3rd party libraries
import numpy
import time

# our libraries
from . import HillClimber
import pybropt.util

class StateHC(HillClimber):
    """docstring for StateHC."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, search_space, name = "State Hill-Climber"):
        super(StateHC, self).__init__(name)

        pybropt.util.check_is_CategoricalSearchSpace(search_space, "search_space")
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
    def optimize(self, objfn, seed = None, nthreads = 0, minimize = True,
            verbose = False, **kwargs):
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

        # seed RNG if a seed is provided, otherwise use previous internal state.
        pybropt.util.cond_seed_rng(seed)

        # print engine state before beginning.
        if verbose:
            print("Steepest Ascent State HillClimber engine state:")
            print()
            print("Objective function")
            print("==================================================")
            print("Objective function =", objfn.__name__)
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

        # extract state array and dimension sizes from search_space
        state = self._search_space.state
        dim_size = self._search_space.dim_size

        # iteration counter
        algoiter = 0

        # create start, stop indices
        dsp = dim_size.cumsum()
        dst = dsp - dim_size

        # copy state array
        gbest_pos = state.copy()

        # randomly shuffle from state options
        # the first index in pos[st:sp] represents the chosen state
        for st,sp in zip(dst,dsp):
            numpy.random.shuffle(gbest_pos[st:sp])

        # score the initial configuration
        # get an array view of only our first indices for each dimension.
        gbest_score = objfn(gbest_pos[dst], **kwargs)

        # stash the initialization result into the history list.
        self.history_add(algoiter, [gbest_score], [gbest_pos[dst]])

        # variables for position, exchange, score matrices
        pos = None
        score = None

        # make loop exit variable
        found_local_optima = False

        # determine if we're threading
        if nthreads > 1:
            raise NotImplementedError("Greater than 1 thread not supported yet.")
        # else we're single threading
        else:
            # begin core loop
            while not found_local_optima:
                # increment iteration
                algoiter += 1

                # create state exchange matrix (local search positions)
                pos = StateHC._state_exchange_matrix(gbest_pos, dim_size)

                # calculate scores
                score = numpy.apply_along_axis(objfn, 1, pos[:,dst], **kwargs)

                # OPTIMIZE: see if any optimizations can be made in this logic chain
                # determine if any improvement has been made
                improvement = None
                if minimize:
                    improvement = (score < gbest_score).any()
                else:
                    improvement = (score > gbest_score).any()

                if improvement:
                    gbest_ix = score.argmin() if minimize else score.argmax()
                    gbest_pos = pos[gbest_ix,:].copy() # OPTIMIZE: see if copy is needed
                    gbest_score = score[gbest_ix]
                else:
                    found_local_optima = True

                # add history
                self.history_add(algoiter, score, pos[:,dst])

                if verbose:
                    print("Iteration", algoiter, "\tscore =", gbest_score)

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def _state_exchange_matrix(vec, size):
        """
        vec : numpy.ndarray
            Packed vector of search space states.
            Example:
                vec = numpy.array([1,2,3,1,2,3,4])
        size : numpy.ndarray
            Vector of dimension sizes of 'vec'
            Example:
                vec = numpy.array([1,2,3,1,2,3,4])
                size = numpy.array([3,4])
        Returns
        -------
        out_vec : numpy.ndarray
            A matrix of state exchange positions.

        Example
        -------
        >>> vec = numpy.array([1,2,3,1,2,3,4])
        >>> size = numpy.array([3,4])
        >>> _state_exchange_matrix(vec, size)
        array([[2, 1, 3, 1, 2, 3, 4],
               [3, 2, 1, 1, 2, 3, 4],
               [1, 2, 3, 2, 1, 3, 4],
               [1, 2, 3, 3, 2, 1, 4],
               [1, 2, 3, 4, 2, 3, 1]])
        """
        # calculate number of rows we'll need
        nrow = (size-1).sum()

        # calculate start indices
        start = size.cumsum() - size

        # tile vec
        out_vec = numpy.tile(vec, (nrow,1))

        # row counter
        k = 0

        # for each dimension group, exchange each element
        for st,l in zip(start, size):
            for i in range(1, l):
                # exchange values
                out_vec[k,st], out_vec[k,st+i] = out_vec[k,st+i], out_vec[k,st]
                # increment counter
                k += 1

        return out_vec
