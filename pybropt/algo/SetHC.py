# 3rd party libraries
import numpy

# our libraries
from . import HillClimber
import pybropt.util

class SetHC(HillClimber):
    """docstring for SetHC."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, k, search_space, name = "Set Hill-Climber"):
        # call super constructor
        super(SetHC, self).__init__(name)

        # check data types
        pybropt.util.check_is_integer(k, "k")
        pybropt.util.check_is_SetSearchSpace(search_space, "search_space")

        # set private variables
        self._k = k
        self._search_space = search_space

    ############################################################################
    ################################ Properties ################################
    ############################################################################
    def k():
        doc = "The k property."
        def fget(self):
            return self._k
        def fset(self, value):
            self._k = value
        def fdel(self):
            del self._k
        return locals()
    k = property(**k())

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
    def optimize(self, objfn, seed = None, nthreads = 0, minimize = True, verbose = False, *args, **kwargs):
        """
        Perform a hillclimbing using steepest ascent strategy. A set of size
        'k' is used. Single states within the set are iteratively exchanged
        with one state not within the set.

        Parameters
        ==========
        objfn : function
            Objective function to optimize.
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
            print("Steepest Ascent Set HillClimber engine state:")
            print()
            print("Objective function")
            print("==================================================")
            print("Objective function =", objfn.__name__)
            print("Objective function arguments:", str(args))
            print("Objective function keyword arguments:", str(kwargs))
            print()
            print("HillClimber configuration")
            print("==================================================")
            print("Subset size =", self._k)
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

        # get states as numpy.ndarray
        states = self._search_space.state_ndarray

        # iteration counter
        algoiter = 0

        tmp = states.copy()         # copy state list
        numpy.random.shuffle(tmp)   # shuffle the copy
        gbest_pos = tmp[:k].copy()  # extract start position
        gbest_exch = tmp[k:].copy() # extract start exchange vector

        # score the initial configuration
        gbest_score = objfn(gbest_pos, *args, **kwargs)

        # add history to self
        self.history_add(algoiter, [gbest_score], gbest_pos)

        # variables for position, exchange, score matrices
        pos = None
        exch = None
        score = None

        # print verbose statements
        if verbose:
            print("Iteration", algoiter, "\tscore =", score)

        # OPTIMIZE: flip this test variable to not need 'not' in front
        # make exit variable
        found_local_optima = False

        # begin exchange process
        while not found_local_optima:
            # increment iteration
            algoiter += 1

            # create state exchange matrix (local search positions)
            pos, exch = SetHC._set_exchange_matrix(gbest_pos, gbest_exch)

            # score positions that have been created (threading applicable)
            if nthreads > 1:
                # https://stackoverflow.com/questions/3033952/threading-pool-similar-to-the-multiprocessing-pool
                # https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python
                raise RuntimeError("Multithreading not supported yet...")
            else:
                # evaluate the objective function
                score = numpy.apply_along_axis(objfn, 1, pos, *args, **kwargs)

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
                gbest_exch = exch[gbest_ix,:].copy() # # OPTIMIZE: see if copy is needed
                gbest_score = score[gbest_ix]
            else:
                found_local_optima = True

            # add history
            self.history_add(algoiter, score, pos)

            if verbose:
                print("Iteration", algoiter, "\tscore =", gbest_score)

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def _set_exchange_matrix(vec, exch):
        """
        Build a matrix of positions
        """
        out_vec = numpy.tile(vec, (len(vec)*len(exch),1))
        for i in range(len(vec)):
            out_vec[i*len(exch):(i+1)*len(exch),i] = exch
        out_exch = numpy.tile(exch, (len(vec)*len(exch),1))
        for i in range(len(exch)):
            out_exch[i::len(exch),i] = vec
        return out_vec, out_exch
