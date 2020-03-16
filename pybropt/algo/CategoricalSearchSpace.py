class CategoricalSearchSpace(SearchSpace):
    """docstring for CategoricalSearchSpace."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################

    @classmethod
    def __init__(self, states = None, dim_sizes = None, *args):
        """
        Constructor for ContinuousSearchSpace class.

        Parameters
        ----------
        states : numpy.ndarray
            An array of search space potential states.
        dim_sizes : numpy.ndarray
            An array of array sizes corresponding to the number of states in
            the dimension.
        *args : iterable
            Additional iterable states.
        """
        # call super constructor
        super(CategoricalSearchSpace, self).__init__()

        # initialize internal variables
        self.reset()

        # test if we have states and
        if (states is not None) and (dim_sizes is not None):
            # check that data are matrices
            check_is_matrix(states, "states")
            check_is_matrix(dim_sizes, "dim_sizes")
            # loop through states and add them to self.state
            stop = dim_sizes.cumsum()
            start = stop - dim_sizes
            for st,sp in zip(start, stop):
                self._state.append(list(states[st:sp]))

        # loop through variable arguments and add them to self.state
        for i,arg in enumerate(args):
            # check that we can make a list out of the argument
            check_is_iterable(arg, "arg %d" % (i+1))
            self._state.append(list(arg))

    @classmethod
    def __len__(self):
        return len(self._state)

    @classmethod
    def __str__(self):
        return str(self._state)

    @classmethod
    def __getitem__(self, key):
        """
        Simple indexing.
        key : int, slice
        """
        return self._state[key]

    @classmethod
    def __setitem__(self, key, value):
        """
        Simple indexing.
        key : int, slice
        """
        self._state[key] = value

    @classmethod
    def __delitem__(self, key):
        """
        Simple indexing.
        key : int, slice
        """
        del self._state[key]

    # TODO: implement me
    # @classmethod
    # def __getitem__(self, key):
    #     """
    #     key : int, slice, tuple, list, numpy.ndarray
    #     """
    #     if isinstance(key, (int, slice)):
    #         return self._state[key]
    #     elif isinstance(key, tuple):
    #         if len(key) == 0:
    #             raise ValueError("Cannot index using empty tuple.")
    #         elif len(key) == 1:
    #             k = *key
    #             if isinstance(k, (int, slice)):
    #                 return self._state[k]
    #             elif isinstance(k, (tuple, list, numpy.ndarray)):
    #                 return [self._state[e] for e in k]
    #             else:
    #                 raise ValueError("Indexing type not supported.")
    #         elif len(key) == 2:
    #             a, b = *key
    #             # get subset list
    #             subset = None
    #             if isinstance(a, (int, slice)):
    #                 subset = self._state[a]
    #             elif isinstance(a, (tuple, list, numpy.ndarray)):
    #                 subset = [self._state[e] for e in a]
    #             else:
    #                 raise ValueError("Indexing type not supported.")
    #             # get values from subset list
    #             if isinstance(b, (int, slice)):
    #                 return [e[b] for e in subset]
    #             elif isinstance(b, (slice, list, numpy.ndarray)):
    #                 return [[e[i] for i in b] for e in subset]
    #             else:
    #                 raise ValueError("Indexing type not supported.")
    #         else:
    #             raise ValueError("tuple length must be 1 or 2.")
    #     elif isinstance(key, (list, numpy.ndarray)):
    #         return [self._state[e] for e in key]
    #     else:
    #         raise ValueError("Indexing type not supported.")
    # TODO: implement me
    # @classmethod
    # def __setitem__(self, key, value):
    #     """
    #     key : int, slice, tuple, list, numpy.ndarray
    #     value : int, tuple, list, numpy.ndarray
    #     """
    #     raise ValueError("Indexing type not supported.")
    # TODO: implement me
    # def __delitem__(self, key):
    #     """
    #     key : int, slice, tuple, list, numpy.ndarray
    #     """
    #     del self._state[key]

    ############################################################################
    ################################ Properties ################################
    ############################################################################

    def state():
        doc = "The state property."
        def fget(self):
            return self._state
        def fset(self, value):
            self._state = value
        def fdel(self):
            del self._state
        return locals()
    state = property(**state())

    def ndim():
        doc = "The ndim property."
        def fget(self):
            return len(self._state)
        def fset(self, value):
            error_readonly("ndim")
        def fdel(self):
            error_readonly("ndim")
        return locals()
    ndim = property(**ndim())

    def size():
        doc = "The size property."
        def fget(self):
            return sum(len(dim) for dim in self._state)
        def fset(self, value):
            error_readonly("size")
        def fdel(self):
            error_readonly("size")
        return locals()
    size = property(**size())

    def dim_size():
        doc = "The dim_size property."
        def fget(self):
            return numpy.array([len(e) for e in self._state])
        def fset(self, value):
            error_readonly("dim_size")
        def fdel(self):
            error_readonly("dim_size")
        return locals()
    dim_size = property(**dim_size())

    def state_flat():
        doc = "The state_flat property."
        def fget(self):
            return numpy.array([ee for ee in e for e in self._state])
        def fset(self, value):
            error_readonly("state_flat")
        def fdel(self):
            error_readonly("state_flat")
        return locals()
    state_flat = property(**state_flat())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    @classmethod
    def reset(self):
        self._state = []
