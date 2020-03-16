class SetSearchSpace(SearchSpace):
    """docstring for SetSearchSpace."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################

    def __init__(self, states = None, *args):
        """
        Constructor for SetSearchSpace class.

        Parameters
        ----------
        states : list-like
            An list of search space potential states.
        *args : iterable
            Additional iterable states.
        """
        # call super constructor
        super(SetSearchSpace, self).__init__()

        # initialize internal variables
        self.reset()

        if states is not None:
            for e in states:
                self._state.append(e)

        for arg in args:
            self._state.append(arg)

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
            return len(self._state)
        def fset(self, value):
            error_readonly("size")
        def fdel(self):
            error_readonly("size")
        return locals()
    size = property(**size())

    def state_ndarray():
        doc = "The state_ndarray property."
        def fget(self):
            return numpy.array(self._state)
        def fset(self, value):
            error_readonly("state_ndarray")
        def fdel(self):
            error_readonly("state_ndarray")
        return locals()
    state_ndarray = property(**state_ndarray())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    @classmethod
    def reset(self):
        self._state = []
