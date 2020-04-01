# 3rd party libraries
import numpy

# our libraries
from . import SearchSpace
import pybropt.util

class CategoricalSearchSpace(SearchSpace):
    """docstring for CategoricalSearchSpace."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, *args, state = None, dim_size = None, name = "Categorical"):
        """
        Constructor for ContinuousSearchSpace class.

        Parameters
        ----------
        state : numpy.ndarray
            An array of search space potential states.
        dim_size : numpy.ndarray
            An array of array sizes corresponding to the number of states in
            the dimension.
        *args : iterable
            Additional iterable states.
        """
        # call super constructor
        super(CategoricalSearchSpace, self).__init__(name)

        # initialize internal variables
        self.reset()

        # test if we have states and
        if all(a is not None for a in (state, dim_size)):
            # check that data are matrices
            pybropt.util.check_is_matrix(state, "state")
            pybropt.util.check_is_matrix(dim_size, "dim_size")

            # matrix start, stop indices
            stop = dim_size.cumsum()
            start = stop - dim_size

            # loop through states and add them to self.state
            for st,sp in zip(start, stop):
                # convert states to list and add to space
                self._space.append(list(state[st:sp]))
        elif any(a is not None for a in (state, dim_size)):
            raise ValueError("Both 'state' and 'dim_size' needed")

        # loop through variable arguments and add them to self.space
        for i,arg in enumerate(args):
            # check that we can make a list out of the argument
            pybropt.util.check_is_iterable(arg, "arg %d" % (i+1))
            self._space.append(list(arg))

    def __len__(self):
        return len(self._space)

    def __str__(self):
        return str(self._space)

    def __getitem__(self, key):
        """
        Simple indexing.
        key : int, slice
        """
        return self._space[key]

    def __setitem__(self, key, value):
        """
        Simple indexing.
        key : int, slice
        """
        self._space[key] = value

    def __delitem__(self, key):
        """
        Simple indexing.
        key : int, slice
        """
        del self._space[key]

    # TODO: implement me
    # @classmethod
    # def __getitem__(self, key):
    #     """
    #     key : int, slice, tuple, list, numpy.ndarray
    #     """
    #     if isinstance(key, (int, slice)):
    #         return self._space[key]
    #     elif isinstance(key, tuple):
    #         if len(key) == 0:
    #             raise ValueError("Cannot index using empty tuple.")
    #         elif len(key) == 1:
    #             k = *key
    #             if isinstance(k, (int, slice)):
    #                 return self._space[k]
    #             elif isinstance(k, (tuple, list, numpy.ndarray)):
    #                 return [self._space[e] for e in k]
    #             else:
    #                 raise ValueError("Indexing type not supported.")
    #         elif len(key) == 2:
    #             a, b = *key
    #             # get subset list
    #             subset = None
    #             if isinstance(a, (int, slice)):
    #                 subset = self._space[a]
    #             elif isinstance(a, (tuple, list, numpy.ndarray)):
    #                 subset = [self._space[e] for e in a]
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
    #         return [self._space[e] for e in key]
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
    #     del self._space[key]

    ############################################################################
    ################################ Properties ################################
    ############################################################################
    def space():
        doc = "The space property."
        def fget(self):
            return self._space
        def fset(self, value):
            self._space = value
        def fdel(self):
            del self._space
        return locals()
    space = property(**space())

    def ndim():
        doc = "The ndim property."
        def fget(self):
            return len(self._space)
        def fset(self, value):
            error_readonly("ndim")
        def fdel(self):
            error_readonly("ndim")
        return locals()
    ndim = property(**ndim())

    def size():
        doc = "The size property."
        def fget(self):
            return sum(len(dim) for dim in self._space)
        def fset(self, value):
            error_readonly("size")
        def fdel(self):
            error_readonly("size")
        return locals()
    size = property(**size())

    def state():
        doc = "The state property."
        def fget(self):
            return numpy.array([ee for e in self._space for ee in e])
        def fset(self, value):
            error_readonly("state")
        def fdel(self):
            error_readonly("state")
        return locals()
    state = property(**state())

    def dim_size():
        doc = "The dim_size property."
        def fget(self):
            return numpy.array([len(e) for e in self._space])
        def fset(self, value):
            error_readonly("dim_size")
        def fdel(self):
            error_readonly("dim_size")
        return locals()
    dim_size = property(**dim_size())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def reset(self):
        self._space = []

    def feasible(self, pos):
        """
        pos : numpy.ndarray
            A position vector.
        """
        if len(pos) != self.__len__():
            raise ValueError("number of dimensions in 'pos' do not align with search space")

        # iterate over all elements in 'pos' and check if it is the dimension
        f = all(e in self._space[i] for i,e in enumerate(pos))

        # return if all are in the defined space
        return f
