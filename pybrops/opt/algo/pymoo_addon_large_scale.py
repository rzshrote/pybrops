"""
Module containing PyMOO addons for large-scale optimization algorithms.
"""

__all__ = [
    "TiledSubsetPoolIterator",
    "TiledSubsetRandomSampling",
    "TiledReducedExchangeMutation",
]

from numbers import Integral
from typing import Iterable, Union
import numpy
from pymoo.core.sampling import Sampling
from pymoo.core.problem import Problem
from pymoo.core.mutation import Mutation

from pybrops.core.error.error_type_numpy import check_is_Integral_or_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_ndim

class TiledSubsetPoolIterator:
    """
    TiledSubsetPoolIterator.
    """
    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self,
            setspace: Union[Integral,numpy.ndarray],
            unique: bool = True,
        ) -> None:
        """
        Constructor for TiledSubsetPoolIterator.
        
        Parameters
        ----------
        setspace : Integral, numpy.ndarray
            Set space from which to sample.
            If ``Integral``, sample from a range.
            If ``numpy.ndarray``, sample from an array of possible values.
        unique : bool
            Whether values in tiled samples using ``getnext`` should be unique.
        """
        # process setspace
        if isinstance(setspace, Integral):
            setspace = numpy.arange(setspace)
        elif isinstance(setspace, numpy.ndarray):
            check_ndarray_ndim(setspace, "setspace", 1)
        else:
            check_is_Integral_or_ndarray(setspace, "setspace")
        
        # set private variables
        self._setspace = setspace
        self._pool = setspace.copy()
        self._poolix = 0
        self._unique = unique

        # shuffle/reset pool
        self.reset_pool()

    ################## Iterator interface ##################
    def __iter__(self) -> Iterable:
        """
        Get an interator object from self.

        Returns
        -------
        out : Iterable
            An iterable object.
        """
        return self

    def __next__(self) -> object:
        """
        Get the next value from self.

        Returns
        -------
        out : object
            An object sampled from the tiled set space pool.
        """
        # make sure counter hasn't exceeded length
        if self._poolix >= len(self._pool):
            self.reset_pool()
        
        # get the next element in the pool
        out = self._pool[self._poolix]

        # increment counter
        self._poolix += 1

        # return output
        return out

    ############################ Object Properties #############################

    @property
    def setspace(self) -> numpy.ndarray:
        """Set space from which to sample elements."""
        return self._setspace

    ############################## Object Methods ##############################

    def reset_pool(self) -> None:
        """
        Reset the pool for the TiledSubsetPoolIterator object.
        """
        # shuffle pool
        numpy.random.shuffle(self._pool)

        # reset pool counter index
        self._poolix = 0

    def getnext(self, n: int) -> numpy.ndarray:
        """
        Get the next ``n`` values.

        Parameters
        ----------
        n : int
            Number of values to extract from the iterator.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n,)`` with the next `n` values.
        """
        # make sure counter hasn't exceeded length
        if self._poolix >= len(self._pool):
            self.reset_pool()
        
        # check inputs
        if self._unique and (n > len(self._pool)):
            raise ValueError("cannot be unique: number of requested values exceeds set space")

        # if we want unique values in each sampling
        if self._unique:
            # calculate number of values currently shuffled and available
            nvalues = min(n, len(self._pool) - self._poolix)

            # if we don't have enough shuffled and available values, then perform tiling
            # this is unbiased with respect
            if nvalues < n:
                out1 = self._pool[(self._poolix):(self._poolix+nvalues)].copy()
                out2 = self._setspace[~numpy.in1d(self._setspace,out1)]
                numpy.random.shuffle(out2)
                self._pool[:nvalues] = out1
                self._pool[nvalues:] = out2
                self._poolix = 0

            # extract next ``n`` values in array
            out = self._pool[(self._poolix):(self._poolix+n)].copy()

            # increment counter pointer
            self._poolix += n

            return out
        else:
            # sample without regard to uniqueness
            out = numpy.array([next(self) for _ in range(n)], dtype=self.setspace.dtype)

            return out

class TiledSubsetRandomSampling(
        Sampling,
    ):
    """
    TiledSubsetRandomSampling.
    """
    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self,
            tspiter: TiledSubsetPoolIterator,
        ) -> None:
        """
        Constructor for TiledSubsetRandomSampling.
        
        Parameters
        ----------
        tspiter : TiledSubsetPoolIterator
            A tiled pool iterator object from which to draw samples.
        """
        # call super constructor
        super(TiledSubsetRandomSampling, self).__init__()
        self.tspiter = tspiter

    ############################ Object Properties #############################
    @property
    def tspiter(self) -> TiledSubsetPoolIterator:
        """tspiter."""
        return self._tspiter
    @tspiter.setter
    def tspiter(self, value: TiledSubsetPoolIterator) -> None:
        """Set tspiter."""
        if not isinstance(value, TiledSubsetPoolIterator):
            raise TypeError("``tspiter`` must be of type ``TiledSubsetPoolIterator``")
        self._tspiter = value

    ############################## Object Methods ##############################

    ########### Abstract method implementations ############
    def _do(
            self,
            problem: Problem,
            n_samples: int,
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform chromosome sampling for a provided problem.

        Parameters
        ----------
        problem : Problem
            A problem for which to construct solution chromosomes.
        n_samples : int
            Number of chromosomes which to construct.
        kwargs : dict
            Additional keywork arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A chromosome sampling.
        """
        # sample subset samples
        out = numpy.array(
            [self.tspiter.getnext(problem.n_var) for _ in range(n_samples)],
            dtype=self.tspiter.setspace.dtype
        )

        # return output array
        return out

class TiledReducedExchangeMutation(
        Mutation,
    ):
    """
    Perform a tiled version of reduced exchange mutation. Given a tiled sampler,
    this ensures that all new alleles are tested before repeats are observed.

    This mutation operator is inspired by Correa et al. (2001).

    Citations
    ---------
    Correa, E. S., Steiner, M. T. A., Freitas, A. A., Carnieri, C. A Genetic 
    Algorithm for the P-median Problem In: Genetic and Evolutionary Computation 
    Conference - GECCO 2001, 2001, San Francisco, California. Proceedings of 
    the Genetic and Evolutionary Computation Conference GECCO 2001. San 
    Francisco, California: Morgan Kaufmann Publishers, 2001.
    """

    ################ Special object methods ################
    def __init__(
            self,
            tspiter: TiledSubsetPoolIterator,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        tspiter : TiledSubsetPoolIterator
            A tiled pool iterator object from which to draw samples.
        kwargs : dict
            Additional keyword arguments.
        """
        super(TiledReducedExchangeMutation, self).__init__(**kwargs)
        self.tspiter = tspiter

    ###################### Properties ######################
    @property
    def tspiter(self) -> TiledSubsetPoolIterator:
        """tspiter."""
        return self._tspiter
    @tspiter.setter
    def tspiter(self, value: TiledSubsetPoolIterator) -> None:
        """Set tspiter."""
        if not isinstance(value, TiledSubsetPoolIterator):
            raise TypeError("``tspiter`` must be of type ``TiledSubsetPoolIterator``")
        self._tspiter = value

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: numpy.ndarray, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform exchange mutation for subsets.

        Parameters
        ----------
        problem : Problem
            An optimization problem.
        X : numpy.ndarray
            An array of shape ``(n_indiv, n_var)`` containing individuals which to mutate.
            This array is unmodified by this function.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n_parents, n_matings, n_var)`` containing progeny chromosomes.
        """
        # get number of individuals, variables from X array
        # (n,p)
        n_indiv, n_var = X.shape

        # copy the individual chromosomes
        # (n,p)
        Xm = X.copy()

        # get the probability that a locus is mutated along the chromosome
        p = self.get_prob_var(problem)
        p = numpy.array(p) if hasattr(p, "__len__") else numpy.repeat(p, n_var)

        # for each individual
        for i in range(n_indiv):                    # for each individual
            for j in range(n_var):                  # for each locus
                if numpy.random.random() < p[j]:    # if mutation at locus
                    indiv = Xm[i,:]                 # get view of individual
                    value = next(self.tspiter)      # sample new allele
                    while value in indiv:           # while new allele is in individual
                        value = next(self.tspiter)  # sample new allele
                    Xm[i,j] = value                 # replace old allele with new allele

        return Xm
