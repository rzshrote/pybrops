"""
Module containing addons for PyMOO needed in optimization algorithms.
"""

__all__ = [
    "SubsetRandomSampling",
    "ReducedExchangeCrossover",
    "ReducedExchangeMutation",
    "IntegerPolynomialMutation",
    "IntegerSimulatedBinaryCrossover"
]

import numpy as np
from pymoo.core.sampling import Sampling
from pymoo.core.problem import Problem
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.operators.crossover.sbx import SimulatedBinaryCrossover
from pymoo.operators.mutation.pm import PolynomialMutation

class SubsetRandomSampling(Sampling):
    """
    Class implementing subset sampling for chromosome construction.
    """

    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            replace: bool = False
        ) -> None:
        """
        Constructor for subset random sampling.

        Parameters
        ----------
        setspace : numpy.ndarray
            Set space from which to sample elements.
        replace : numpy.ndarray
            Whether to replace elements when subset sampling.
        """
        # call super constructor
        super(SubsetRandomSampling, self).__init__()

        # perform error checks and make value assignments.
        self.setspace = setspace
        self.replace = replace
    
    ###################### Properties ######################
    @property
    def setspace(self) -> np.ndarray:
        """Set space from which to sample elements."""
        return self._setspace
    @setspace.setter
    def setspace(self, value: np.ndarray) -> None:
        """Set the set space from which to sample elements."""
        if not isinstance(value, np.ndarray):
            raise TypeError("'setspace' must be of type numpy.ndarray")
        self._setspace = value
    
    @property
    def replace(self) -> bool:
        """Whether to replace elements when subset sampling."""
        return self._replace
    @replace.setter
    def replace(self, value: bool) -> None:
        """Set whether to replace elements when subset sampling."""
        if not isinstance(value, bool):
            raise TypeError("'replace' must be of type bool")
        self._replace = value
    
    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            n_samples: int, 
            **kwargs: dict
        ) -> np.ndarray:
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
        # allocate empty output array
        out = np.empty(
            (n_samples, problem.n_var),     # array shape
            dtype = self._setspace.dtype    # set dtype as same as setspace
        )

        # fill array with subset samples
        for i in range(n_samples):
            out[i,:] = np.random.choice(
                self._setspace,             # set space from which to sample
                problem.n_var,              # number of elements to sample 
                replace = self._replace     # whether to replace elements for the subset sampling
            )

        # return output array
        return out

class ReducedExchangeCrossover(Crossover):
    """
    Perform a subset crossover according to Correa et al. (2001).

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
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetCrossover.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(ReducedExchangeCrossover, self).__init__(
            n_parents = 2,
            n_offsprings = 2,
            prob = 0.9,
            **kwargs
        )

    ###################### Properties ######################

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform reduced exchange crossover for subsets.

        Parameters
        ----------
        problem : Problem
            An optimization problem.
        X : numpy.ndarray
            An array of shape ``(n_parents, n_matings, n_var)`` containing parents which to mate.
            This array is unmodified by this function.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n_parents, n_matings, n_var)`` containing progeny chromosomes.
        """
        # get number of parents, matings, variables from X array
        # assume there are only two parents
        n_parents, n_matings, n_var = X.shape

        # copy the parent chromosomes to create progeny chromosomes
        Xp = np.copy(X)

        # for each mate pair, perform exchanges
        for i in range(n_matings):
            mab = ~np.isin(Xp[0,i,:],Xp[1,i,:]) # get mask for indiv 1 not in indiv 2
            mba = ~np.isin(Xp[1,i,:],Xp[0,i,:]) # get mask for indiv 2 not in indiv 1
            ap = Xp[0,i,mab]                    # get reduced indiv 1 chromosome
            bp = Xp[1,i,mba]                    # get reduced indiv 2 chromosome
            clen = min(len(ap), len(bp))        # get minimum chromosome length
            nex = 0 if clen < 2 else np.random.randint(1, clen) # get the number of element exchanges
            mex = np.random.choice(clen, nex)   # get a mask of the elements to exchange
            ap[mex], bp[mex] = bp[mex], ap[mex] # exchange alleles
            Xp[0,i,mab] = ap                    # copy over exchanges to indiv 1
            Xp[1,i,mba] = bp                    # copy over exchanges to indiv 2

        return Xp

class ReducedExchangeMutation(Mutation):
    """
    Perform a subset exchange mutation according to Correa et al. (2001).

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
            setspace: np.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(ReducedExchangeMutation, self).__init__(**kwargs)
        self.setspace = setspace

    ###################### Properties ######################
    @property
    def setspace(self) -> np.ndarray:
        """Set space from which to sample elements."""
        return self._setspace
    @setspace.setter
    def setspace(self, value: np.ndarray) -> None:
        """Set the set space from which to sample elements."""
        if not isinstance(value, np.ndarray):
            raise TypeError("'setspace' must be of type numpy.ndarray")
        self._setspace = value
    
    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
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
        n_indiv, n_var = X.shape

        # copy the individual chromosomes
        Xm = X.copy()

        # get the probability that a locus is mutated along the chromosome
        p = self.get_prob_var(problem)
        p = np.array(p) if hasattr(p, "__len__") else np.repeat(p, n_var)

        # for each individual
        for i in range(n_indiv):
            mab = ~np.isin(Xm[i,:],self.setspace)   # get mask for indiv not in set space
            mba = ~np.isin(self.setspace,Xm[i,:])   # get mask for set space not in indiv
            ap = Xm[i,mab]                          # get reduced indiv chromosome
            pp = p[mab]                             # get reduced indiv mutation probability
            bp = self.setspace[mba]                 # get reduced setspace chromosome
            mex = np.random.random(len(pp)) < pp    # get mutation locations
            nex = mex.sum()                         # get number of mutations
            ap[mex] = np.random.choice(bp, nex)     # randomly assign new alleles to reduced indiv chromosome
            Xm[i,mab] = ap                          # overwrite mutations to indiv
        
        return Xm

class IntegerSimulatedBinaryCrossover(SimulatedBinaryCrossover):
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform simulated binary crossover, but round results to nearest integer.

        Parameters
        ----------
        problem : Problem
            An optimization problem.
        X : numpy.ndarray
            An array of shape ``(n_parents, n_matings, n_var)`` containing parents which to mate.
            This array is unmodified by this function.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n_parents, n_matings, n_var)`` containing progeny chromosomes.
        """
        # perform SBX
        out = super(IntegerSimulatedBinaryCrossover, self)._do(
            problem = problem,
            X = X,
            **kwargs
        )

        # round results and convert type
        out = out.round(0).astype(X.dtype)

        return out

class IntegerPolynomialMutation(PolynomialMutation):
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform polynomial muation, but round results to nearest integer.

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
        # perform PM
        out = super(IntegerPolynomialMutation, self)._do(
            problem = problem,
            X = X,
            **kwargs
        )

        # round results and convert type
        out = out.round(0).astype(X.dtype)

        return out
