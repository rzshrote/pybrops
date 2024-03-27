"""
Module containing addons for PyMOO needed in optimization algorithms.
"""

__all__ = [
    "dominates",
    "SubsetRandomSampling",
    "ReducedExchangeCrossover",
    "ReducedExchangeMutation",
    "MultiObjectiveStochasticHillClimberMutation",
    "IntegerPolynomialMutation",
    "IntegerSimulatedBinaryCrossover",
    "MultiObjectiveSteepestDescentHillClimberMutation",
    "MultiObjectiveStochasticDescentHillClimberMutation",
]

from typing import Optional
import numpy as np
import copy
from pymoo.core.variable import get
from pymoo.core.sampling import Sampling
from pymoo.core.problem import Problem
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.population import Population
from pymoo.core.individual import Individual
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from pymoo.util.dominator import get_relation
from pymoo.operators.crossover.sbx import SimulatedBinaryCrossover
from pymoo.operators.mutation.pm import PolynomialMutation

def dominates(obj1: np.ndarray, cv1: float, obj2: np.ndarray, cv2: float) -> bool:
    """
    Test whether solution 1 dominates solution 2 using a modification for constraints.
    Preference is given for solutions with the least constraint violation.

    Parameters
    ----------
    obj1 : numpy.ndarray
        A vector of objective function evaluations for solution 1.
    cv1 : float
        The constaint violation score for solution 1.
        A score of <= 0.0 is assumed to be no constraint violation.
    obj2 : numpy.ndarray
        A vector of objective function evaluations for solution 2.
    cv2 : float
        The constaint violation score for solution 2.
        A score of <= 0.0 is assumed to be no constraint violation.
    
    Returns
    -------
    out : bool
        Whether solution 1 dominates solution 2 (True), or whether solution 1 
        is non-dominated or dominated by solution 2 (False).
    """
    # if neither solution has any constraint violations (constraints are less than zero)
    if cv1 <= 0.0 and cv2 <= 0.0:
        # if soln 1 <= soln 2 in all objectives and soln 1 < soln 2 in any objective, then soln 1 dominates soln 2
        return np.all(obj1 <= obj2) and np.any(obj1 < obj2)
    # otherwise return whether solution 1 has less constraint violations than solution 2
    # if cv1 < cv2, then soln 1 dominates soln 2
    # if cv1 >= cv2, then soln 1 may be dominated or non-dominated by soln 2
    return cv1 < cv2

def tiled_choice(a: int, size: int):
    out = np.empty(size, int)
    ndiv = size // a
    nrem = size % a
    for i in range(ndiv):
        out[a*i:a*(i+1)] = np.random.choice(a, a, replace = False)
    out[a*ndiv:] = np.random.choice(a, nrem, replace = False)
    return out

class SubsetRandomSampling(
        Sampling,
    ):
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

class ReducedExchangeCrossover(
        Crossover,
    ):
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

class ReducedExchangeMutation(
        Mutation,
    ):
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

class IntegerSimulatedBinaryCrossover(
        SimulatedBinaryCrossover,
    ):
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

class IntegerPolynomialMutation(
        PolynomialMutation,
    ):
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

class MultiObjectiveStochasticHillClimberMutation(
        Mutation,
    ):
    """
    Perform a memetic subset exchange mutation.
    """

    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            p_hillclimb: float,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(MultiObjectiveStochasticHillClimberMutation, self).__init__(**kwargs)
        self.setspace = setspace
        self.p_hillclimb = p_hillclimb

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
    
    ### Helper methods ###
    def hillclimb(self, problem: Problem, x: np.ndarray, *args, **kwargs):
        """
        Hillclimb mutation for a single individual.

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            Problem object.
        x : numpy.ndarray
            Chromosome of a single individual.
        
        Returns
        -------
        out : numpy.ndarray
            Chromosome of a single individual.
        """
        # initialize the leader solution
        lead = x.copy()

        # get search space elements not in current leader solution
        wrkss = self.setspace[np.logical_not(np.in1d(self.setspace, lead))]

        # evaluate the lead solution
        # use 2d array view 
        lead_eval = {}
        problem._evaluate(lead[None,:], lead_eval, *args, **kwargs)

        # unpack evaluation for leader solution and calculate constraint violation
        lead_obj    = lead_eval["F"][0] if "F" in lead_eval else np.zeros(0, float)
        lead_ineqcv = lead_eval["G"][0] if "G" in lead_eval else np.zeros(0, float)
        lead_eqcv   = lead_eval["H"][0] if "H" in lead_eval else np.zeros(0, float)
        lead_cv     = lead_ineqcv.sum() + lead_eqcv.sum()

        # create lists to house non-dominated solutions
        ndom     = []
        ndom_obj = []
        ndom_cv  = []

        # for each element in the solution vector
        for i in range(len(lead)):
            # randomly choose an index to exchange
            j = np.random.choice(len(wrkss))
            # exchange indices to propose new solution
            lead[i], wrkss[j] = wrkss[j], lead[i]
            # evaluate proposed solution
            prop_eval = {}
            problem._evaluate(lead[None,:], prop_eval, *args, **kwargs)
            # unpack evaluation for proposed solution and calculate constraint violation
            prop_obj    = prop_eval["F"][0] if "F" in prop_eval else np.zeros(0, float)
            prop_ineqcv = prop_eval["G"][0] if "G" in prop_eval else np.zeros(0, float)
            prop_eqcv   = prop_eval["H"][0] if "H" in prop_eval else np.zeros(0, float)
            prop_cv     = prop_ineqcv.sum() + prop_eqcv.sum()
            
            # if the leader solution dominates the proposed solution, then 
            # reject the proposed solution and continue to the next locus.
            if dominates(lead_obj, lead_cv, prop_obj, prop_cv):
                lead[i], wrkss[j] = wrkss[j], lead[i]   # exchange values back to original solution
            
            # if the leader solution is dominated by the proposed solution,
            # then overwrite the leader with the proposed and check all
            # non-dominated solutions against the new leader.
            elif dominates(prop_obj, prop_cv, lead_obj, lead_cv):
                lead_obj    = prop_obj                  # overwrite objective values
                lead_cv     = prop_cv                   # overwrite constraint violation value
                k = 0
                while k < len(ndom):
                    # if new leader dominates non-dominated in list
                    if dominates(lead_obj, lead_cv, ndom_obj[k], ndom_cv[k]):
                        del ndom[k]                     # remove non-dominated chromosome value
                        del ndom_obj[k]                 # remove non-dominated objective values
                        del ndom_cv[k]                  # remove non-dominated constraint violation
                    else:
                        k += 1
            
            # otherwise the leader is non-dominated by the proposed solution, so
            # stash non-dominated solution into solutions lists and return to 
            # leader sequence
            else:
                ndom.append(lead.copy())                # stash copy of chromosome
                ndom_obj.append(lead_obj)               # stash objective values
                ndom_cv.append(lead_cv)                 # stash constraint violation value
                lead[i], wrkss[j] = wrkss[j], lead[i]   # exchange values back to original solution

        # leader is already non-dominated for all other individuals, so 
        # just remove dominated individuals for all other non-dominated pairs
        # use naive algorithm for removing dominated individiuals
        isndom = [True] * len(ndom)
        for i in range(len(ndom)-1):
            for j in range(i+1,len(ndom)):
                # if i dominates j, set isdom[j] to False, otherwise keep to True
                isndom[j] &= not dominates(ndom_obj[i],ndom_cv[i],ndom_obj[j],ndom_cv[j])
                # if j dominates i, set isdom[j] to False, otherwise keep to True
                isndom[i] &= not dominates(ndom_obj[j],ndom_cv[j],ndom_obj[i],ndom_cv[i])

        # filter out dominated solutions, leaving nondominated solutions
        ndom     = [ndom[i]     for i in range(len(ndom))     if isndom[i]]
        ndom_obj = [ndom_obj[i] for i in range(len(ndom_obj)) if isndom[i]]
        ndom_cv  = [ndom_cv[i]  for i in range(len(ndom_cv))  if isndom[i]]

        # append the leader to the ndom list
        ndom.insert(0, lead.copy())     # stash copy of chromosome
        ndom_obj.insert(0, lead_obj)    # stash objective values
        ndom_cv.insert(0, lead_cv)      # stash constraint violation value

        # determine the solution furthest from the leader for maximum mutation
        ix = 0
        dist = np.linalg.norm(lead_obj - ndom_obj[ix], 2)
        for i in range(1, len(ndom)):
            tmp = np.linalg.norm(lead_obj - ndom_obj[i], 2)
            dist = tmp if tmp > dist else dist

        return ndom[ix]

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform exchange mutation for subsets + hillclimb

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

        # apply hillclimber to mutated individuals
        for i in range(n_indiv):
            if np.random.random() < self.p_hillclimb:
                Xm[i,:] = self.hillclimb(problem, Xm[i,:])

        return Xm

class MultiObjectiveSteepestDescentHillClimberMutation(
        Mutation,
    ):
    """
    Perform a memetic subset exchange mutation.
    """

    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            p_hillclimb: float,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(MultiObjectiveSteepestDescentHillClimberMutation, self).__init__(**kwargs)
        self.setspace = setspace
        self.p_hillclimb = p_hillclimb

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
    
    ### Helper methods ###
    def hillclimb(
            self, 
            problem: Problem, 
            indiv: Individual,
            *args: tuple, 
            **kwargs: dict
        ) -> Population:
        """
        Hillclimb mutation for a single individual.
        For a single individual, pick a single allele locus and mutate it with
        all alternative alleles. Remove individuals that are dominated.
        Return a population of individuals that are non-dominated.

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            Problem object.
        indiv : pymoo.core.individual.Individual
            A single individual to mutate via a hillclimb.
        
        Returns
        -------
        out : pymoo.core.population.Population
            Chromosome of a single individual.
        """
        # get chromosome from individual
        x = indiv.get("X")
        
        # get search space elements not in current individual
        wrkss = self.setspace[np.logical_not(np.in1d(self.setspace, x))]

        # create X matrix of holding mutated individuals
        X = np.empty((len(wrkss),len(x)), dtype = x.dtype)
        X[:,:] = x[None,:]
        
        # randomly select a locus to be mutated
        locus = np.random.choice(len(x))

        # set alleles at selected locus to those not found in the starting chromosome
        X[:,locus] = wrkss
        
        # evaluate mutated individuals
        evals = problem.evaluate(X, *args, return_as_dictionary = True, **kwargs)

        # create a new population from the individuals
        pop = Population.new("X", X)

        # set "F", "G", "H" fields for individuals in population
        for key, value in evals.items():
            pop.set(key, value)

        # get objective function values
        F = pop.get("F")

        # perform non-dominated sorting on populations
        nds = NonDominatedSorting()

        # get non-dominated front indices
        ix = nds.do(F, only_non_dominated_front = True)

        # subset population for individuals along the non-dominated front
        pop = pop[ix]
        
        return pop

    def do(
            self, 
            problem: Problem, 
            pop: Population, 
            inplace: bool = True, 
            **kwargs: dict
        ) -> Population:
        """
        Perform steepest descent hillclimber for subset chromosomes

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            An optimization problem.
        pop : pymoo.core.population.Population
            A population of individuals to mutate.
        inplace : bool
            Whether to modify a population in-place.
        
        Returns
        -------
        out : pymoo.core.population.Population
            A mutated population.
        """
        # if not inplace copy the population first
        if not inplace:
            pop = copy.deepcopy(pop)

        n_mut = len(pop)

        # get the variables to be mutated
        X = pop.get("X")

        # retrieve the mutation variables
        Xp = self._do(problem, X, **kwargs)

        # the likelihood for a mutation on the individuals
        prob = get(self.prob, size=n_mut)
        mut = np.random.random(size=n_mut) <= prob

        # store the mutated individual back to the population
        pop[mut].set("X", Xp[mut])

        # sample from a binomial distribution the number of individuals to hillclimb
        nhc = np.random.binomial(len(pop), self.p_hillclimb)

        # sample indices for individuals to hillclimb
        hcix = np.random.choice(len(pop), nhc, replace = False)

        # apply hillclimber to select individuals
        hcm = [self.hillclimb(problem, pop[ix]) for ix in hcix]

        # merge populations
        out = Population.merge(pop, *hcm) if len(hcm) > 0 else pop
        return out

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform exchange mutation for subsets + hillclimb

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

class MultiObjectiveStochasticDescentHillClimberMutation(
        Mutation,
    ):
    """
    Perform a memetic subset exchange mutation.
    """

    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            phc: float,
            nhc: Optional[int] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for MultiObjectiveStochasticDescentHillClimberMutation.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(MultiObjectiveStochasticDescentHillClimberMutation, self).__init__(**kwargs)
        self.setspace = setspace
        self.phc = phc
        self.nhc = nhc

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
    
    ### Helper methods ###
    def hillclimb(
            self, 
            prob: Problem, 
            indiv: Individual,
            *args: tuple, 
            **kwargs: dict
        ) -> Population:
        """
        Hillclimb mutation for a single individual.
        For a single individual, pick a single allele locus and mutate it with
        all alternative alleles. Remove individuals that are dominated.
        Return a population of individuals that are non-dominated.

        Parameters
        ----------
        prob : pymoo.core.problem.Problem
            Problem object.
        indiv : pymoo.core.individual.Individual
            A single individual to mutate via a hillclimb.
        
        Returns
        -------
        out : pymoo.core.population.Population
            Chromosome of a single individual.
        """
        # copy and evaluate input individual (current leader)
        lead = copy.deepcopy(indiv)
        lead_evals = prob.evaluate(lead.X, *args, return_as_dictionary = True, **kwargs)
        lead.set_by_dict(**lead_evals)

        # get search space elements not in original leader individual
        alleles = self.setspace[np.logical_not(np.in1d(self.setspace, lead.X))]

        # calculate the number of loci and number of alternative alleles
        nloci = len(lead.X)
        nalleles = len(alleles)

        # calculate the number of hillclimb/mutation attempts
        nhc = nloci if self.nhc is None else self.nhc

        # indicator variable for whether leader has been improved
        lead_improved = False

        # calculate the loci to mutate using a tiled sampling method
        lociix = np.empty(nhc, int)
        for i in range(nhc//nloci):
            lociix[nloci*i:nloci*(i+1)] = np.random.choice(nloci, nloci, replace = False)
        lociix[nloci*(nhc//nloci):] = np.random.choice(nloci, nhc % nloci, replace = False)

        # create lists to house non-dominated individuals
        ndom = []

        # for each locus index
        for i in lociix:
            # get index of allele to exchange
            j = np.random.choice(nalleles)

            # copy the leader, modify the copy, evaluate new solution
            prop = Individual()
            prop = copy.deepcopy(lead)
            prop.X[i], alleles[j] = alleles[j], prop.X[i]
            prop_evals = prob.evaluate(prop.X, *args, return_as_dictionary = True, **kwargs)
            prop.set_by_dict(**prop_evals)
            
            # get dominance relationship between two individuals
            rel = get_relation(lead, prop)
            
            # if lead is dominated by proposed solution, keep the proposed as the new leader
            if rel == -1:
                lead_improved = True
                lead = prop
            
            # if the lead is non-dominated by the proposed solution, stash the 
            # proposed and return lead vector back to its original state.
            elif rel == 0:
                ndom.append(prop)
                prop.X[i], alleles[j] = alleles[j], prop.X[i]

            # if the lead dominates the proposed solution, 
            # return the lead vector back to its original state
            elif rel == 1:
                prop.X[i], alleles[j] = alleles[j], prop.X[i]

        # if the leader was improved, add to non-dominated individuals, otherwise keep out
        if lead_improved:
            ndom.append(lead)
        
        # create a new population from non-dominated individual pool
        pop = Population(individuals = ndom)

        # perform non-dominated sorting on populations
        nds = NonDominatedSorting()

        # get non-dominated front indices
        ix = nds.do(pop.get("F"), only_non_dominated_front = True)

        # subset population for individuals along the non-dominated front
        pop = pop[ix]

        return pop

    def do(
            self, 
            problem: Problem, 
            pop: Population, 
            inplace: bool = True, 
            **kwargs: dict
        ) -> Population:
        """
        Perform steepest descent hillclimber for subset chromosomes

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            An optimization problem.
        pop : pymoo.core.population.Population
            A population of individuals to mutate.
        inplace : bool
            Whether to modify a population in-place.
        
        Returns
        -------
        out : pymoo.core.population.Population
            A mutated population.
        """
        print("len(pop)=",len(pop))
        ### Step 1) Perform hillclimb on (presumeably) the best individuals

        # sample from a binomial distribution the number of individuals to hillclimb
        nhc = np.random.binomial(len(pop), self.phc)

        # sample indices for individuals to hillclimb
        hcix = np.random.choice(len(pop), nhc, replace = False)

        # apply hillclimber to select individuals; this is a list of populations
        hcm = [self.hillclimb(problem, pop[ix]) for ix in hcix]

        print("nhc=",nhc, "len(hcm)=",len(hcm))

        # if not inplace copy the population first
        if not inplace:
            pop = copy.deepcopy(pop)

        # mutate population chromosomes (as copy)
        Xp = self._do(problem, pop.get("X"), **kwargs)

        # the likelihood for a mutation on the individuals
        prob = get(self.prob, size = len(pop))
        mut = np.random.random(size = len(pop)) <= prob

        # store select mutated individuals back to the population
        pop[mut].set("X", Xp[mut])

        # merge populations
        out = Population.merge(pop, *hcm) if len(hcm) > 0 else pop
        

        return out

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform exchange mutation for subsets + hillclimb

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

class StochasticHillClimberMutation(
        Mutation,
    ):
    """
    Perform a memetic subset exchange mutation.
    """

    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            phc: float,
            nhcstep: int,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        setspace : numpy.ndarray
            Set space of which individuals are a discrete subset.
        phc : float
            Probability that an individiual is hillclimbed.
        nhc : int
            Number of steps to take during hillclimbing.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(StochasticHillClimberMutation, self).__init__(**kwargs)
        self.setspace = setspace
        self.phc = phc
        self.nhcstep = nhcstep

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
    
    ### Helper methods ###
    def reduced_exchange(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Perform a reduced exchange mutation
        """
        # get the probability that a locus is mutated along the chromosome
        p = self.get_prob_var(problem)
        p = np.array(p) if hasattr(p, "__len__") else np.repeat(p, len(x))

        # mutate the individual
        mab = ~np.isin(x,self.setspace)      # get mask for indiv not in set space
        mba = ~np.isin(self.setspace,x)      # get mask for set space not in indiv
        ap = x[mab]                          # get reduced indiv chromosome
        pp = p[mab]                          # get reduced indiv mutation probability
        bp = self.setspace[mba]              # get reduced setspace chromosome
        mex = np.random.random(len(pp)) < pp # get mutation locations
        nex = mex.sum()                      # get number of mutations
        ap[mex] = np.random.choice(bp, nex)  # randomly assign new alleles to reduced indiv chromosome
        x[mab] = ap                          # overwrite mutations to indiv

        return x

    def hillclimb(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Hillclimb mutation for a single individual.

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            Problem object.
        x : numpy.ndarray
            Chromosome of a single individual.
        
        Returns
        -------
        out : numpy.ndarray
            Chromosome of a single individual.
        """
        # copy and evaluate input individual (current leader)
        lead = Individual()
        lead.X = x.copy()
        lead_evals = problem.evaluate(lead.X, *args, return_as_dictionary = True, **kwargs)
        lead.set_by_dict(**lead_evals)

        # get search space elements not in original leader individual
        alleles = self.setspace[np.logical_not(np.in1d(self.setspace, lead.X))]

        # calculate the number of loci and number of alternative alleles
        nloci = len(lead.X)
        nalleles = len(alleles)

        # calculate the number of hillclimb/mutation attempts
        nhcstep = nloci if self.nhcstep is None else self.nhcstep

        # indicator variable for whether leader has been improved
        lead_improved = False

        # calculate the loci to mutate using a tiled sampling method
        lociix = np.empty(nhcstep, int)
        for i in range(nhcstep//nloci):
            lociix[nloci*i:nloci*(i+1)] = np.random.choice(nloci, nloci, replace = False)
        lociix[nloci*(nhcstep//nloci):] = np.random.choice(nloci, nhcstep % nloci, replace = False)

        # create lists to house non-dominated individuals
        ndom = []

        # for each locus index
        for i in lociix:
            # get index of allele to exchange
            j = np.random.choice(nalleles)

            # copy the leader, modify the copy, evaluate new solution
            prop = Individual()
            prop = copy.deepcopy(lead)
            prop.X[i], alleles[j] = alleles[j], prop.X[i]
            prop_evals = problem.evaluate(prop.X, *args, return_as_dictionary = True, **kwargs)
            prop.set_by_dict(**prop_evals)
            
            # get dominance relationship between two individuals
            rel = get_relation(lead, prop)
            
            # if lead is dominated by proposed solution, keep the proposed as the new leader
            if rel == -1:
                lead_improved = True
                lead = prop
            
            # if the lead is non-dominated by the proposed solution, stash the 
            # proposed and return lead vector back to its original state.
            elif rel == 0:
                ndom.append(prop)
                prop.X[i], alleles[j] = alleles[j], prop.X[i]

            # if the lead dominates the proposed solution, 
            # return the lead vector back to its original state
            elif rel == 1:
                prop.X[i], alleles[j] = alleles[j], prop.X[i]

        # if the leader was improved, add to non-dominated individuals, otherwise keep out
        if lead_improved:
            ndom.append(lead)
        
        # if there is nothing to select from, revert to simple mutation
        if len(ndom) == 0:
            return self.reduced_exchange(problem, x, *args, **kwargs)

        # create a new population from non-dominated individual pool
        pop = Population(individuals = ndom)

        # perform non-dominated sorting on populations
        nds = NonDominatedSorting()

        # get non-dominated front indices
        ix = nds.do(pop.get("F"), only_non_dominated_front = True)

        # subset population for individuals along the non-dominated front
        pop = pop[ix]

        # get F values for population
        F = pop.get("F")

        # identify the indices of the extreme individuals
        Famax = np.argmin(F, axis = 0)

        # randomly select an index of an extreme individual to select
        Fix = np.random.choice(Famax)

        # select that individual
        indiv = pop[Fix]

        # get chromosome of individual
        out = indiv.X

        return out

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform exchange mutation for subsets + hillclimb

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

        for i in range(n_indiv):
            rnd = np.random.random()
            if rnd < self.phc:
                Xm[i,:] = self.hillclimb(problem, Xm[i,:])
            else:
                Xm[i,:] = self.reduced_exchange(problem, Xm[i,:])
        
        return Xm

class MutatorA(
        Mutation,
    ):
    """
    Randomly test nearby solutions and randomly choose one solution that 
    non-dominates from the tested mutations.
    """
    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            phc: float,
            nhcstep: int,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        setspace : numpy.ndarray
            Set space of which individuals are a discrete subset.
        phc : float
            Probability that an individiual is hillclimbed.
        nhc : int
            Number of steps to take during hillclimbing.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(MutatorA, self).__init__(**kwargs)
        self.setspace = setspace
        self.phc = phc
        self.nhcstep = nhcstep

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
    
    ### Helper methods ###
    def reduced_exchange(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Perform a reduced exchange mutation
        """
        # get the probability that a locus is mutated along the chromosome
        p = self.get_prob_var(problem)
        p = np.array(p) if hasattr(p, "__len__") else np.repeat(p, len(x))

        # mutate the individual
        mab = ~np.isin(x,self.setspace)      # get mask for indiv not in set space
        mba = ~np.isin(self.setspace,x)      # get mask for set space not in indiv
        ap = x[mab]                          # get reduced indiv chromosome
        pp = p[mab]                          # get reduced indiv mutation probability
        bp = self.setspace[mba]              # get reduced setspace chromosome
        mex = np.random.random(len(pp)) < pp # get mutation locations
        nex = mex.sum()                      # get number of mutations
        ap[mex] = np.random.choice(bp, nex)  # randomly assign new alleles to reduced indiv chromosome
        x[mab] = ap                          # overwrite mutations to indiv

        return x

    def hillclimb(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Hillclimb mutation for a single individual.

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            Problem object.
        x : numpy.ndarray
            Chromosome of a single individual.
        
        Returns
        -------
        out : numpy.ndarray
            Chromosome of a single individual.
        """
        # copy and evaluate input individual (current leader)
        origin = Individual()
        origin.X = x.copy()
        lead_evals = problem.evaluate(origin.X, *args, return_as_dictionary = True, **kwargs)
        origin.set_by_dict(**lead_evals)

        # get search space elements not in original leader individual
        alleles = self.setspace[np.logical_not(np.in1d(self.setspace, origin.X))]

        # calculate the:
        # 1) Number of loci
        # 2) Number of hillclimb mutation attempts
        # 3) Number of alternative alleles
        nloci = len(origin.X)
        nhcstep = nloci if self.nhcstep is None else self.nhcstep
        nalleles = len(alleles)

        # calculate the loci to mutate using a tiled sampling method
        lociix = tiled_choice(nloci, nhcstep)
        alleleix = tiled_choice(nalleles, nhcstep)

        # create X matrix holding mutations for testing mutations
        Xhc = np.empty((nhcstep,nloci), dtype = x.dtype)
        Xhc[:,:] = x[None,:]

        # apply mutations identify neighbors
        Xhc[:,lociix] = alleles[alleleix]

        # create hillclimber population
        pophc = Population.new(X = Xhc)

        # evaluate individuals
        pophc_eval = problem.evaluate(pophc.get("X"), *args, return_as_dictionary = True, **kwargs)
        for key, value in pophc_eval.items():
            pophc.set(key, value)

        # perform non-dominated sorting on populations
        nds = NonDominatedSorting()

        # get F values for population
        F = pophc.get("F")

        # get non-dominated front indices
        ndix = nds.do(F, only_non_dominated_front = True)

        # subset population for individuals along the non-dominated front
        pophc = pophc[ndix]

        # randomly select an individual from the frontier
        selix = np.random.choice(len(pophc))

        # select that individual
        indiv = pophc[selix]

        # get chromosome of individual
        out = indiv.X

        return out

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform exchange mutation for subsets + hillclimb

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

        for i in range(n_indiv):
            rnd = np.random.random()
            if rnd < self.phc:
                Xm[i,:] = self.hillclimb(problem, Xm[i,:])
            else:
                Xm[i,:] = self.reduced_exchange(problem, Xm[i,:])
        
        return Xm

class MutatorB(
        Mutation,
    ):
    """
    Randomly test nearby solutions and randomly choose one solution that 
    exhibits extreme values for one of the objectives.
    """
    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            phc: float,
            nhcstep: int,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        setspace : numpy.ndarray
            Set space of which individuals are a discrete subset.
        phc : float
            Probability that an individiual is hillclimbed.
        nhc : int
            Number of steps to take during hillclimbing.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(MutatorB, self).__init__(**kwargs)
        self.setspace = setspace
        self.phc = phc
        self.nhcstep = nhcstep

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
    
    ### Helper methods ###
    def reduced_exchange(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Perform a reduced exchange mutation
        """
        # get the probability that a locus is mutated along the chromosome
        p = self.get_prob_var(problem)
        p = np.array(p) if hasattr(p, "__len__") else np.repeat(p, len(x))

        # mutate the individual
        mab = ~np.isin(x,self.setspace)      # get mask for indiv not in set space
        mba = ~np.isin(self.setspace,x)      # get mask for set space not in indiv
        ap = x[mab]                          # get reduced indiv chromosome
        pp = p[mab]                          # get reduced indiv mutation probability
        bp = self.setspace[mba]              # get reduced setspace chromosome
        mex = np.random.random(len(pp)) < pp # get mutation locations
        nex = mex.sum()                      # get number of mutations
        ap[mex] = np.random.choice(bp, nex)  # randomly assign new alleles to reduced indiv chromosome
        x[mab] = ap                          # overwrite mutations to indiv

        return x

    def hillclimb(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Hillclimb mutation for a single individual.

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            Problem object.
        x : numpy.ndarray
            Chromosome of a single individual.
        
        Returns
        -------
        out : numpy.ndarray
            Chromosome of a single individual.
        """
        # copy and evaluate input individual (current leader)
        origin = Individual()
        origin.X = x.copy()
        lead_evals = problem.evaluate(origin.X, *args, return_as_dictionary = True, **kwargs)
        origin.set_by_dict(**lead_evals)

        # get search space elements not in original leader individual
        alleles = self.setspace[np.logical_not(np.in1d(self.setspace, origin.X))]

        # calculate the:
        # 1) Number of loci
        # 2) Number of hillclimb mutation attempts
        # 3) Number of alternative alleles
        nloci = len(origin.X)
        nhcstep = nloci if self.nhcstep is None else self.nhcstep
        nalleles = len(alleles)

        # calculate the loci to mutate using a tiled sampling method
        lociix = tiled_choice(nloci, nhcstep)
        alleleix = tiled_choice(nalleles, nhcstep)

        # create X matrix holding mutations for testing mutations
        Xhc = np.empty((nhcstep,nloci), dtype = x.dtype)
        Xhc[:,:] = x[None,:]

        # apply mutations identify neighbors
        Xhc[:,lociix] = alleles[alleleix]

        # create hillclimber population
        pophc = Population.new(X = Xhc)

        # evaluate individuals
        pophc_eval = problem.evaluate(pophc.get("X"), *args, return_as_dictionary = True, **kwargs)
        for key, value in pophc_eval.items():
            pophc.set(key, value)

        # perform non-dominated sorting on populations
        nds = NonDominatedSorting()

        # get F values for population
        F = pophc.get("F")

        # get non-dominated front indices
        ndix = nds.do(F, only_non_dominated_front = True)

        # subset population for individuals along the non-dominated front
        pophc = pophc[ndix]

        # find the indices of individuals with extreme (minimum) values for each objective
        minix = np.argmin(F, axis = 0)

        # randomly select one individual with an extreme objective
        selix = np.random.choice(minix)

        # select that individual
        indiv = pophc[selix]

        # get chromosome of individual
        out = indiv.X

        return out

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform exchange mutation for subsets + hillclimb

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

        for i in range(n_indiv):
            rnd = np.random.random()
            if rnd < self.phc:
                Xm[i,:] = self.hillclimb(problem, Xm[i,:])
            else:
                Xm[i,:] = self.reduced_exchange(problem, Xm[i,:])
        
        return Xm

class MutatorF(
        Mutation,
    ):
    """
    Perform a memetic subset exchange mutation (Mutation F).
    """

    ################ Special object methods ################
    def __init__(
            self,
            setspace: np.ndarray,
            phc: float,
            maxhc: Optional[int] = None,
            maxhcstep: Optional[int] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMutation.
        
        Parameters
        ----------
        setspace : numpy.ndarray
            Set space of which individuals are a discrete subset.
        phc : float
            Probability that an individiual is hillclimbed.
        maxhc : int
            Maximum number of solutions to test in a hillclimb.
            If None, then set to the number of loci in an individual's chromosome.
        maxhcstep : int
            Maximum number of steps (assignment of new leader solutions) to take.
            If None, then set to the maximum number of solutions tested.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(StochasticHillClimberMutation, self).__init__(**kwargs)
        self.setspace = setspace
        self.phc = phc
        self.maxhc = maxhc
        self.maxhcstep = maxhcstep

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
    
    ### Helper methods ###
    def reduced_exchange(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Perform a reduced exchange mutation
        """
        # get the probability that a locus is mutated along the chromosome
        p = self.get_prob_var(problem)
        p = np.array(p) if hasattr(p, "__len__") else np.repeat(p, len(x))

        # mutate the individual
        mab = ~np.isin(x,self.setspace)      # get mask for indiv not in set space
        mba = ~np.isin(self.setspace,x)      # get mask for set space not in indiv
        ap = x[mab]                          # get reduced indiv chromosome
        pp = p[mab]                          # get reduced indiv mutation probability
        bp = self.setspace[mba]              # get reduced setspace chromosome
        mex = np.random.random(len(pp)) < pp # get mutation locations
        nex = mex.sum()                      # get number of mutations
        ap[mex] = np.random.choice(bp, nex)  # randomly assign new alleles to reduced indiv chromosome
        x[mab] = ap                          # overwrite mutations to indiv

        return x

    def hillclimb(self, problem: Problem, x: np.ndarray, *args, **kwargs) -> np.ndarray:
        """
        Hillclimb mutation for a single individual.

        Parameters
        ----------
        problem : pymoo.core.problem.Problem
            Problem object.
        x : numpy.ndarray
            Chromosome of a single individual.
        
        Returns
        -------
        out : numpy.ndarray
            Chromosome of a single individual.
        """
        # store the origin individual
        ori = Individual()
        ori.X = x.copy()
        ori_evals = problem.evaluate(ori.X, *args, return_as_dictionary=True, **kwargs)
        ori.set_by_dict(**ori_evals)

        # copy and evaluate input individual (current leader)
        lead = Individual()
        lead.X = x.copy()
        lead_evals = problem.evaluate(lead.X, *args, return_as_dictionary = True, **kwargs)
        lead.set_by_dict(**lead_evals)

        # get search space elements not in original leader individual
        alleles = self.setspace[np.logical_not(np.in1d(self.setspace, lead.X))]

        # calculate the number of loci and number of alternative alleles
        nloci = len(lead.X)
        nalleles = len(alleles)

        # calculate the number of hillclimb/mutation attempts
        maxhc = nloci if self.maxhc is None else self.maxhc
        maxhcstep = self.maxhc if self.maxhcstep is None else self.maxhcstep

        # calculate the loci to mutate using a tiled sampling method
        lociix = tiled_choice(nloci, maxhc)

        # create lists to house non-dominated individuals
        ndom = []

        # add leader to the non-dominated list
        ndom.append(lead)

        # get minimums of each objective for non-dominated individuals
        ndom_min = lead.X.copy()

        # create counter varables
        i = 0           # counter for iterating over lociix
        nhc = 0         # counter for number of hillclimb attempts (solutions tested)
        nhcstep = 0     # counter for number of hillclimb steps made

        while nhc < maxhc and nhcstep < maxhcstep:
            # get index of allele to exchange
            j = np.random.choice(nalleles)

            # copy the leader, modify the copy, evaluate new solution
            prop = Individual()
            prop = copy.deepcopy(lead)
            prop.X[i], alleles[j] = alleles[j], prop.X[i]
            prop_evals = problem.evaluate(prop.X, *args, return_as_dictionary = True, **kwargs)
            prop.set_by_dict(**prop_evals)

            # get dominance relationships between proposed and non-dominated pool
            ndrel = np.array([get_relation(prop, indiv) for indiv in ndom])

            # determine if the proposed solution beats any of the minimums
            prop_isextreme = np.any(prop.X < ndom_min)

            # if proposed is dominated by at least one non-dominated individual, reject and revert
            if np.any(ndrel == -1):
                prop.X[i], alleles[j] = alleles[j], prop.X[i]

            elif None:
                pass

            # for each individual in the non-dominated list
            for indiv in ndom:
                # get dominance relationship between proposed and individual in non-dominated set
                rel = get_relation(prop, indiv)
            
                # if proposed is dominated by a non-dominated solution, reject and revert
                if rel == -1:
                    break
            
                # if the lead is non-dominated by the proposed solution, stash the 
                # proposed and return lead vector back to its original state.
                elif rel == 0:
                    # if proposed solution 
                    if np.any(prop.X < ndom_min):
                        ndom.append(prop)
                        lead = prop
                    prop.X[i], alleles[j] = alleles[j], prop.X[i]

                # if the lead dominates the proposed solution, 
                # return the lead vector back to its original state
                elif rel == 1:
                    prop.X[i], alleles[j] = alleles[j], prop.X[i]

            # increment counters
            nhcstep += 1
            i += 1      # increment to next lociix value
            nhc += 1    # increment to next hillclimb attempt

        # for each locus index
        for i in lociix:
            # get index of allele to exchange
            j = np.random.choice(nalleles)

            # copy the leader, modify the copy, evaluate new solution
            prop = Individual()
            prop = copy.deepcopy(lead)
            prop.X[i], alleles[j] = alleles[j], prop.X[i]
            prop_evals = problem.evaluate(prop.X, *args, return_as_dictionary = True, **kwargs)
            prop.set_by_dict(**prop_evals)
            
            # get dominance relationship between two individuals
            rel = get_relation(lead, prop)
            
            # if lead is dominated by proposed solution, keep the proposed as the new leader
            if rel == -1:
                lead_improved = True
                lead = prop
            
            # if the lead is non-dominated by the proposed solution, stash the 
            # proposed and return lead vector back to its original state.
            elif rel == 0:
                ndom.append(prop)
                prop.X[i], alleles[j] = alleles[j], prop.X[i]

            # if the lead dominates the proposed solution, 
            # return the lead vector back to its original state
            elif rel == 1:
                prop.X[i], alleles[j] = alleles[j], prop.X[i]

        # if the leader was improved, add to non-dominated individuals, otherwise keep out
        if lead_improved:
            ndom.append(lead)
        
        # if there is nothing to select from, revert to simple mutation
        if len(ndom) == 0:
            return self.reduced_exchange(problem, x, *args, **kwargs)

        # create a new population from non-dominated individual pool
        pop = Population(individuals = ndom)

        # perform non-dominated sorting on populations
        nds = NonDominatedSorting()

        # get non-dominated front indices
        ix = nds.do(pop.get("F"), only_non_dominated_front = True)

        # subset population for individuals along the non-dominated front
        pop = pop[ix]

        # get F values for population
        F = pop.get("F")

        # identify the indices of the extreme individuals
        Famax = np.argmin(F, axis = 0)

        # randomly select an index of an extreme individual to select
        Fix = np.random.choice(Famax)

        # select that individual
        indiv = pop[Fix]

        # get chromosome of individual
        out = indiv.X

        return out

    ########### Abstract method implementations ############
    def _do(
            self, 
            problem: Problem, 
            X: np.ndarray, 
            **kwargs: dict
        ) -> np.ndarray:
        """
        Perform exchange mutation for subsets + hillclimb

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

        for i in range(n_indiv):
            rnd = np.random.random()
            if rnd < self.phc:
                Xm[i,:] = self.hillclimb(problem, Xm[i,:])
            else:
                Xm[i,:] = self.reduced_exchange(problem, Xm[i,:])
        
        return Xm

