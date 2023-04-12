"""
Module implementing generalized weighted genomic selection protocols.
"""

from numbers import Integral, Number, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.ConstrainedSelectionProtocol import ConstrainedSelectionProtocol
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.SubsetGeneralizedWeightedGenomicSelectionProblem import SubsetGeneralizedWeightedGenomicSelectionProblem
from pybrops.breed.prot.sel.prob.trans import trans_empty, trans_identity
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Callable, check_is_Integral, check_is_Real, check_is_dict, check_is_str
from pybrops.core.error.error_value_numpy import check_ndarray_is_1d, check_ndarray_len_eq
from pybrops.core.error.error_value_python import check_Number_in_interval, check_is_gt, check_is_gteq
from pybrops.core.random import global_prng
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.opt.algo.ConstrainedNSGA2SubsetGeneticAlgorithm import ConstrainedNSGA2SubsetGeneticAlgorithm
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm, check_is_ConstrainedOptimizationAlgorithm
from pybrops.opt.algo.ConstrainedSteepestDescentSubsetHillClimber import ConstrainedSteepestDescentSubsetHillClimber
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class ConstrainedGeneralizedWeightedGenomicSelection(ConstrainedSelectionProtocol):
    """
    docstring for ConstrainedGeneralizedWeightedGenomicSelection.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
            alpha: Real,
            method: str,
            nobj: Integral,
            obj_wt: numpy.ndarray,
            obj_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            eqcv_trans_kwargs: Optional[dict],
            ndset_wt: Real = 1.0,
            ndset_trans: Optional[Callable] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Union[Generator,RandomState] = global_prng, 
            soalgo: Optional[ConstrainedOptimizationAlgorithm] = None,
            moalgo: Optional[ConstrainedOptimizationAlgorithm] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for ConstrainedGeneralizedWeightedGenomicSelection.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(ConstrainedGeneralizedWeightedGenomicSelection, self).__init__(
            **kwargs
        )

        # make value assignments (order dependent)
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.alpha = alpha
        self.method = method
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.obj_trans = obj_trans
        self.obj_trans_kwargs = obj_trans_kwargs
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.ineqcv_trans = ineqcv_trans
        self.ineqcv_trans_kwargs = ineqcv_trans_kwargs
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt
        self.eqcv_trans = eqcv_trans
        self.eqcv_trans_kwargs = eqcv_trans_kwargs
        self.ndset_wt = ndset_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.rng = rng
        self.soalgo = soalgo
        self.moalgo = moalgo

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nparent(self) -> int:
        """Number of parents to select."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set number of parents to select."""
        check_is_Integral(value, "nparent") # must be integer
        check_is_gt(value, "nparent", 0)    # integer must be >0
        self._nparent = value
    @nparent.deleter
    def nparent(self) -> None:
        """Delete number of parents to select."""
        del self._nparent

    @property
    def ncross(self) -> int:
        """Number of crosses per configuration."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: int) -> None:
        """Set number of crosses per configuration."""
        check_is_Integral(value, "ncross")  # must be integer
        check_is_gt(value, "ncross", 0)     # integer must be >0
        self._ncross = value
    @ncross.deleter
    def ncross(self) -> None:
        """Delete number of crosses per configuration."""
        del self._ncross

    @property
    def nprogeny(self) -> int:
        """Number of progeny to derive from each cross configuration."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: int) -> None:
        """Set number of progeny to derive from each cross configuration."""
        check_is_Integral(value, "nprogeny")    # must be integer
        check_is_gt(value, "nprogeny", 0)       # integer must be >0
        self._nprogeny = value
    @nprogeny.deleter
    def nprogeny(self) -> None:
        """Delete number of progeny to derive from each cross configuration."""
        del self._nprogeny

    @property
    def alpha(self) -> Real:
        """Exponent to which to raise the favorable allele frequency. Must be in the range [0,1]."""
        return self._alpha
    @alpha.setter
    def alpha(self, value: Real) -> None:
        """Set exponent to which to raise the favorable allele frequency."""
        check_is_Real(value, "alpha")
        check_Number_in_interval(value, "alpha", 0, 1)
        self._alpha = value
    @alpha.deleter
    def alpha(self) -> None:
        """Delete exponent to which to raise the favorable allele frequency."""
        del self._alpha

    @property
    def method(self) -> str:
        """Selection method."""
        return self._method
    @method.setter
    def method(self, value: str) -> None:
        """Set selection method."""
        check_is_str(value, "method")       # must be string
        value = value.lower()               # convert to lowercase
        options = ("single", "pareto")      # method options
        # if not method supported raise ValueError
        if value not in options:
            raise ValueError("Unsupported 'method'. Options are: " + ", ".join(map(str, options)))
        self._method = value
    @method.deleter
    def method(self) -> None:
        """Delete selection method."""
        del self._method

    @property
    def nobj(self) -> Integral:
        """Number of objectives."""
        return self._nobj
    @nobj.setter
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        check_is_Integral(value, "nobj")
        check_is_gt(value, "nobj", 0)     # cannot have 0 objectives
        self._nobj = value
        self._n_obj = value # for easy separation from PyMOO
    @nobj.deleter
    def nobj(self) -> None:
        """Delete number of objectives."""
        del self._nobj
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        return self._obj_wt
    @obj_wt.setter
    def obj_wt(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set objective function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "obj_wt")
            check_ndarray_len_eq(value, "obj_wt", self.nobj)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.nobj)
        else:
            raise TypeError("'obj_wt' must be of type numpy.ndarray or a numeric type")
        self._obj_wt = value
    @obj_wt.deleter
    def obj_wt(self) -> None:
        """Delete objective function weights."""
        del self._obj_wt

    @property
    def obj_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to objective function values."""
        return self._obj_trans
    @obj_trans.setter
    def obj_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to objective space transformation function."""
        if value is None:
            value = trans_identity
        check_is_Callable(value, "obj_trans")
        self._obj_trans = value
    @obj_trans.deleter
    def obj_trans(self) -> None:
        """Delete latent space to objective space transformation function."""
        del self._obj_trans
    
    @property
    def obj_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to objective space transformation function."""
        return self._obj_trans_kwargs
    @obj_trans_kwargs.setter
    def obj_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to objective space transformation function."""
        if value is None:
            value = {}
        check_is_dict(value, "obj_trans_kwargs")
        self._obj_trans_kwargs = value
    @obj_trans_kwargs.deleter
    def obj_trans_kwargs(self) -> None:
        """Delete keyword arguments for the latent space to objective space transformation function."""
        del self._obj_trans_kwargs
    
    @property
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violation functions."""
        return self._nineqcv
    @nineqcv.setter
    def nineqcv(self, value: Integral) -> None:
        """Set number of inequality constraint violation functions."""
        check_is_Integral(value, "nineqcv")
        check_is_gteq(value, "nineqcv", 0)  # possible to have 0 inequality constraints
        self._nineqcv = value
        self._n_ieq_constr = value # for easy separation from PyMOO
    @nineqcv.deleter
    def nineqcv(self) -> None:
        """Delete number of inequality constraint violation functions."""
        del self._nineqcv

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        return self._ineqcv_wt
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: Union[numpy.ndarray,Number]) -> None:
        """Set inequality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "ineqcv_wt")
            check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.nineqcv)
        else:
            raise TypeError("'ineqcv_wt' must be of type numpy.ndarray or a numeric type")
        self._ineqcv_wt = value
    @ineqcv_wt.deleter
    def ineqcv_wt(self) -> None:
        """Delete inequality constraint violation function weights."""
        del self._ineqcv_wt

    @property
    def ineqcv_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to inequality constraint violation values."""
        return self._ineqcv_trans
    @ineqcv_trans.setter
    def ineqcv_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to inequality constraint violation transformation function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "ineqcv_trans")
        self._ineqcv_trans = value
    @ineqcv_trans.deleter
    def ineqcv_trans(self) -> None:
        """Delete latent space to inequality constraint violation transformation function."""
        del self._ineqcv_trans
    
    @property
    def ineqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to inequality constraint violation transformation function."""
        return self._ineqcv_trans_kwargs
    @ineqcv_trans_kwargs.setter
    def ineqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to inequality constraint violation transformation function."""
        if value is None:
            value = {}
        check_is_dict(value, "ineqcv_trans_kwargs")
        self._ineqcv_trans_kwargs = value
    @ineqcv_trans_kwargs.deleter
    def ineqcv_trans_kwargs(self) -> None:
        """Delete keyword arguments for the latent space to inequality constraint violation transformation function."""
        del self._ineqcv_trans_kwargs
    
    @property
    def neqcv(self) -> Integral:
        """Number of equality constraint violations."""
        return self._neqcv
    @neqcv.setter
    def neqcv(self, value: Integral) -> None:
        """Set number of equality constraint violations."""
        check_is_Integral(value, "neqcv")
        check_is_gteq(value, "neqcv", 0)    # possible to have 0 equality constraints
        self._neqcv = value
        self._n_eq_constr = value # for easy separation from PyMOO
    @neqcv.deleter
    def neqcv(self) -> None:
        """Delete number of equality constraint violations."""
        del self._neqcv
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        return self._eqcv_wt
    @eqcv_wt.setter
    def eqcv_wt(self, value: Union[numpy.ndarray,Number]) -> None:
        """Set equality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "eqcv_wt")
            check_ndarray_len_eq(value, "eqcv_wt", self.neqcv)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.neqcv)
        else:
            raise TypeError("'eqcv_wt' must be of type numpy.ndarray or a numeric type")
        self._eqcv_wt = value
    @eqcv_wt.deleter
    def eqcv_wt(self) -> None:
        """Delete equality constraint violation function weights."""
        del self._eqcv_wt

    @property
    def eqcv_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to equality constraint violation values."""
        return self._eqcv_trans
    @eqcv_trans.setter
    def eqcv_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to equality constraint violation transformation function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "eqcv_trans")
        self._eqcv_trans = value 
    @eqcv_trans.deleter
    def eqcv_trans(self) -> None:
        """Delete latent space to equality constraint violation transformation function."""
        del self._eqcv_trans
    
    @property
    def eqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to equality constraint violation transformation function."""
        return self._eqcv_trans_kwargs
    @eqcv_trans_kwargs.setter
    def eqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to equality constraint violation transformation function."""
        if value is None:
            value = {}
        check_is_dict(value, "eqcv_trans_kwargs")
        self._eqcv_trans_kwargs = value
    @eqcv_trans_kwargs.deleter
    def eqcv_trans_kwargs(self) -> None:
        """Delete keyword arguments for the latent space to equality constraint violation transformation function."""
        del self._eqcv_trans_kwargs

    @property
    def ndset_wt(self) -> Union[float,numpy.ndarray]:
        """Nondominated set weights."""
        return self._ndset_wt
    @ndset_wt.setter
    def ndset_wt(self, value: Union[float,numpy.ndarray]) -> None:
        """Set nondominated set weights."""
        self._ndset_wt = value
    @ndset_wt.deleter
    def ndset_wt(self) -> None:
        """Delete nondominated set weights."""
        del self._ndset_wt

    @property
    def ndset_trans(self) -> Union[Callable,None]:
        """Nondominated set transformation function."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable,None]) -> None:
        """Set nondominated set transformation function."""
        if value is not None:                       # if given object
            check_is_Callable(value, "ndset_trans") # must be callable
        self._ndset_trans = value
    @ndset_trans.deleter
    def ndset_trans(self) -> None:
        """Delete nondominated set transformation function."""
        del self._ndset_trans

    @property
    def ndset_trans_kwargs(self) -> dict:
        """Nondominated set transformation function keyword arguments."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.setter
    def ndset_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set nondominated set transformation function keyword arguments."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "ndset_trans_kwargs")  # check is dict
        self._ndset_trans_kwargs = value
    @ndset_trans_kwargs.deleter
    def ndset_trans_kwargs(self) -> None:
        """Delete nondominated set transformation function keyword arguments."""
        del self._ndset_trans_kwargs

    @property
    def rng(self) -> Union[Generator,RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
        self._rng = value
    @rng.deleter
    def rng(self) -> None:
        """Delete random number generator source."""
        del self._rng

    @property
    def soalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Single-objective optimization algorithm."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: ConstrainedOptimizationAlgorithm) -> None:
        """Set single-objective optimization algorithm."""
        if value is None:
            # construct default hillclimber
            value = ConstrainedSteepestDescentSubsetHillClimber(
                self.rng
            )
        check_is_ConstrainedOptimizationAlgorithm(value, "soalgo")
        self._soalgo = value
    @soalgo.deleter
    def soalgo(self) -> None:
        """Delete single-objective optimization algorithm."""
        del self._soalgo

    @property
    def moalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[ConstrainedOptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = ConstrainedNSGA2SubsetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100, # number of parents in population
                rng = self.rng  # PRNG source
            )
        check_is_ConstrainedOptimizationAlgorithm(value, "moalgo")
        self._moalgo = value
    @moalgo.deleter
    def moalgo(self) -> None:
        """Delete multi-objective opimization algorithm."""
        del self._moalgo


    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: AdditiveLinearGenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> SelectionProblem:
        """
        Create an optimization problem definition using provided inputs.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SelectionProblem
            An optimization problem definition.
        """
        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(ntaxa-1, ntaxa)

        # construct problem
        prob = SubsetGeneralizedWeightedGenomicSelectionProblem(
            Z_a = gmat.mat,
            u_a = gpmod.u_a,
            fafreq = gpmod.fafreq(gmat),
            alpha = self.alpha,
            ndecn = ntaxa,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = self.nobj,
            obj_wt = self.obj_wt,
            obj_trans = self.obj_trans,
            obj_trans_kwargs = self.obj_trans_kwargs,
            nineqcv = self.nineqcv,
            ineqcv_wt = self.ineqcv_wt,
            ineqcv_trans = self.ineqcv_trans,
            ineqcv_trans_kwargs = self.ineqcv_trans_kwargs,
            neqcv = self.neqcv,
            eqcv_wt = self.eqcv_wt,
            eqcv_trans = self.eqcv_trans,
            eqcv_trans_kwargs = self.eqcv_trans_kwargs,
            **kwargs
        )

        return prob

    ############## Pareto Frontier Functions ###############
    def pareto(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: AdditiveLinearGenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> tuple:
        """
        Calculate a Pareto frontier for objectives.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing two objects ``(frontier, sel_config)``.

            Where:

            - ``frontier`` is a ``numpy.ndarray`` of shape ``(q,v)`` containing
              Pareto frontier points.
            - ``sel_config`` is a ``numpy.ndarray`` of shape ``(q,k)`` containing
              parent selection decisions for each corresponding point in the
              Pareto frontier.

            Where:

            - ``q`` is the number of points in the frontier.
            - ``v`` is the number of objectives for the frontier.
            - ``k`` is the number of search space decision variables.
        """
        # construct the problem
        prob = self.problem(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max,
            **kwargs
        )

        # optimize the problem
        soln = self.moalgo.minimize(
            prob = prob,
            miscout = miscout
        )

        if miscout is not None:
            miscout["soln"] = soln

        frontier = soln.soln_obj
        sel_config = soln.soln_decn

        return frontier, sel_config

    ################# Selection Functions ##################
    def select(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: AdditiveLinearGenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> tuple:
        """
        Select individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing four objects: ``(pgmat, sel, ncross, nprogeny)``.

            Where:

            - ``pgmat`` is a PhasedGenotypeMatrix of parental candidates.
            - ``sel`` is a ``numpy.ndarray`` of indices specifying a cross
              pattern. Each index corresponds to an individual in ``pgmat``.
            - ``ncross`` is a ``numpy.ndarray`` specifying the number of
              crosses to perform per cross pattern.
            - ``nprogeny`` is a ``numpy.ndarray`` specifying the number of
              progeny to generate per cross.
        """
        # single-objective method: objfn_trans returns a single value for each
        # selection configuration
        if self.method == "single":
            # construct the problem
            prob = self.problem(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                **kwargs
            )

            # optimize the problem
            soln = self.soalgo.minimize(
                prob = prob,
                miscout = miscout
            )

            if miscout is not None:
                miscout["soln"] = soln

            # extract decision variables
            sel = soln.soln_decn[0]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

            return pgmat, sel, self.ncross, self.nprogeny

        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif self.method == "pareto":
            # get the pareto frontier
            frontier, sel_config = self.pareto(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                miscout = miscout,
                **kwargs
            )

            # get scores for each of the points along the pareto frontier
            score = self.ndset_wt * self.ndset_trans(frontier, **self.ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to miscout
            if miscout is not None:
                miscout["frontier"] = frontier
                miscout["sel_config"] = sel_config

            return pgmat, sel_config[ix], self.ncross, self.nprogeny
        else:
            raise ValueError("argument 'method' must be either 'single' or 'pareto'")
