"""
Module implementing selection protocols for maximum mean expected heterozygosity selection.
"""

__all__ = [
    "OptimalContributionBaseSelection",
    "OptimalContributionBinarySelection",
    "OptimalContributionIntegerSelection",
    "OptimalContributionRealSelection",
    "OptimalContributionSubsetSelection"
]

from numbers import Integral, Real
import numpy
from numpy.random import Generator, RandomState
from typing import Optional, Union
from typing import Callable
from pybrops.breed.prot.sel.prob.OptimalContributionSelectionProblem import OptimalContributionBinarySelectionProblem, OptimalContributionIntegerSelectionProblem, OptimalContributionRealSelectionProblem, OptimalContributionSubsetSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.opt.algo.ConstrainedSteepestDescentSubsetHillClimber import ConstrainedSteepestDescentSubsetHillClimber
from pybrops.opt.algo.MemeticNSGA2SetGeneticAlgorithm import MemeticNSGA2SetGeneticAlgorithm
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm, check_is_ConstrainedOptimizationAlgorithm
from pybrops.breed.prot.sel.ConstrainedSelectionProtocol import ConstrainedSelectionProtocol
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Real, check_is_int
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.random.prng import global_prng
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix, check_is_BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix, check_is_GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory, check_is_CoancestryMatrixFactory
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class OptimalContributionBaseSelection(ConstrainedSelectionProtocol):
    """
    Semi-abstract class for Optimal Contribution Selection (OCS) with constraints.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
            cmatfcty: CoancestryMatrixFactory,
            descale: bool,
            method: str,
            nobj: Integral,
            obj_wt: Optional[numpy.ndarray] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[numpy.ndarray] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[numpy.ndarray] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[ConstrainedOptimizationAlgorithm] = None,
            moalgo: Optional[ConstrainedOptimizationAlgorithm] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for ConventionalGenomicSelection.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(OptimalContributionBaseSelection, self).__init__(
            method = method,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            ndset_wt = ndset_wt,
            ndset_trans = ndset_trans,
            ndset_trans_kwargs = ndset_trans_kwargs
            **kwargs
        )
        # order dependent assignments
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.cmatfcty = cmatfcty
        self.descale = descale
        self.rng = rng
        self.soalgo = soalgo
        self.moalgo = moalgo

    ############################ Object Properties #############################
    @property
    def nparent(self) -> int:
        """Number of parents to select."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set number of parents to select."""
        check_is_int(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value

    @property
    def ncross(self) -> int:
        """Number of crosses per configuration."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: int) -> None:
        """Set number of crosses per configuration."""
        check_is_int(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
        self._ncross = value

    @property
    def nprogeny(self) -> int:
        """Number of progeny to derive from each cross configuration."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: int) -> None:
        """Set number of progeny to derive from each cross configuration."""
        check_is_int(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value

    @property
    def cmatfcty(self) -> CoancestryMatrixFactory:
        """Factory for creating a coancestry matrix for use in the optimization."""
        return self._cmatfcty
    @cmatfcty.setter
    def cmatfcty(self, value: CoancestryMatrixFactory) -> None:
        """Set coancestry matrix factory."""
        check_is_CoancestryMatrixFactory(value, "cmatfcty")
        self._cmatfcty = value

    @property
    def descale(self) -> bool:
        """Whether to use descaled and decentered breeding values for optimization."""
        return self._descale
    @descale.setter
    def descale(self, value: bool) -> None:
        """Set descale."""
        self._descale = value

    @property
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[numpy.random.Generator,numpy.random.RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
        self._rng = value

    @property
    def soalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Description for property soalgo."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: ConstrainedOptimizationAlgorithm) -> None:
        """Set data for property soalgo."""
        if value is None:
            value = ConstrainedSteepestDescentSubsetHillClimber(rng = self.rng)
        check_is_ConstrainedOptimizationAlgorithm(value, "soalgo")
        self._soalgo = value

    @property
    def moalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Description for property moalgo."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: ConstrainedOptimizationAlgorithm) -> None:
        """Set data for property moalgo."""
        if value is None:
            value = MemeticNSGA2SetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                mu = 100,       # number of parents in population
                lamb = 100,     # number of progeny to produce
                M = 1.5,        # algorithm crossover genetic map length
                mememu = 15,    # number to local search in parent population
                memelamb = 15,  # number to local search in progeny population
                rng = self.rng  # PRNG source
            )
        check_is_ConstrainedOptimizationAlgorithm(value, "moalgo")
        self._moalgo = value

    ######################### Private Object Methods ###########################
    def _calc_bv(self, bvmat: BreedingValueMatrix) -> numpy.ndarray:
        """
        Construct a breeding value matrix for use by a Problem specification.

        Parameters
        ----------
        bvmat : BreedingValueMatrix
            Input breeding value matrix.
        
        Returns
        -------
        out : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. May be scaled or descaled.
        """
        # get breeding value matrix (n,t)
        out = bvmat.descale() if self.descale else bvmat.mat

        return out
    
    def _calc_C(self, gmat: GenotypeMatrix) -> numpy.ndarray:
        """
        Construct a Cholesky decomposition of a coancestry matrix.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix from which to calculate a coancestry matrix.
        
        Returns
        -------
        out : numpy.ndarray
            A Cholesky decomposition of a coancestry matrix of shape ``(n,n)``.
        """
        # get genomic relationship matrix: (n,n)
        G = self.cmatfcty.from_gmat(gmat)

        # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
        # if we are unable to fix, then raise value error
        if not G.apply_jitter():
            raise ValueError(
                "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                "    This could be caused by lack of genetic diversity.\n"
            )

        # convert G to (1/2)G (kinship analogue): (n,n)
        K = G.mat_asformat("kinship")

        # cholesky decomposition of K matrix: (n,n)
        out = numpy.linalg.cholesky(K).T

        return out

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    # leave problem() abstract

    ############## Pareto Frontier Functions ###############
    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
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
        # construct the optimization problem
        prob = self.problem(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max
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
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotype matrix from which to calculate genomic relationships.
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
        miscout : dict, None, default = None
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
        # check inputs
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_BreedingValueMatrix(bvmat, "bvmat")
        check_is_Real(t_cur, "t_cur")
        check_is_Real(t_max, "t_max")

        # Solve problem using a single objective method
        if self.method == "single":
            # create optimization problem
            prob = self.problem(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max
            )

            # miscellaneous output
            misc = {}

            # use single-objective optimization to get a solution.
            soln = self.soalgo.minimize(
                prob = prob,
                miscout = miscout
            )

            # get first and only solution
            obj = soln.soln_obj[0]
            sel = soln.soln_decn[0]

            # shuffle selection to ensure random mating
            numpy.random.shuffle(sel)

            # add optimization details to miscellaneous output
            if miscout is not None:
                miscout["obj"] = obj
                miscout["sel"] = sel
                miscout.update(misc) # add dict to dict

            return pgmat, sel, self.ncross, self.nprogeny

        # estimate Pareto frontier, then choose from non-dominated points.
        elif self.method == "pareto":
            # raises error
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

class OptimalContributionSubsetSelection(OptimalContributionBaseSelection):
    """
    Class defining Optimal Contribution Selection (OCS) for a subset search space.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
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
        # calculate BVs and coancestry for individuals
        bv = self._calc_bv(bvmat)
        C = self._calc_C(gmat)

        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nparent)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nparent)

        # construct problem
        prob = OptimalContributionSubsetSelectionProblem(
            bv = bv,
            C = C,
            ndecn = self.nparent,
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
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ############## Pareto Frontier Functions ###############
    # inherit pareto() implementation

    ################# Selection Functions ##################
    # inherit select() implementation

class OptimalContributionRealSelection(OptimalContributionBaseSelection):
    """
    Class defining Optimal Contribution Selection (OCS) for a subset search space.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
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
        # calculate BVs and coancestry for individuals
        bv = self._calc_bv(bvmat)
        C = self._calc_C(gmat)

        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, ntaxa)
        decn_space_upper = numpy.repeat(1.0, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = OptimalContributionRealSelectionProblem(
            bv = bv,
            C = C,
            ndecn = self.nparent,
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
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ############## Pareto Frontier Functions ###############
    # inherit pareto() implementation

    ################# Selection Functions ##################
    # inherit select() implementation

class OptimalContributionIntegerSelection(OptimalContributionBaseSelection):
    """
    Class defining Optimal Contribution Selection (OCS) for a subset search space.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
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
        # calculate BVs and coancestry for individuals
        bv = self._calc_bv(bvmat)
        C = self._calc_C(gmat)

        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(ntaxa, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = OptimalContributionIntegerSelectionProblem(
            bv = bv,
            C = C,
            ndecn = self.nparent,
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
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ############## Pareto Frontier Functions ###############
    # inherit pareto() implementation

    ################# Selection Functions ##################
    # inherit select() implementation

class OptimalContributionBinarySelection(OptimalContributionBaseSelection):
    """
    Class defining Optimal Contribution Selection (OCS) for a subset search space.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
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
        # calculate BVs and coancestry for individuals
        bv = self._calc_bv(bvmat)
        C = self._calc_C(gmat)

        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(1, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = OptimalContributionBinarySelectionProblem(
            bv = bv,
            C = C,
            ndecn = self.nparent,
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
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ############## Pareto Frontier Functions ###############
    # inherit pareto() implementation

    ################# Selection Functions ##################
    # inherit select() implementation
