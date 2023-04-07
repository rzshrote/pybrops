"""
Module implementing selection protocols for maximum mean expected heterozygosity selection.
"""

from numbers import Integral, Real
import warnings
import numpy
from typing import Optional, Tuple, Union
from typing import Callable
from pybrops.breed.prot.sel.prob.BinaryOptimalContributionSelectionProblem import BinaryOptimalContributionSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.core.error.error_value_numpy import check_ndarray_align
from pybrops.opt.algo.ConstrainedSteepestDescentSubsetHillClimber import ConstrainedSteepestDescentSubsetHillClimber
from pybrops.opt.algo.MemeticNSGA2SetGeneticAlgorithm import MemeticNSGA2SetGeneticAlgorithm
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm, check_is_ConstrainedOptimizationAlgorithm
from pybrops.breed.prot.sel.ConstrainedSelectionProtocol import ConstrainedSelectionProtocol
from pybrops.core.error.error_attr_python import check_is_callable
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Real, check_is_dict, check_is_int, check_is_str
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.random.prng import global_prng
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix, check_is_BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix, check_is_GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix
from pybrops.breed.prot.sel.transfn import trans_identity_unconstrained
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory, check_is_CoancestryMatrixFactory
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class BinaryOptimalContributionSelection(ConstrainedSelectionProtocol):
    """
    Class implementing selection protocols for optimal mean expected heterozygosity selection.

    Maximum Mean Expected Heterozygosity Selection (MMEHS) is defined as:

    .. math::
        \\max_{\\mathbf{x}} f_{MMEHS}(\\mathbf{x}) = 1 - \\mathbf{x'Kx}

    With constraints:

    .. math::
        \\mathbf{1_{n}'x} = 1

        \\mathbf{x} \\in \\mathbb{R}^n

    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, 
            nparent: int, 
            ncross: int, 
            nprogeny: int,
            cmatfcty: CoancestryMatrixFactory,
            method: str = "single",
            descale: bool = False,
            encode_trans: Optional[Callable[[numpy.ndarray,dict],Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]]] = None,
            encode_trans_kwargs: Optional[dict] = None,
            nobj: Optional[Integral] = None,
            obj_wt: Optional[numpy.ndarray] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[numpy.ndarray] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[numpy.ndarray] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            ndset_wt: Optional[Real] = None,
            rng = global_prng, 
            soalgo: Optional[ConstrainedOptimizationAlgorithm] = None,
            moalgo: Optional[ConstrainedOptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for Optimal Contribution Selection (OCS).

        Parameters
        ----------
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        cmatfcty : CoancestryMatrixFactory
            Factory object from which to generate a coancestry matrix.
        method : str
            Optimization strategy.

            +--------+---------------------------------------------------------+
            | Option | Description                                             |
            +========+=========================================================+
            | single | Transform all breeding values into a single overall     |
            |        | breeding value using the function ``objfn_trans``. Then |
            |        | solve for OCS with a diversity constraint using         |
            |        | transformed breeding values.                            |
            +--------+---------------------------------------------------------+
            | pareto | Treat inbreeding and each trait as different            |
            |        | objectives. Transform this list of objectives using     |
            |        | ``objfn_trans`` to get a list of transformed            |
            |        | objectives. Approximate the Pareto by identifying a set |
            |        | of non-dominated points along each transformed          |
            |        | objective. Then apply ``ndset_trans`` to score the      |
            |        | non-dominated points.                                   |
            +--------+---------------------------------------------------------+
        objfn_trans : function or callable
            Function to transform the OCS objective function.

            If method = "single", this function must accept an array of length
            ``(t,)`` and return a scalar, where ``t`` is the number of trait
            breeding values for an individual.

            If method = "pareto", this function must accept an array of length
            ``(1+t,)`` and return a numpy.ndarray, where ``t`` is the number of
            trait breeding values for an individual. The first element of the
            input array is the mean inbreeding coefficient for the selection.

            General function definition::

                objfn_trans(obj, **kwargs: dict):
                    Parameters
                        obj : scalar, numpy.ndarray
                            Objective scalar or vector to be transformed
                        kwargs : dict
                            Additional keyword arguments
                    Returns
                        out : scalar, numpy.ndarray
                            Transformed objective scalar or vector.
        objfn_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'objfn_trans'.
        objfn_wt : float, numpy.ndarray
            Weight applied to transformed objective function. Indicates whether
            a function is maximizing or minimizing:

            - ``1.0`` for maximizing function.
            - ``-1.0`` for minimizing function.
        ndset_trans : numpy.ndarray
            Function to transform nondominated points along the Pareto frontier
            into a single score for each point.

            Function definition::

                ndset_trans(ndset, **kwargs: dict):
                    Parameters
                        ndset : numpy.ndarray
                            Array of shape (j,o) containing nondominated points.
                            Where 'j' is the number of nondominated points and
                            'o' is the number of objectives.
                        kwargs : dict
                            Additional keyword arguments.
                    Returns
                        out : numpy.ndarray
                            Array of shape (j,) containing transformed Pareto
                            frontier points.
        ndset_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'ndset_trans'.
        ndset_wt : float
            Weight applied to transformed nondominated points along Pareto
            frontier. Indicates whether a function is maximizing or minimizing.
                1.0 for maximizing function.
                -1.0 for minimizing function.
        soalgo : ConstrainedOptimizationAlgorithm
            Single objective optimization algorithm with which to optimize the
            optimization problem.
        moalgo : ConstrainedOptimizationAlgorithm
            Multi-objective optimization algorithm to optimize the objective
            functions. If ``None``, use a NSGA3UnityConstraintGeneticAlgorithm
            with the following parameters::

                moalgo = NSGA3UnityConstraintGeneticAlgorithm(
                    ngen = 250,             # number of generations to evolve
                    mu = 100,               # number of parents in population
                    lamb = 100,             # number of progeny to produce
                    cxeta = 30.0,           # crossover variance parameter
                    muteta = 20.0,          # mutation crossover parameter
                    refpnts = None,         # hyperplane reference points
                    save_logbook = False,   # whether to save logs or not
                    rng = self.rng          # PRNG source
                )
        rng : numpy.random.Generator or None
            A random number generator source. Used for optimization algorithms.
        """
        super(BinaryOptimalContributionSelection, self).__init__(**kwargs)
        
        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.cmatfcty = cmatfcty
        self.method = method
        self.descale = descale

        self.encode_trans = encode_trans
        self.encode_trans_kwargs = encode_trans_kwargs
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt

        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs # property replaces None with {}
        self.ndset_wt = ndset_wt
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
        check_is_int(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
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
        check_is_int(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
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
        check_is_int(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value
    @nprogeny.deleter
    def nprogeny(self) -> None:
        """Delete number of progeny to derive from each cross configuration."""
        del self._nprogeny

    @property
    def cmatfcty(self) -> CoancestryMatrixFactory:
        """Factory for creating a coancestry matrix for use in the optimization."""
        return self._cmatfcty
    @cmatfcty.setter
    def cmatfcty(self, value: CoancestryMatrixFactory) -> None:
        """Set coancestry matrix factory."""
        check_is_CoancestryMatrixFactory(value, "cmatfcty")
        self._cmatfcty = value
    @cmatfcty.deleter
    def cmatfcty(self) -> None:
        """Delete coancestry matrix factory."""
        del self._cmatfcty

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
    def descale(self) -> bool:
        """Whether to use descaled and decentered breeding values for optimization."""
        return self._descale
    @descale.setter
    def descale(self, value: bool) -> None:
        """Set descale."""
        self._descale = value
    @descale.deleter
    def descale(self) -> None:
        """Delete descale."""
        del self._descale

    @property
    def ndset_trans(self) -> Union[Callable,None]:
        """Nondominated set transformation function."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable,None]) -> None:
        """Set nondominated set transformation function."""
        if value is not None:                       # if given object
            check_is_callable(value, "ndset_trans") # must be callable
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
    @rng.deleter
    def rng(self) -> None:
        """Delete random number generator source."""
        del self._rng

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
    @soalgo.deleter
    def soalgo(self) -> None:
        """Delete data for property soalgo."""
        del self._soalgo

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
    @moalgo.deleter
    def moalgo(self) -> None:
        """Delete data for property moalgo."""
        del self._moalgo

    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
    def _calc_inbmax(self, K, t_cur, t_max):
        inb_constraint = self.inbfn(t_cur, t_max)           # get target inbreeding level
        inb_lower_bound = 1.0 / numpy.linalg.inv(K).sum()   # 1 / (1'K^(-1)1)
        inb_upper_bound = numpy.max(numpy.trace(K))         # max(trace(K))

        # test for infeasibility
        if inb_constraint < inb_lower_bound:
            warnings.warn(
                "Provided inbreeding target of {0} is infeasible.\n".format(inb_constraint) +
                "Increasing inbreeding target to {0} which is in feasible region...".format(inb_lower_bound + 1e-6)
            )
            # give some wiggle room for numerical inaccuracies due to matrix
            # inversion, which is by nature unstable.
            inb_constraint = inb_lower_bound + 1e-6
        if inb_constraint > inb_upper_bound:
            warnings.warn(
                "Provided inbreeding target of {0} is infeasible.\n".format(inb_constraint) +
                "Decreasing inbreeding target to {0} which is in feasible region...".format(inb_upper_bound)
            )
            inb_constraint = inb_upper_bound
        
        return inb_constraint

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
        # check inputs
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_BreedingValueMatrix(bvmat, "bvmat")

        # check dimensions
        if gmat.ntaxa != bvmat.ntaxa:
            raise ValueError("'gmat' and 'bvmat' do not have same number of taxa")

        # make sure gmat and bvmat are aligned
        if (gmat.taxa is not None) and (bvmat.taxa is not None):
            check_ndarray_align(gmat.taxa, "gmat.taxa", bvmat.taxa, "bvmat.taxa")

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
        C = numpy.linalg.cholesky(K).T
        
        # get breeding value matrix (n,t)
        bv = bvmat.descale() if self.descale else bvmat.mat

        # the number of parents is the number of decision variables
        ndecn = self.nparent

        # the decision space is determined by the number of taxa
        decn_space = numpy.arange(bvmat.ntaxa)

        # if encode_trans is None, provide a default
        if self.encode_trans is None:
            encode_trans = trans_identity_unconstrained
            encode_trans_kwargs = {}
            # set objectives to minimize for inbreeding, maximize for breeding values
            nobj = 1 + bvmat.ntrait
            obj_wt = numpy.array([-1.0] + bvmat.ntrait * [1.0], dtype = float)
            # no inequality constraints
            nineqcv = 0
            ineqcv_wt = numpy.array([], dtype = float)
            # no equality constraints
            neqcv = 0
            eqcv_wt = numpy.array([], dtype = float)
        else:
            encode_trans = self.encode_trans
            encode_trans_kwargs = self.encode_trans_kwargs
            nobj = self.nobj
            obj_wt = self.obj_wt
            nineqcv = self.nineqcv
            ineqcv_wt = self.ineqcv_wt
            neqcv = self.neqcv
            eqcv_wt = self.eqcv_wt
            
        # construct problem
        prob = BinaryOptimalContributionSelectionProblem(
            bv = bv,
            C = C,
            ndecn = ndecn,
            decn_space = decn_space,
            encode_trans = encode_trans,
            encode_trans_kwargs = encode_trans_kwargs,
            nobj = nobj,
            obj_wt = obj_wt,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt
        )

        return prob

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

        # use multi-objective optimization to approximate Pareto front.
        soln = self.moalgo.minimize(prop = prob, miscout = misc)

        # handle miscellaneous output
        if miscout is not None:     # if miscout is provided
            miscout.update(misc)    # add 'misc' to 'miscout', overwriting as needed

        return soln.soln_obj, soln.soln_decn

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
            soln = self.soalgo.minimize(prop = prob, miscout = misc)

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

