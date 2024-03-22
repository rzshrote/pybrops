"""
Module implementing selection protocols for multi-objective genomic mating.
"""

from numbers import Real
from typing import Callable
from typing import Union
import numpy
import types

from pybrops.opt.algo.UnconstrainedOptimizationAlgorithm import UnconstrainedOptimizationAlgorithm
from pybrops.opt.algo.UnconstrainedOptimizationAlgorithm import check_is_OptimizationAlgorithm
from pybrops.breed.prot.sel.targetfn import target_negative
from pybrops.breed.prot.sel.targetfn import target_positive
from pybrops.breed.prot.sel.targetfn import target_stabilizing
from pybrops.breed.prot.sel.weightfn import weight_absolute
from pybrops.breed.prot.sel.weightfn import weight_one
from pybrops.core.error.error_type_python import check_is_int_or_inf
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.core.random.prng import global_prng
from pybrops.opt.algo.UnconstrainedNSGA2SetGeneticAlgorithm import UnconstrainedNSGA2SetGeneticAlgorithm
from pybrops.opt.algo.UnconstrainedSteepestAscentSetHillClimber import UnconstrainedSteepestAscentSetHillClimber
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import UnconstrainedSelectionProtocol
from pybrops.core.error.error_type_python import check_isinstance
from pybrops.core.error.error_type_python import check_is_bool
from pybrops.core.error.error_attr_python import check_is_callable
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.util.arrayix import triudix
from pybrops.core.util.arrayix import triuix
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import check_is_GeneticVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class MultiObjectiveGenomicMating(UnconstrainedSelectionProtocol):
    """
    Class implementing selection protocols for multi-objective genomic mating.

    # TODO: add formulae for methodology.
    """

    ########################## Special Object Methods ##########################
    def __init__(self,
            nconfig: int, 
            nparent: int, 
            ncross: int, 
            nprogeny: int, 
            vmatfcty: GeneticVarianceMatrixFactory, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            weight: Union[numpy.ndarray,Callable,str] = weight_absolute,
            target: Union[numpy.ndarray,Callable,str] = target_positive,
            unique_parents: bool = True, 
            mem: int = 1024,
            method: str = "single",
            objfn_trans = None, 
            objfn_trans_kwargs = None, 
            objfn_wt = 1.0,
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = 1.0,
            rng = None, 
            soalgo = None, 
            moalgo = None,
            **kwargs: dict
        ):
        """
        Constructor for MultiObjectiveGenomicSelection class.

        Parameters
        ----------
        nconfig : int
            Number of cross configurations to consider.

            Examples:

            - 20 two-way crosses would be: ``nconfig = 20``
            - 20 three way crosses would be: ``nconfig = 20``
        nparent : int
            Number of parents to per configuration.

            Example:

            - 20 two-way crosses would be: ``nparent = 2``
            - 20 three-way crosses would be: ``nparent = 3``
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        vmatcls : class type
            Variance matrix class name from which to construct additive
            variance matrices from
        s : int
            Used for 'vmatcls' matrix construction.
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-----------------+-------------------------+
            | Example         | Description             |
            +=================+=========================+
            | ``nself = 0``   | Derive gametes from F1  |
            +-----------------+-------------------------+
            | ``nself = 1``   | Derive gametes from F2  |
            +-----------------+-------------------------+
            | ``nself = 2``   | Derive gametes from F3  |
            +-----------------+-------------------------+
            | ``...``         | etc.                    |
            +-----------------+-------------------------+
            | ``nself = inf`` | Derive gametes from SSD |
            +-----------------+-------------------------+
        gmapfn : GeneticMapFunction
            Used for 'vmatcls' matrix construction.
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : int, default = 1024
            Used for 'vmatcls' matrix construction.
            Memory chunk size to use during matrix operations. If ``None``,
            then memory chunk size is not limited.

            WARNING: Setting ``mem = None`` might result in memory allocation
            errors! For reference, ``mem = 1024`` refers to a matrix of size
            1024x1024, which needs about 8.5 MB of storage. Matrices of course
            need a quadratic amount of memory: :math:`O(n^2)`.
        unique_parents : bool, default = True
            Whether to allow force unique parents or not.
            If ``True``, all parents in the mating configuration must be unique.
            If ``False``, non-unique parents are allowed. In this scenario,
            self-fertilization is considered as a viable option.
        method : str
            Method of selecting parents.

            +--------------+---------------------------------------------------+
            | Method       | Description                                       |
            +==============+===================================================+
            | ``"single"`` | MOGM is transformed to a single objective and     |
            |              | optimization is done on the transformed function. |
            |              | This is done using the ``trans`` function         |
            |              | provided::                                        |
            |              |                                                   |
            |              |    optimize : objfn_trans(MOGM)                   |
            +--------------+---------------------------------------------------+
            | ``"pareto"`` | MOGM is transformed by a transformation function, |
            |              | but NOT reduced to a single objective. The Pareto |
            |              | frontier for this transformed function is mapped  |
            |              | using a multi-objective genetic algorithm.        |
            |              |                                                   |
            |              | Objectives are scaled to :math:`[0,1]` and a      |
            |              | vector orthogonal to the hyperplane defined by    |
            |              | the extremes of the front is drawn starting at    |
            |              | the point defined by ``ndset_trans``. The closest |
            |              | point on the Pareto frontier to the orthogonal    |
            |              | vector is selected.                               |
            +--------------+---------------------------------------------------+
        target : str or numpy.ndarray
            If target is a string, check value and follow these rules:

            +-------------------+----------------------------------------------+
            | Value             | Description                                  |
            +===================+==============================================+
            | ``"positive"``    | Select alleles with the most positive        |
            |                   | effect.                                      |
            +-------------------+----------------------------------------------+
            | ``"negative"``    | Select alleles with the most negate effect.  |
            +-------------------+----------------------------------------------+
            | ``"stabilizing"`` | Set target allele frequency to ``0.5``.      |
            +-------------------+----------------------------------------------+
            | ``numpy.ndarray`` | Use frequency values in ``target`` as is.    |
            +-------------------+----------------------------------------------+
        weight : str or numpy.ndarray
            If weight is a string, check value and follow these rules:

            +-----------------+------------------------------------------------+
            | Value           | Description                                    |
            +=================+================================================+
            | ``"magnitude"`` | Assign weights using the magnitudes of         |
            |                 | regression coefficients.                       |
            +-----------------+------------------------------------------------+
            | ``"equal"``     | Assign weights equally.                        |
            +-----------------+------------------------------------------------+
        objfn_trans : function, callable
            Function to transform the MOGM function. If method = "single", this
            function must return a scalar. If method = "pareto", this function
            must return a ``numpy.ndarray``.

            Function definition::

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
        soalgo : OptimizationAlgorithm
            Single-objective optimization algorithm to optimize the objective
            function. If ``None``, use a SteepestAscentSetHillClimber with the
            following parameters::

                soalgo = SteepestAscentSetHillClimber(
                    rng = self.rng  # PRNG source
                )
        moalgo : OptimizationAlgorithm
            Multi-objective optimization algorithm to optimize the objective
            functions. If ``None``, use a NSGA2SetGeneticAlgorithm with the
            following parameters::

                moalgo = NSGA2SetGeneticAlgorithm(
                    ngen = 250,     # number of generations to evolve
                    mu = 100,       # number of parents in population
                    lamb = 100,     # number of progeny to produce
                    M = 1.5,        # algorithm crossover genetic map length
                    rng = self.rng  # PRNG source
                )
        rng : numpy.random.Generator or None
            A random number generator source. Used for optimization algorithms.
            If ``rng`` is ``None``, use ``pybrops.core.random`` module
            (NOT THREAD SAFE!).
        """
        super(MultiObjectiveGenomicMating, self).__init__(**kwargs)

        # error checks and assignments (ORDER DEPENDENT!!!)
        self.nconfig = nconfig
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.vmatcls = vmatfcty
        self.nself = nself
        self.gmapfn = gmapfn
        self.mem = mem
        self.unique_parents = unique_parents
        self.method = method
        self.target = target
        self.weight = weight
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs # property replaces None with {}
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs # property replaces None with {}
        self.ndset_wt = ndset_wt
        self.rng = rng  # property replaces None with pybrops.core.random
        # soalgo, moalgo MUST GO AFTER 'rng'; properties provide default if None
        self.soalgo = soalgo
        self.moalgo = moalgo

    ############################ Object Properties #############################
    @property
    def nconfig(self) -> int:
        """Description for property nconfig."""
        return self._nconfig
    @nconfig.getter
    def nconfig(self) -> int:
        """Get data for property nconfig."""
        return self._nconfig
    @nconfig.setter
    def nconfig(self, value: int) -> None:
        """Set data for property nconfig."""
        check_is_int(value, "nconfig")      # must be int
        check_is_gt(value, "nconfig", 0)    # int must be >0
        self._nconfig = value
    
    @property
    def nparent(self) -> int:
        """Description for property nparent."""
        return self._nparent
    @nparent.getter
    def nparent(self) -> int:
        """Get data for property nparent."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set data for property nparent."""
        check_is_int(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value
    
    @property
    def ncross(self) -> int:
        """Description for property ncross."""
        return self._ncross
    @ncross.getter
    def ncross(self) -> int:
        """Get data for property ncross."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: int) -> None:
        """Set data for property ncross."""
        check_is_int(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
        self._ncross = value
    
    @property
    def nprogeny(self) -> int:
        """Description for property nprogeny."""
        return self._nprogeny
    @nprogeny.getter
    def nprogeny(self) -> int:
        """Get data for property nprogeny."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: int) -> None:
        """Set data for property nprogeny."""
        check_is_int(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value
    
    @property
    def vmatfcty(self) -> GeneticVarianceMatrixFactory:
        """Description for property vmatfcty."""
        return self._vmatfcty
    @vmatfcty.getter
    def vmatfcty(self) -> GeneticVarianceMatrixFactory:
        """Get data for property vmatfcty."""
        return self._vmatfcty
    @vmatfcty.setter
    def vmatfcty(self, value: GeneticVarianceMatrixFactory) -> None:
        """Set data for property vmatfcty."""
        check_is_GeneticVarianceMatrixFactory(value, "vmatfcty")
        self._vmatfcty = value

    @property
    def nself(self) -> Union[int,Real]:
        """Description for property nself."""
        return self._nself
    @nself.getter
    def nself(self) -> Union[int,Real]:
        """Get data for property nself."""
        return self._nself
    @nself.setter
    def nself(self, value: Union[int,Real]) -> None:
        """Set data for property nself."""
        check_is_int_or_inf(value, "nself") # must be int or inf
        check_is_gteq(value, "nself", 0)    # must be >= 0; cannot be negative
        self._nself = value
        
    @property
    def gmapfn(self) -> GeneticMapFunction:
        """Description for property gmapfn."""
        return self._gmapfn
    @gmapfn.getter
    def gmapfn(self) -> GeneticMapFunction:
        """Get data for property gmapfn."""
        return self._gmapfn
    @gmapfn.setter
    def gmapfn(self, value: GeneticMapFunction) -> None:
        """Set data for property gmapfn."""
        check_isinstance(value, "gmapfn", GeneticMapFunction)
        self._gmapfn = value
    
    @property
    def mem(self) -> int:
        """Description for property mem."""
        return self._mem
    @mem.getter
    def mem(self) -> int:
        """Get data for property mem."""
        return self._mem
    @mem.setter
    def mem(self, value: int) -> None:
        """Set data for property mem."""
        check_is_int(value, "mem")
        check_is_gt(value, "mem", 0)
        self._mem = value
    
    @property
    def unique_parents(self) -> bool:
        """Description for property unique_parents."""
        return self._unique_parents
    @unique_parents.getter
    def unique_parents(self) -> bool:
        """Get data for property unique_parents."""
        return self._unique_parents
    @unique_parents.setter
    def unique_parents(self, value: bool) -> None:
        """Set data for property unique_parents."""
        check_is_bool(value, "unique_parents")
        self._unique_parents = value
    
    @property
    def method(self) -> str:
        """Description for property method."""
        return self._method
    @method.getter
    def method(self) -> str:
        """Get data for property method."""
        return self._method
    @method.setter
    def method(self, value: str) -> None:
        """Set data for property method."""
        check_is_str(value, "method")       # must be string
        value = value.lower()               # convert to lowercase
        options = ("single", "pareto")      # method options
        # if not method supported raise ValueError
        if value not in options:
            raise ValueError("Unsupported method. Options are: " + ", ".join(map(str, options)))
        self._method = value
    
    @property
    def weight(self) -> Union[numpy.ndarray,Callable,str]:
        """Description for property weight."""
        return self._weight
    @weight.getter
    def weight(self) -> Union[numpy.ndarray,Callable,str]:
        """Get data for property weight."""
        return self._weight
    @weight.setter
    def weight(self, value: Union[numpy.ndarray,Callable,str]) -> None:
        """Set data for property weight."""
        check_isinstance(value, "weight", (numpy.ndarray,Callable,str))
        if isinstance(value, str):
            # convert to lowercase
            value = value.lower()
            # convert string to function
            if value == 'magnitude':
                value = weight_absolute
            elif value == 'equal':
                value = weight_one
            else:
                raise ValueError("Unsupported weight. Options are: 'magnitude', 'equal'")
        self._weight = value
    
    @property
    def target(self) -> Union[numpy.ndarray,Callable,str]:
        """Description for property target."""
        return self._target
    @target.getter
    def target(self) -> Union[numpy.ndarray,Callable,str]:
        """Get data for property target."""
        return self._target
    @target.setter
    def target(self, value: Union[numpy.ndarray,Callable,str]) -> None:
        """Set data for property target."""
        check_isinstance(value, "target", (numpy.ndarray,Callable,str))
        if isinstance(value, str):
            # convert to lowercase
            value = value.lower()
            # convert string to function
            if value == 'positive':
                value = target_positive
            elif value == 'negative':
                value = target_negative
            elif value == 'stabilizing':
                value = target_stabilizing
            else:
                raise ValueError("Unsupported weight. Options are: 'positive', 'negative', 'stabilizing'")
        self._target = value
    
    @property
    def objfn_trans(self) -> Union[Callable,None]:
        """Description for property objfn_trans."""
        return self._objfn_trans
    @objfn_trans.getter
    def objfn_trans(self) -> Union[Callable,None]:
        """Get data for property objfn_trans."""
        return self._objfn_trans
    @objfn_trans.setter
    def objfn_trans(self, value: Union[Callable,None]) -> None:
        """Set data for property objfn_trans."""
        if value is not None:                       # if given object
            check_is_callable(value, "objfn_trans") # must be callable
        self._objfn_trans = value
    
    @property
    def objfn_trans_kwargs(self) -> dict:
        """Description for property objfn_trans_kwargs."""
        return self._objfn_trans_kwargs
    @objfn_trans_kwargs.getter
    def objfn_trans_kwargs(self) -> dict:
        """Get data for property objfn_trans_kwargs."""
        return self._objfn_trans_kwargs
    @objfn_trans_kwargs.setter
    def objfn_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set data for property objfn_trans_kwargs."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "objfn_trans_kwargs")  # check is dict
        self._objfn_trans_kwargs = value
    
    @property
    def objfn_wt(self) -> numpy.ndarray:
        """Description for property objfn_wt."""
        return self._objfn_wt
    @objfn_wt.getter
    def objfn_wt(self) -> numpy.ndarray:
        """Get data for property objfn_wt."""
        return self._objfn_wt
    @objfn_wt.setter
    def objfn_wt(self, value: numpy.ndarray) -> None:
        """Set data for property objfn_wt."""
        self._objfn_wt = value
    
    @property
    def ndset_trans(self) -> Union[Callable,None]:
        """Description for property ndset_trans."""
        return self._ndset_trans
    @ndset_trans.getter
    def ndset_trans(self) -> Union[Callable,None]:
        """Get data for property ndset_trans."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable,None]) -> None:
        """Set data for property ndset_trans."""
        if value is not None:                       # if given object
            check_is_callable(value, "ndset_trans") # must be callable
        self._ndset_trans = value
    
    @property
    def ndset_trans_kwargs(self) -> dict:
        """Description for property ndset_trans_kwargs."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.getter
    def ndset_trans_kwargs(self) -> dict:
        """Get data for property ndset_trans_kwargs."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.setter
    def ndset_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set data for property ndset_trans_kwargs."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "ndset_trans_kwargs")  # check is dict
        self._ndset_trans_kwargs = value
    
    @property
    def ndset_wt(self) -> numpy.ndarray:
        """Description for property ndset_wt."""
        return self._ndset_wt
    @ndset_wt.getter
    def ndset_wt(self) -> numpy.ndarray:
        """Get data for property ndset_wt."""
        return self._ndset_wt
    @ndset_wt.setter
    def ndset_wt(self, value: numpy.ndarray) -> None:
        """Set data for property ndset_wt."""
        self._ndset_wt = value
    
    @property
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Description for property rng."""
        return self._rng
    @rng.getter
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Get data for property rng."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[numpy.random.Generator,numpy.random.RandomState,None]) -> None:
        """Set data for property rng."""
        if value is None:       # if None
            value = global_prng # use default random number generator
        check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
        self._rng = value
    
    @property
    def soalgo(self) -> UnconstrainedOptimizationAlgorithm:
        """Description for property soalgo."""
        return self._soalgo
    @soalgo.getter
    def soalgo(self) -> UnconstrainedOptimizationAlgorithm:
        """Get data for property soalgo."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[UnconstrainedOptimizationAlgorithm,None]) -> None:
        """Set data for property soalgo."""
        # if value is None, use a default hillclimber
        if value is None:
            value = UnconstrainedSteepestAscentSetHillClimber(rng = self.rng)
        check_is_OptimizationAlgorithm(value, "soalgo")
        self._soalgo = value
    
    @property
    def moalgo(self) -> UnconstrainedOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.getter
    def moalgo(self) -> UnconstrainedOptimizationAlgorithm:
        """Get data for property moalgo."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[UnconstrainedOptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        # if value is None, use a default nsga-ii algorithm
        if value is None:
            value = UnconstrainedNSGA2SetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                mu = 100,       # number of parents in population
                lamb = 100,     # number of progeny to produce
                M = 1.5,        # algorithm crossover genetic map length
                rng = self.rng  # PRNG source
            )
        check_is_OptimizationAlgorithm(value, "moalgo")
        self._moalgo = value
    
    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
    def _calc_mkrwt(self, gpmod: AdditiveLinearGenomicModel) -> numpy.ndarray:
        if callable(self.weight):
            return self.weight(gpmod.u_a)
        elif isinstance(self.weight, numpy.ndarray):
            return self.weight
        else:
            raise TypeError("variable 'weight' must be a callable function or numpy.ndarray")
    
    def _calc_tfreq(self, gpmod: AdditiveLinearGenomicModel) -> numpy.ndarray:
        if callable(self.target):
            return self.target(gpmod.u_a)
        elif isinstance(self.target, numpy.ndarray):
            return self.target
        else:
            raise TypeError("variable 'target' must be a callable function or numpy.ndarray")

    def _calc_tminor(self, tfreq: numpy.ndarray) -> numpy.ndarray:
        return (tfreq == 0.0)
    
    def _calc_tmajor(self, tfreq: numpy.ndarray) -> numpy.ndarray:
        return (tfreq == 1.0)
    
    def _calc_thet(self, tminor: numpy.ndarray, tmajor: numpy.ndarray):
        return numpy.logical_not(numpy.logical_or(tminor, tmajor))

    def _calc_xmap(self, ntaxa: int) -> numpy.ndarray:
        """
        Calculate the cross map.

        Parameters
        ----------
        ntaxa : int
            Number of taxa.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(s,d)`` containing cross map indices.

            Where:

            - ``s`` is the number of elements in the upper triangle, including
              or not including the diagonal (depending on ``unique_parents``).
            - ``d`` is the number of parents in the cross.
        """
        if self.unique_parents:         # if we want unique parents
            return numpy.array(         # create a numpy.ndarray
                list(                   # convert to list
                    triudix(            # generator for indices without diagonal
                        ntaxa,          # number of taxa
                        self.nparent    # number of parents
                    )
                )
            )
        else:                           # otherwise we don't want unique parents
            return numpy.array(         # create a numpy.ndarray
                list(                   # convert to list
                    triuix(             # generator for indices with diagonal
                        ntaxa,          # number of taxa
                        self.nparent    # number of parents
                    )
                )
            )

    def _calc_vmat(
            self, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix
        ) -> numpy.ndarray:
        # calculate variance matrix
        vmat = self.vmatfcty.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            ncross = self.ncross, 
            nprogeny = self.nprogeny, 
            nself = self.nself, 
            gmapfn = self.gmapfn, 
            mem = self.mem
        )
        return vmat.mat

    ############################## Object Methods ##############################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
        """
        Select individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes (unphased most likely)
        ptdf : pandas.DataFrame
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
        # get selection parameters
        nparent = self.nparent
        ncross = self.ncross
        nprogeny = self.nprogeny
        objfn_wt = self.objfn_wt
        ndset_trans = self.ndset_trans
        ndset_trans_kwargs = self.ndset_trans_kwargs
        ndset_wt = self.ndset_wt
        method = self.method

        # single-objective method: objfn_trans returns a single value for each
        # selection configuration
        if method == "single":
            # get number of taxa
            ntaxa = pgmat.ntaxa

            # get vectorized objective function
            objfn = self.objfn(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                **kwargs
            )

            # optimize using hill-climber algorithm
            opt = self.soalgo.optimize(
                k = nparent,                    # number of parents to select
                sspace = numpy.arange(ntaxa),   # parental indices
                rng = self.rng,                 # PRNG source
                objwt = objfn_wt                # maximizing function
            )

            # get best solution
            sel = opt["soln"]

            # add optimization details to miscellaneous output
            if miscout is not None:     # if miscout was provided
                miscout.update(opt)     # add dict to dict

            return pgmat, sel, ncross, nprogeny

        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif method == "pareto":
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
                nparent = nparent,
                objfn_trans = self.objfn_trans,
                objfn_trans_kwargs = self.objfn_trans_kwargs,
                objfn_wt = objfn_wt,
                weight = self.weight,
                target = self.target
            )

            # get scores for each of the points along the pareto frontier
            score = ndset_wt * ndset_trans(frontier, **ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to miscout
            if miscout is not None:
                miscout["frontier"] = frontier
                miscout["sel_config"] = sel_config

            return pgmat, sel_config[ix], ncross, nprogeny

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a selection objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Phased genotype matrix.
        gmat : GenotypeMatrix
            Input genotype matrix.
        ptdf : pandas.DataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Not used by this function.
        gpmod : AdditiveLinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A selection objective function for the specified problem.
        """
        # get selection parameters
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # calculate default function parameters
        xmap = self._calc_xmap(gmat.ntaxa)
        mat = gmat.mat                      # (n,p) get genotype matrix
        ploidy = gmat.ploidy                # (scalar) get number of phases
        mkrwt = self._calc_mkrwt(gpmod)     # (p,t) get marker weights
        tfreq = self._calc_tfreq(gpmod)     # (p,t) get target allele frequencies
        tminor = self._calc_tminor(tfreq)
        tmajor = self._calc_tmajor(tfreq)
        thet = self._calc_thet(tminor, tmajor)
        vmat = self._calc_vmat(gpmod, pgmat)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,         # byte code pointer
            self.objfn_static.__globals__,      # global variables
            None,                               # new name for the function
            (xmap, mat, ploidy, mkrwt, tfreq, 
             tminor, thet, tmajor, vmat, trans, 
             trans_kwargs),                     # default values for arguments
            self.objfn_static.__closure__       # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a vectorized selection objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Not used by this function.
        gmat : GenotypeMatrix
            Input genotype matrix.
        ptdf : pandas.DataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Not used by this function.
        gpmod : AdditiveLinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A vectorized selection objective function for the specified problem.
        """
        # get selection parameters
        weight = self.weight
        target = self.target

        # calculate default function parameters
        mat = gmat.mat                      # (n,p) get genotype matrix
        ntaxa = pgmat.ntaxa                 # get number of taxa
        ploidy = gmat.ploidy                # (scalar) get number of phases
        u = gpmod.u_a                       # (p,t) get regression coefficients
        xmap = self._calc_xmap(ntaxa)       # (s,p) get the cross map
        mkrwt = self._calc_mkrwt(weight, u) # (p,t) get marker weights
        tfreq = self._calc_tfreq(target, u) # (p,t) get target allele frequencies
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # generate variance matrix
        if AdditiveGeneticVarianceMatrix in self.vmatcls.__mro__:
            vmat = self.vmatcls.from_algmod(
                algmod = gpmod,
                pgmat = pgmat,
                ncross = self.ncross,
                nprogeny = self.nprogeny,
                s = self.nself,
                gmapfn = self.gmapfn,
                mem = self.mem
            )
        elif AdditiveGenicVarianceMatrix in self.vmatcls.__mro__:
            vmat = self.vmatcls.from_algmod(
                algmod = gpmod,
                pgmat = pgmat,
                nprogeny = self.nprogeny,
                mem = self.mem
            )

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (xmap, mat, ploidy, tfreq, mkrwt,
            vmat, trans, trans_kwargs),         # default values for arguments
            self.objfn_vec_static.__closure__   # closure byte code pointer
        )

        return outfn

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
        """
        Calculate a Pareto frontier for objectives.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : pandas.DataFrame
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
        # process inputs, apply defaults as needed.
        nparent = self.nparent
        objfn_wt = self.objfn_wt

        # get number of taxa
        ntaxa = gmat.ntaxa

        # create objective function
        objfn = self.objfn(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max,
            **kwargs
        )

        # use multi-objective optimization to approximate Pareto front.
        frontier, sel_config, misc = self.moalgo.optimize(
            objfn = objfn,                  # objective function
            k = nparent,                    # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = objfn_wt             # weights to apply to each objective
        )

        # handle miscellaneous output
        if miscout is not None:     # if miscout is provided
            miscout.update(misc)    # add 'misc' to 'miscout', overwriting as needed

        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(
            xsel: numpy.ndarray, 
            xmap: numpy.ndarray, 
            mat: numpy.ndarray, 
            ploidy: int, 
            mkrwt: numpy.ndarray, 
            tfreq: numpy.ndarray, 
            tminor: numpy.ndarray, 
            thet: numpy.ndarray, 
            tmajor: numpy.ndarray, 
            vmat: numpy.ndarray, 
            trans: Callable, 
            kwargs: dict
        ) -> numpy.ndarray:
        """
        Multi-objective genomic mating objective function.

        - The goal is to minimize all objectives for this function.
        - This is a bare bones function. Minimal error checking is done.

        Objectives: :math:`F(\\textbf{x})`

        .. math::

            F(\\textbf{x}) = {[f^{\\textup{PAU}}(\\textbf{x}), f^{\\textup{PAFD}}(\\textbf{x})]}'

        Population Allele Unavailability (PAU): :math:`f^{\\textup{PAU}}(\\textbf{x})`

        Formal PAU definition:

        .. math::

            f^{\\textup{PAU}}(\\textbf{x}) = \\textbf{w} \\cdot \\textbf{u}

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency.
        From the selection allele frequencies and the target allele frequencies
        ``tfreq``, determine if the target frequencies can be attained after
        unlimited generations of selection. If the target allele frequency at a
        locus cannot be attained, score locus as ``1``, otherwise score as
        ``0``. Store this into a binary score vector :math:`\\textbf{u}`.
        Take the dot product between the binary score vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAU}}(\\textbf{x})` and return the result.

        Population Allele Frequency Distance (PAFD): :math:`f^{\\textup{PAFD}}(\\textbf{x})`

        Formal PAFD definition:

        .. math::
            f^{\\textup{PAFD}}(\\textbf{x}) = \\textbf{w} \\cdot \\left | \\textbf{p}_{x} - \\textbf{p}_{t} \\right |

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency
        :math:`\\textbf{p}_{x}`. From the selection allele frequencies and the
        target allele frequencies :math:`\\textbf{p}_{t} =` ``tfreq``,
        calculate the absolute value of the difference between the two vectors.
        Finally, take the dot product between the difference vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAFD}}(\\textbf{x})` and return the result.

        Sum of Progeny Standard Deviations of Additive Variance (SPstdA): :math:`f^{\\textup{SPstdA}}(\\textbf{x})`

        Formal SPstdA definition:

        .. math::

            f^{\\textup{SPstdA}}(\\textbf{x}) = \\sum_{c \\in S} \\sigma_{A,c}

        Given a progeny variance matrix :math:`\\Sigma_{A} =` ``vmat`` and a
        selection indices vector :math:`\\textbf{x} =` ``sel``, take the sum of
        the square root of the progeny variance
        :math:`\\sigma_{A,c} = \\sqrt{\\Sigma_{A,c}}` for each cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A cross selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of crosses to select.

            Each index indicates which cross specified by ``xmap`` to select.
        xmap : numpy.ndarray
            A cross selection index map array of shape ``(s,d)``.

            Where:

            - ``s`` is the size of the sample space (number of cross
              combinations for ``d`` parents).
            - ``d`` is the number of parents.
        mat : numpy.ndarray
            A genotype matrix of shape ``(n,p)`` representing only biallelic
            loci. One of the two alleles at a locus is coded using a ``1``. The
            other allele is coded as a ``0``. ``mat`` holds the counts of the
            allele coded by ``1``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.

            Example::

                # matrix of shape (n = 3, p = 4)
                mat = numpy.array([[0,2,1,0],
                                   [2,2,1,1],
                                   [0,1,0,2]])
        ploidy : int
            Number of phases that the genotype matrix ``mat`` represents.
        tfreq : floating, numpy.ndarray
            A target allele frequency matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Example::

                tfreq = numpy.array([0.2, 0.6, 0.7, 0.5])
        mkrwt : numpy.ndarray
            A marker weight coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Remarks:

            - All values in ``mkrwt`` must be non-negative.
        vmat : numpy.ndarray, Matrix
            A variance matrix of shape ``(n,...,n,t)``. Can be a
            ``numpy.ndarray`` or a Matrix of some sort. Must be have the ``[]``
            operator to access elements of the matrix.

            Where:

            - ``n`` is the number of parental candidates.
            - ``t`` is the number of traits.
            - ``(n,...,n,t)`` is a tuple of length ``d + 1``.
            - ``d`` is the number of parents for a cross.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single numpy.ndarray argument.
            - Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        mogm : numpy.ndarray
            A MOGM score matrix of shape ``(t + t + t,)`` if ``trans`` is
            ``None``. Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.

            Matrix element ordering for un-transformed MOGM score matrix:

            - The first set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` PAU outputs for each trait.
            - The second set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` PAFD outputs for each trait.
            - The third set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` SPstdA outputs for each trait.
        """
        # get selected individuals from cross selection indices
        sel = xmap[xsel].flatten()

        # calculate the allele frequency of the selected subset
        # (n,p)[(k,),:,None] -> (p,1)
        pfreq = (1.0 / (ploidy * len(sel))) * mat[sel,:,None].sum(0)

        # determine where allele frequencies are < 1.0
        # (p,1)
        p_ltmajor = (pfreq < 1.0)

        # determine where allele frequencies are > 0.0
        # (p,1)
        p_gtminor = (pfreq > 0.0)

        # determine where allele frequencies are < 1.0 and > 0.0
        # (p,1)
        p_het = numpy.logical_and(p_ltmajor, p_gtminor)

        # determine where alleles are unavailable using precomputed arrays
        # (p,t)
        allele_unavail = numpy.logical_not(
            numpy.logical_or(
                numpy.logical_and(p_ltmajor, tminor), 
                numpy.logical_or(
                    numpy.logical_and(p_het, thet), 
                    numpy.logical_and(p_gtminor, tmajor)
                )
            )
        )

        # get selection tuple for the variance matrix
        # ([...],...,[...],:)
        vsel = tuple(xmap[xsel].T) + (slice(None),)

        # calculate the manhattan distance and PAFD
        # (p,t) -> (t,)
        pafd = (mkrwt * numpy.absolute(tfreq - pfreq)).sum(0)
        
        # calculate the allele unavailability
        # (p,t) -> (t,)
        pau = (mkrwt * allele_unavail).sum(0)

        # calculate the sum of standard deviations of additive variance
        spstda = -numpy.sqrt(vmat[vsel]).sum(0)

        # concatenate to make MOGS vector
        # (t,) and (t,) and (t,) -> (t + t + t,)
        mogs = numpy.concatenate([pau, pafd, spstda])

        # apply transformations
        if trans is not None:
            mogs = trans(mogs, **kwargs)

        return mogs

    @staticmethod
    def objfn_static(sel, xmap, mat, ploidy, tfreq, mkrwt, vmat, trans, kwargs):
        """
        Multi-objective genomic mating objective function.

        - The goal is to minimize all objectives for this function.
        - This is a bare bones function. Minimal error checking is done.

        Objectives: :math:`F(\\textbf{x})`

        .. math::

            F(\\textbf{x}) = {[f^{\\textup{PAU}}(\\textbf{x}), f^{\\textup{PAFD}}(\\textbf{x})]}'

        Population Allele Unavailability (PAU): :math:`f^{\\textup{PAU}}(\\textbf{x})`

        Formal PAU definition:

        .. math::

            f^{\\textup{PAU}}(\\textbf{x}) = \\textbf{w} \\cdot \\textbf{u}

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency.
        From the selection allele frequencies and the target allele frequencies
        ``tfreq``, determine if the target frequencies can be attained after
        unlimited generations of selection. If the target allele frequency at a
        locus cannot be attained, score locus as ``1``, otherwise score as
        ``0``. Store this into a binary score vector :math:`\\textbf{u}`.
        Take the dot product between the binary score vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAU}}(\\textbf{x})` and return the result.

        Population Allele Frequency Distance (PAFD): :math:`f^{\\textup{PAFD}}(\\textbf{x})`

        Formal PAFD definition:

        .. math::
            f^{\\textup{PAFD}}(\\textbf{x}) = \\textbf{w} \\cdot \\left | \\textbf{p}_{x} - \\textbf{p}_{t} \\right |

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency
        :math:`\\textbf{p}_{x}`. From the selection allele frequencies and the
        target allele frequencies :math:`\\textbf{p}_{t} =` ``tfreq``,
        calculate the absolute value of the difference between the two vectors.
        Finally, take the dot product between the difference vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAFD}}(\\textbf{x})` and return the result.

        Sum of Progeny Standard Deviations of Additive Variance (SPstdA): :math:`f^{\\textup{SPstdA}}(\\textbf{x})`

        Formal SPstdA definition:

        .. math::

            f^{\\textup{SPstdA}}(\\textbf{x}) = \\sum_{c \\in S} \\sigma_{A,c}

        Given a progeny variance matrix :math:`\\Sigma_{A} =` ``vmat`` and a
        selection indices vector :math:`\\textbf{x} =` ``sel``, take the sum of
        the square root of the progeny variance
        :math:`\\sigma_{A,c} = \\sqrt{\\Sigma_{A,c}}` for each cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A cross selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of crosses to select.

            Each index indicates which cross specified by ``xmap`` to select.
        xmap : numpy.ndarray
            A cross selection index map array of shape ``(s,d)``.

            Where:

            - ``s`` is the size of the sample space (number of cross
              combinations for ``d`` parents).
            - ``d`` is the number of parents.
        mat : numpy.ndarray
            A genotype matrix of shape ``(n,p)`` representing only biallelic
            loci. One of the two alleles at a locus is coded using a ``1``. The
            other allele is coded as a ``0``. ``mat`` holds the counts of the
            allele coded by ``1``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.

            Example::

                # matrix of shape (n = 3, p = 4)
                mat = numpy.array([[0,2,1,0],
                                   [2,2,1,1],
                                   [0,1,0,2]])
        ploidy : int
            Number of phases that the genotype matrix ``mat`` represents.
        tfreq : floating, numpy.ndarray
            A target allele frequency matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Example::

                tfreq = numpy.array([0.2, 0.6, 0.7, 0.5])
        mkrwt : numpy.ndarray
            A marker weight coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Remarks:

            - All values in ``mkrwt`` must be non-negative.
        vmat : numpy.ndarray, Matrix
            A variance matrix of shape ``(n,...,n,t)``. Can be a
            ``numpy.ndarray`` or a Matrix of some sort. Must be have the ``[]``
            operator to access elements of the matrix.

            Where:

            - ``n`` is the number of parental candidates.
            - ``t`` is the number of traits.
            - ``(n,...,n,t)`` is a tuple of length ``d + 1``.
            - ``d`` is the number of parents for a cross.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single numpy.ndarray argument.
            - Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        mogm : numpy.ndarray
            A MOGM score matrix of shape ``(t + t + t,)`` if ``trans`` is
            ``None``. Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.

            Matrix element ordering for un-transformed MOGM score matrix:

            - The first set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` PAU outputs for each trait.
            - The second set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` PAFD outputs for each trait.
            - The third set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` SPstdA outputs for each trait.
        """
        # get cross configurations
        # (s,d)[(k,),:] -> (k,d)
        sel = xmap[sel,:]

        # flatten cross selections for PAU and PAFD calculations
        # (k,d) -> (kd,)
        fsel = sel.ravel()

        ####################################################
        ######### PAU and PAFD shared calculations #########
        ####################################################

        # generate a view of the genotype matrix that only contains 'sel' rows.
        # (n,p)[(kd,),:] -> (kd,p)
        sgeno = mat[fsel,:]

        # calculate reciprocal number of phases
        # ploidy * number of individuals in 'sgeno'
        rphase = 1.0 / (ploidy * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        # (k,p).sum(0) -> (p,)
        # (p,) * scalar -> (p,)
        # (p,None) -> (p,1)
        # Remark: we need (p,1) for broadcasting with (p,t) arrays
        pfreq = (sgeno.sum(0) * rphase)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

        # calculate allele unavailability
        # (p,t)
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # then set True if sel has allele freq == 0
            numpy.where(            # else
                tfreq > 0.0,        # if 0.0 < target freq < 1.0
                numpy.logical_or(   # then set True if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,   # mask for whether population freq <= 0.0
                    pfreq_gteq_1    # mask for whether population freq >= 1.0
                ),
                pfreq_gteq_1        # else set True if pop freq is >= 1.0
            )
        )

        # calculate distance between target and population
        # (p,t)-(p,1) -> (p,t)
        dist = numpy.absolute(tfreq - pfreq)

        # compute f_PAU(x)
        # (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        pau = (mkrwt * allele_unavail).sum(0)

        # compute f_PAFD(x)
        # (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        pafd = (mkrwt * dist).sum(0)

        ####################################################
        ############### SPstdA calculations ################
        ####################################################

        # contruct a matrix element selection tuple
        # (k,d).T -> (d,k)
        # (d,k) --transform--> ((k,),...,(k,)) (d elements in outer tuple)
        # ((k,),...,(k,)) + (:) -> ((k,),...,(k,),:)
        vmatsel = tuple(tuple(e) for e in sel.T) + (slice(None),)

        # select variance elements
        # (n,...,n,t)[((k,),...,(k,),:)] -> (k,t)
        velem = vmat[vmatsel]

        # compute -f_SPstdA
        # sqrt((k,t)) -> (k,t)
        # (k,t).sum(0) -> (t,)
        # -1 * (t,) -> (t,)
        spstda = -numpy.sqrt(velem).sum(0)

        ####################################################
        ######### output preparation calculations ##########
        ####################################################

        # concatenate to make MOGM matrix
        # (t,) cat (t,) cat (t,) -> (t + t + t,)
        mogm = numpy.concatenate([pau, pafd, spstda])

        # apply transformations
        if trans is not None:
            mogm = trans(mogm, **kwargs)

        return mogm

    @staticmethod
    def objfn_vec_static(sel, xmap, mat, ploidy, tfreq, mkrwt, vmat, trans, kwargs):
        """
        A vectorized multi-objective genomic selection objective function.

        - The goal is to minimize all objectives for this function.
        - This is a bare bones function. Minimal error checking is done.

        Objectives: :math:`F(\\textbf{x})`

        .. math::

            F(\\textbf{x}) = {[f^{\\textup{PAU}}(\\textbf{x}), f^{\\textup{PAFD}}(\\textbf{x})]}'

        Population Allele Unavailability (PAU): :math:`f^{\\textup{PAU}}(\\textbf{x})`

        .. math::

            f^{\\textup{PAU}}(\\textbf{x}) = \\textbf{w} \\cdot \\textbf{u}

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency.
        From the selection allele frequencies and the target allele frequencies
        ``tfreq``, determine if the target frequencies can be attained after
        unlimited generations of selection. If the target allele frequency at a
        locus cannot be attained, score locus as ``1``, otherwise score as
        ``0``. Store this into a binary score vector :math:`\\textbf{u}`.
        Take the dot product between the binary score vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAU}}(\\textbf{x})` and return the result.

        Population Allele Frequency Distance (PAFD): :math:`f^{\\textup{PAFD}}(\\textbf{x})`

        .. math::
            f^{\\textup{PAFD}}(\\textbf{x}) = \\textbf{w} \\cdot \\left | \\textbf{p}_{x} - \\textbf{p}_{t} \\right |

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency
        :math:`\\textbf{p}_{x}`. From the selection allele frequencies and the
        target allele frequencies :math:`\\textbf{p}_{t} =` ``tfreq``,
        calculate the absolute value of the difference between the two vectors.
        Finally, take the dot product between the difference vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAFD}}(\\textbf{x})` and return the result.

        Sum of Progeny Standard Deviations of Additive Variance (SPstdA): :math:`f^{\\textup{SPstdA}}(\\textbf{x})`

        Formal SPstdA definition:

        .. math::

            f^{\\textup{SPstdA}}(\\textbf{x}) = \\sum_{c \\in S} \\sigma_{A,c}

        Given a progeny variance matrix :math:`\\Sigma_{A} =` ``vmat`` and a
        selection indices vector :math:`\\textbf{x} =` ``sel``, take the sum of
        the square root of the progeny variance
        :math:`\\sigma_{A,c} = \\sqrt{\\Sigma_{A,c}}` for each cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of configurations to score.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            ``sel`` cannot be ``None``.
        xmap : numpy.ndarray
            A cross selection index map array of shape ``(s,d)``.

            Where:

            - ``s`` is the size of the sample space (number of cross
              combinations for ``d`` parents).
            - ``d`` is the number of parents.
        mat : numpy.ndarray
            A genotype matrix of shape ``(n,p)`` representing only biallelic
            loci. One of the two alleles at a locus is coded using a ``1``. The
            other allele is coded as a ``0``. ``mat`` holds the counts of the
            allele coded by ``1``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.

            Example::

                # matrix of shape (n = 3, p = 4)
                mat = numpy.array([[0,2,1,0],
                                   [2,2,1,1],
                                   [0,1,0,2]])
        ploidy : int
            Number of phases that the genotype matrix ``mat`` represents.
        tfreq : floating, numpy.ndarray
            A target allele frequency matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Example::

                tfreq = numpy.array([0.2, 0.6, 0.7, 0.5])
        mkrwt : numpy.ndarray
            A marker weight coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Remarks:

            - All values in ``mkrwt`` must be non-negative.
        vmat : numpy.ndarray, Matrix
            A variance matrix of shape ``(n,...,n,t)``. Can be a
            ``numpy.ndarray`` or a Matrix of some sort. Must be have the ``[]``
            operator to access elements of the matrix.

            Where:

            - ``n`` is the number of parental candidates.
            - ``t`` is the number of traits.
            - ``(n,...,n,t)`` is a tuple of length ``d + 1``.
            - ``d`` is the number of parents for a cross.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        mogm : numpy.ndarray
            A MOGM score matrix of shape ``(j,t + t + t)`` if ``trans`` is
            ``None``. Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.

            Matrix element ordering for un-transformed MOGM score matrix:

            - The first set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` PAU outputs for each trait.
            - The second set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` PAFD outputs for each trait.
            - The third set of ``t`` elements in the ``mogm`` output correspond
              to the ``t`` SPstdA outputs for each trait.
        """
        # get cross configurations
        # (s,d)[(j,k,),:] -> (j,k,d)
        sel = xmap[sel,:]

        # flatten cross selections for PAU and PAFD calculations
        # (j,k,d) -> (j,kd)
        fsel = numpy.empty(                                 # create empty array
            (sel.shape[0], sel.shape[1] * sel.shape[2]),    # (j,kd)
            dtype = sel.dtype                               # array data type
        )
        for i in range(sel.shape[0]):   # for each row
            fsel[i] = sel[i].ravel()    # copy column information

        ####################################################
        ######### PAU and PAFD shared calculations #########
        ####################################################

        # generate a view of the genotype matrix that only contains 'sel' rows.
        # (n,p)[(j,kd),:] -> (j,kd,p)
        sgeno = mat[fsel,:]

        # calculate reciprocal number of phases
        # ploidy * number of individuals in 'sgeno'
        rphase = 1.0 / (ploidy * sgeno.shape[2])

        # calculate population frequencies; add axis for correct broadcast
        # (j,kd,p).sum(1) -> (j,p)
        # (j,p) * scalar -> (j,p)
        # (j,p)[:,None] -> (j,p,1)
        # Remark: we need (j,p,1) for broadcasting with (p,t) arrays
        pfreq = (sgeno.sum(1) * rphase)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # (j,p,1) is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # (j,p,1) is population freq >= 1.0

        # calculate allele unavailability
        # (j,p,t)
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # (p,t) if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # (j,p,1) then set True if sel has allele freq == 0
            numpy.where(            # (j,p,t) else
                tfreq > 0.0,        # (p,t) if 0.0 < target freq < 1.0
                numpy.logical_or(   # (j,p,1) then set True if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,   # (j,p,1) mask for whether population freq <= 0.0
                    pfreq_gteq_1    # (j,p,1) mask for whether population freq >= 1.0
                ),
                pfreq_gteq_1        # (j,p,1) else set True if pop freq is >= 1.0
            )
        )

        # calculate distance between target and population
        # (p,t)-(j,p,1) -> (j,p,t)
        dist = numpy.absolute(tfreq - pfreq)

        # compute f_PAU(x)
        # (p,t) * (j,p,t) -> (j,p,t)
        # (j,p,t).sum[1] -> (j,t)
        pau = (mkrwt * allele_unavail).sum(1)

        # compute f_PAFD(x)
        # (p,t) * (j,p,t) -> (j,p,t)
        # (j,p,t).sum[1] -> (j,t)
        pafd = (mkrwt * dist).sum(1)

        ####################################################
        ############### SPstdA calculations ################
        ####################################################

        # contruct a matrix element selection tuple
        # (j,k,d) --transform--> ( ((k,),...,(k,)), ..., ((k,),...,(k,)) )
        # (((k,),...,(k,)),...,((k,),...,(k,))) + (:) -> (((k,),...,(k,)),...,((k,),...,(k,)),:)
        vmatsel = tuple(tuple(tuple(e) for e in x.T) for x in sel) + (slice(None),)

        # select variance elements
        # (n,...,n,t)[(((k,),...,(k,)),...,((k,),...,(k,)),:)] -> (j,k,t)
        velem = vmat[vmatsel]

        # compute -f_SPstdA
        # sqrt((j,k,t)) -> (j,k,t)
        # (k,t).sum(1) -> (k,t)
        # -1 * (k,t) -> (k,t)
        spstda = -numpy.sqrt(velem).sum(1)

        ####################################################
        ######### output preparation calculations ##########
        ####################################################

        # concatenate to make MOGM matrix
        # (j,t) cat (j,t) cat (j,t) -> (j,t + t + t)
        mogm = numpy.concatenate([pau, pafd, spstda], axis = 1)

        # apply transformations
        if trans is not None:
            mogm = trans(mogm, **kwargs)

        return mogm
