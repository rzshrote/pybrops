"""
Module implementing selection protocols for maximum mean expected heterozygosity selection.
"""

__all__ = [
    "OptimalHaploidValueSelectionMixin",
    "OptimalHaploidValueBinarySelection",
    "OptimalHaploidValueIntegerSelection",
    "OptimalHaploidValueRealSelection",
    "OptimalHaploidValueSubsetSelection"
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
import pandas

from pybrops.breed.prot.sel.prob.OptimalHaploidValueSelectionProblem import OptimalHaploidValueBinarySelectionProblem, OptimalHaploidValueIntegerSelectionProblem, OptimalHaploidValueRealSelectionProblem, OptimalHaploidValueSubsetSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.SubsetSelectionProtocol import SubsetSelectionProtocol
from pybrops.breed.prot.sel.RealSelectionProtocol import RealSelectionProtocol
from pybrops.breed.prot.sel.IntegerSelectionProtocol import IntegerSelectionProtocol
from pybrops.breed.prot.sel.BinarySelectionProtocol import BinarySelectionProtocol
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.core.error.error_type_python import check_is_Integral, check_is_bool
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class OptimalHaploidValueSelectionMixin(metaclass=ABCMeta):
    """
    Semi-abstract class for Optimal Haploid Value (OHV) Selection with constraints.
    """

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def ntrait(self) -> Integral:
        """Number of traits to expect from matrix inputs."""
        return self._ntrait
    @ntrait.setter
    def ntrait(self, value: Integral) -> None:
        """Set number of traits to expect."""
        check_is_Integral(value, "ntrait")
        check_is_gt(value, "ntrait", 0)
        self._ntrait = value

    @property
    def nhaploblk(self) -> Integral:
        """Number of haplotype blocks to consider."""
        return self._nhaploblk
    @nhaploblk.setter
    def nhaploblk(self, value: Integral) -> None:
        """Set number of haplotype blocks to consider."""
        check_is_Integral(value, "nhaploblk")
        check_is_gt(value, "nhaploblk", 0)
        self._nhaploblk = value
    
    @property
    def unique_parents(self) -> bool:
        """Whether parents should be unique."""
        return self._unique_parents
    @unique_parents.setter
    def unique_parents(self, value: bool) -> None:
        """Set whether parents should be unique."""
        check_is_bool(value, "unique_parents")
        self._unique_parents = value

class OptimalHaploidValueSubsetSelection(OptimalHaploidValueSelectionMixin,SubsetSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nhaploblk: Integral,
            unique_parents: bool,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[OptimizationAlgorithm] = None,
            moalgo: Optional[OptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nhaploblk = nhaploblk
        self.unique_parents = unique_parents
        # make assignments from SubsetSelectionProtocol second
        super(OptimalHaploidValueSubsetSelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map (inefficient)
        xmap = OptimalHaploidValueSubsetSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space = numpy.arange(len(xmap))
        decn_space_lower = numpy.repeat(0, self.ncross)
        decn_space_upper = numpy.repeat(len(xmap)-1, self.ncross)

        # construct problem
        prob = OptimalHaploidValueSubsetSelectionProblem.from_pgmat_gpmod(
            nparent = self.nparent,
            nhaploblk = self.nhaploblk,
            unique_parents = self.unique_parents,
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = self.ncross,
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

class OptimalHaploidValueRealSelection(OptimalHaploidValueSelectionMixin,RealSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for real search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nhaploblk: Integral,
            unique_parents: bool,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[OptimizationAlgorithm] = None,
            moalgo: Optional[OptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nhaploblk = nhaploblk
        self.unique_parents = unique_parents
        # make assignments from RealSelectionProtocol second
        super(OptimalHaploidValueRealSelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map (inefficient)
        xmap = OptimalHaploidValueRealSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, len(xmap))
        decn_space_upper = numpy.repeat(1.0, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = OptimalHaploidValueRealSelectionProblem.from_pgmat_gpmod(
            nparent = self.nparent,
            nhaploblk = self.nhaploblk,
            unique_parents = self.unique_parents,
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = len(xmap),
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

class OptimalHaploidValueIntegerSelection(OptimalHaploidValueSelectionMixin,IntegerSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for a integer search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nhaploblk: Integral,
            unique_parents: bool,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[OptimizationAlgorithm] = None,
            moalgo: Optional[OptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nhaploblk = nhaploblk
        self.unique_parents = unique_parents
        # make assignments from IntegerSelectionProtocol second
        super(OptimalHaploidValueIntegerSelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map (inefficient)
        xmap = OptimalHaploidValueIntegerSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(numpy.sum(self.nmating * self.nprogeny), len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = OptimalHaploidValueIntegerSelectionProblem.from_pgmat_gpmod(
            nparent = self.nparent,
            nhaploblk = self.nhaploblk,
            unique_parents = self.unique_parents,
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = len(xmap),
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

class OptimalHaploidValueBinarySelection(OptimalHaploidValueSelectionMixin,BinarySelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for a binary search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nhaploblk: Integral,
            unique_parents: bool,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[OptimizationAlgorithm] = None,
            moalgo: Optional[OptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nhaploblk = nhaploblk
        self.unique_parents = unique_parents
        # make assignments from BinarySelectionProtocol second
        super(OptimalHaploidValueBinarySelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map (inefficient)
        xmap = OptimalHaploidValueBinarySelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(1, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = OptimalHaploidValueBinarySelectionProblem.from_pgmat_gpmod(
            nparent = self.nparent,
            nhaploblk = self.nhaploblk,
            unique_parents = self.unique_parents,
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = len(xmap),
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
