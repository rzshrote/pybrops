"""
Module implementing selection protocols for Usefulness Criterion (UC) selection.
"""

__all__ = [
    "UsefulnessCriterionSelectionMixin",
    "UsefulnessCriterionBinarySelection",
    "UsefulnessCriterionIntegerSelection",
    "UsefulnessCriterionRealSelection",
    "UsefulnessCriterionSubsetSelection"
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Optional, Union
from typing import Callable

import numpy
from numpy.random import Generator, RandomState
import pandas
import scipy.stats

from pybrops.breed.prot.sel.BinaryMateSelectionProtocol import BinaryMateSelectionProtocol
from pybrops.breed.prot.sel.IntegerMateSelectionProtocol import IntegerMateSelectionProtocol
from pybrops.breed.prot.sel.RealMateSelectionProtocol import RealMateSelectionProtocol
from pybrops.breed.prot.sel.SubsetMateSelectionProtocol import SubsetMateSelectionProtocol
from pybrops.breed.prot.sel.prob.BinaryMateSelectionProblem import BinaryMateSelectionProblem
from pybrops.breed.prot.sel.prob.IntegerMateSelectionProblem import IntegerMateSelectionProblem
from pybrops.breed.prot.sel.prob.RealMateSelectionProblem import RealMateSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetMateSelectionProblem import SubsetMateSelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionBinarySelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionIntegerSelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionRealSelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionSubsetSelectionProblem
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory, check_is_GeneticVarianceMatrixFactory
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.core.error.error_type_python import check_is_Integral, check_is_Real, check_is_bool
from pybrops.core.error.error_value_python import check_is_gt, check_is_gteq, check_is_in_interval_exclusive
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction, check_is_GeneticMapFunction
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class UsefulnessCriterionSelectionMixin(metaclass=ABCMeta):
    """
    Semi-abstract class for Optimal Haploid Value (OHV) Selection with constraints.
    """

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def ntrait(self) -> Integral:
        """Number of traits to expect from EBV matrix inputs."""
        return self._ntrait
    @ntrait.setter
    def ntrait(self, value: Integral) -> None:
        """Set number of traits to expect."""
        check_is_Integral(value, "ntrait")
        check_is_gt(value, "ntrait", 0)
        self._ntrait = value

    @property
    def nself(self) -> Integral:
        """Description for property nself."""
        return self._nself
    @nself.setter
    def nself(self, value: Integral) -> None:
        """Set data for property nself."""
        check_is_Integral(value, "nself")     # must be int
        check_is_gteq(value, "nself", 0)   # int must be >=0
        self._nself = value

    @property
    def upper_percentile(self) -> Real:
        """Description for property upper_percentile."""
        return self._upper_percentile
    @upper_percentile.setter
    def upper_percentile(self, value: Real) -> None:
        """Set data for property upper_percentile."""
        check_is_Real(value, "upper_percentile")  # must be a number
        check_is_in_interval_exclusive(value, "upper_percentile", 0.0, 1.0)
        self._upper_percentile = value
    
    @property
    def selection_intensity(self) -> Real:
        """Get selection intensity."""
        return scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - self._upper_percentile)) / self._upper_percentile

    @property
    def vmatfcty(self) -> GeneticVarianceMatrixFactory:
        """Description for property vmatfcty."""
        return self._vmatfcty
    @vmatfcty.setter
    def vmatfcty(self, value: GeneticVarianceMatrixFactory) -> None:
        """Set data for property vmatfcty."""
        check_is_GeneticVarianceMatrixFactory(value, "vmatfcty")
        self._vmatfcty = value
    
    @property
    def gmapfn(self) -> GeneticMapFunction:
        """Description for property gmapfn."""
        return self._gmapfn
    @gmapfn.setter
    def gmapfn(self, value: GeneticMapFunction) -> None:
        """Set data for property gmapfn."""
        check_is_GeneticMapFunction(value, "gmapfn")
        self._gmapfn = value

    @property
    def unique_parents(self) -> bool:
        """Whether parents should be unique."""
        return self._unique_parents
    @unique_parents.setter
    def unique_parents(self, value: bool) -> None:
        """Set whether parents should be unique."""
        check_is_bool(value, "unique_parents")
        self._unique_parents = value

class UsefulnessCriterionSubsetSelection(UsefulnessCriterionSelectionMixin,SubsetMateSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
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
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from SubsetMateSelectionProtocol second
        super(UsefulnessCriterionSubsetSelection, self).__init__(
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
        ) -> SubsetMateSelectionProblem:
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
        
        # get the cross map
        xmap = UsefulnessCriterionSubsetSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space = numpy.arange(len(xmap))
        decn_space_lower = numpy.repeat(0, self.ncross)
        decn_space_upper = numpy.repeat(len(xmap)-1, self.ncross)

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionSubsetSelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median,
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
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

    ################ Single Objective Solve ################
    # inherit sosolve() from SubsetMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from SubsetMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from SubsetMateSelectionProtocol

class UsefulnessCriterionRealSelection(UsefulnessCriterionSelectionMixin,RealMateSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for real search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
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
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from RealMateSelectionProtocol second
        super(UsefulnessCriterionRealSelection, self).__init__(
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
        ) -> RealMateSelectionProblem:
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
        xmap = UsefulnessCriterionRealSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, len(xmap))
        decn_space_upper = numpy.repeat(1.0, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionRealSelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median, 
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
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

    ################ Single Objective Solve ################
    # inherit sosolve() from RealMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from RealMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from RealMateSelectionProtocol

class UsefulnessCriterionIntegerSelection(UsefulnessCriterionSelectionMixin,IntegerMateSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for a integer search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
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
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from IntegerMateSelectionProtocol second
        super(UsefulnessCriterionIntegerSelection, self).__init__(
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
        ) -> IntegerMateSelectionProblem:
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
        xmap = UsefulnessCriterionIntegerSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(self.ncross * self.nparent * self.nmating, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionIntegerSelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median, 
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
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

    ################ Single Objective Solve ################
    # inherit sosolve() from IntegerMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from IntegerMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from IntegerMateSelectionProtocol

class UsefulnessCriterionBinarySelection(UsefulnessCriterionSelectionMixin,BinaryMateSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for a binary search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
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
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from BinaryMateSelectionProtocol second
        super(UsefulnessCriterionBinarySelection, self).__init__(
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
        ) -> BinaryMateSelectionProblem:
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
        xmap = UsefulnessCriterionBinarySelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(1, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionBinarySelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median, 
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
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

    ################ Single Objective Solve ################
    # inherit sosolve() from BinaryMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from BinaryMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from BinaryMateSelectionProtocol
