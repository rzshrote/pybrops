"""
Module implementing selection protocols for Optimal Population Value selection.
"""

__all__ = [
    "OptimalPopulationValueBaseSelection",
    "OptimalPopulationValueSubsetSelection"
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.SubsetSelectionProtocol import SubsetSelectionProtocol
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.OptimalPopulationValueSelectionProblem import OptimalPopulationValueSubsetSelectionProblem
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class OptimalPopulationValueSelectionMixin(metaclass=ABCMeta):
    """
    Semi-abstract class for Optimal Population Value (OPV) Selection with constraints.
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
    def nhaploblk(self) -> Integral:
        """Number of haplotype blocks to consider."""
        return self._nhaploblk
    @nhaploblk.setter
    def nhaploblk(self, value: Integral) -> None:
        """Set number of haplotype blocks to consider."""
        check_is_Integral(value, "nhaploblk")
        check_is_gt(value, "nhaploblk", 0)
        self._nhaploblk = value

class OptimalPopulationValueSubsetSelection(OptimalPopulationValueSelectionMixin,SubsetSelectionProtocol):
    """
    Class defining Optimal Haploid Value (OHV) Selection for subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nhaploblk: Integral,
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
        # make assignments from SubsetSelectionProtocol second
        super(OptimalPopulationValueSubsetSelection, self).__init__(
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
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_AdditiveLinearGenomicModel(gpmod, "gpmod")
        
        # get decision space parameters
        ntaxa = pgmat.ntaxa
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nparent)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nparent)

        # construct problem
        prob = OptimalPopulationValueSubsetSelectionProblem.from_pgmat_gpmod(
            nhaploblk = self.nhaploblk,
            pgmat = pgmat,
            gpmod = gpmod,
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

    ################ Single Objective Solve ################
    # inherit sosolve() from SubsetSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from SubsetSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from SubsetSelectionProtocol
