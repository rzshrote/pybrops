"""
Module implementing Optimal Population Value (OPV) Selection optimization problems.
"""

__all__ = [
    "OptimalPopulationValueSubsetSelectionProblem"
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
# from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
# from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.util.haplo import haplobin, haplobin_bounds, nhaploblk_chrom
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class OptimalPopulationValueSelectionProblem(SelectionProblem,metaclass=ABCMeta):
    """Helper class to implement properties common to OPV."""

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            haplomat: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
            **kwargs: dict
        ) -> None:
        """
        Constructor for OptimalPopulationValueSelectionProblem.
        
        Parameters
        ----------
        haplomat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,h,t)``.

            Where:

            - ``m`` is the number of chromosome phases (2 for diploid, etc.).
            - ``n`` is the number of individuals.
            - ``h`` is the number of haplotype blocks.
            - ``t`` is the number of traits.
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a upper limit for the decision variables.
        nobj: Integral
            Number of objectives.
        obj_wt: numpy.ndarray
            Objective function weights.
        obj_trans: Callable, None
            A transformation function transforming a latent space vector to an objective space vector.
            The transformation function must be of the form: ``obj_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the identity transformation function: copy the latent space vector to the objective space vector.
        obj_trans_kwargs: dict, None
            Keyword arguments for the latent space to objective space transformation function.
            If None, an empty dictionary is used.
        nineqcv: Integral,
            Number of inequality constraints.
        ineqcv_wt: numpy.ndarray,
            Inequality constraint violation weights.
        ineqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an inequality constraint violation vector.
            The transformation function must be of the form: ``ineqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        ineqcv_trans_kwargs: Optional[dict],
            Keyword arguments for the latent space to inequality constraint violation space transformation function.
            If None, an empty dictionary is used.
        neqcv: Integral
            Number of equality constraints.
        eqcv_wt: numpy.ndarray
            Equality constraint violation weights.
        eqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an equality constraint violation vector.
            The transformation function must be of the form: ``eqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        eqcv_trans_kwargs: dict, None
            Keyword arguments for the latent space to equality constraint violation space transformation function.
            If None, an empty dictionary is used.
        kwargs : dict
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(OptimalPopulationValueSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
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
            **kwargs
        )
        # assignments
        self.haplomat = haplomat

    ############################ Object Properties #############################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in haplotype matrix
        return self._haplomat.shape[3]

    ################### Haplotype matrix ###################
    @property
    def haplomat(self) -> numpy.ndarray:
        """Haplotype effect matrix of shape ``(m,n,h,t)``."""
        return self._haplomat
    @haplomat.setter
    def haplomat(self, value: numpy.ndarray) -> None:
        """Set haplotype effect matrix."""
        check_is_ndarray(value, "haplomat")
        check_ndarray_ndim(value, "haplomat", 4)
        self._haplomat = value
    
    ##################### Ploidy level #####################
    @property
    def ploidy(self) -> Integral:
        """ploidy."""
        return self._haplomat.shape[0]

    ######################### Private Object Methods ###########################
    @staticmethod
    def _calc_haplomat(
            pgmat: PhasedGenotypeMatrix, 
            algpmod: AdditiveLinearGenomicModel, 
            nhaploblk: Integral
        ) -> numpy.ndarray:
        """
        Calculate a haplotype matrix from a genome matrix and model.

        Parameters
        ----------
        gmat : PhasedGenotypeMatrix
            A genome matrix.
        mod : DenseAdditiveLinearGenomicModel
            A genomic prediction model.

        Returns
        -------
        hmat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,b,t)``.
        """
        mat         = pgmat.mat              # get genotypes
        genpos      = pgmat.vrnt_genpos      # get genetic positions
        chrgrp_stix = pgmat.vrnt_chrgrp_stix # get chromosome start indices
        chrgrp_spix = pgmat.vrnt_chrgrp_spix # get chromosome stop indices
        chrgrp_len  = pgmat.vrnt_chrgrp_len  # get chromosome marker lengths
        u           = algpmod.u_a              # get regression coefficients

        if (chrgrp_stix is None) or (chrgrp_spix is None):
            raise RuntimeError("markers are not sorted by chromosome position")

        # get number of chromosomes
        nchr = len(chrgrp_stix)

        if nhaploblk < nchr:
            raise RuntimeError("number of haplotype blocks is less than the number of chromosomes")

        # calculate number of marker blocks to assign to each chromosome
        nblk = nhaploblk_chrom(nhaploblk, genpos, chrgrp_stix, chrgrp_spix)

        # ensure there are enough markers per chromosome
        if numpy.any(nblk > chrgrp_len):
            raise RuntimeError(
                "number of haplotype blocks assigned to a chromosome greater than number of available markers"
            )

        # calculate haplotype bins
        hbin = haplobin(nblk, genpos, chrgrp_stix, chrgrp_spix)

        # define shape
        # (m,n,b,t)
        s = (mat.shape[0], mat.shape[1], nhaploblk, u.shape[1])

        # allocate haplotype matrix
        # (m,n,b,t)
        hmat = numpy.empty(s, dtype = u.dtype)

        # get boundary indices
        hstix, hspix, hlen = haplobin_bounds(hbin)

        # OPTIMIZE: perhaps eliminate one loop using dot function
        # fill haplotype matrix
        for i in range(hmat.shape[3]):                          # for each trait
            for j,(st,sp) in enumerate(zip(hstix,hspix)):       # for each haplotype block
                hmat[:,:,j,i] = mat[:,:,st:sp].dot(u[st:sp,i])  # take dot product and fill

        return hmat

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod(
            cls,
            nhaploblk: Integral,
            pgmat: PhasedGenotypeMatrix,
            gpmod: AdditiveLinearGenomicModel,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
            **kwargs: dict
        ) -> "OptimalPopulationValueSelectionProblem":
        # type checks
        check_is_Integral(nhaploblk, "nhaploblk")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_AdditiveLinearGenomicModel(gpmod, "gpmod")
        
        # calculate estimated breeding values and relationships
        haplomat = cls._calc_haplomat(pgmat, gpmod, nhaploblk)

        # construct class
        out = cls(
            haplomat = haplomat,
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
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
            **kwargs
        )

        return out

class OptimalPopulationValueSubsetSelectionProblem(SubsetSelectionProblem,OptimalPopulationValueSelectionProblem):
    """
    docstring for SubsetOptimalPopulationValueSelectionProblem.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            haplomat: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetOptimalPopulationValueSelectionProblem.
        
        Parameters
        ----------
        haplomat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,h,t)``.

            Where:

            - ``m`` is the number of chromosome phases (2 for diploid, etc.).
            - ``n`` is the number of individuals.
            - ``h`` is the number of haplotype blocks.
            - ``t`` is the number of traits.
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a upper limit for the decision variables.
        nobj: Integral
            Number of objectives.
        obj_wt: numpy.ndarray
            Objective function weights.
        obj_trans: Callable, None
            A transformation function transforming a latent space vector to an objective space vector.
            The transformation function must be of the form: ``obj_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the identity transformation function: copy the latent space vector to the objective space vector.
        obj_trans_kwargs: dict, None
            Keyword arguments for the latent space to objective space transformation function.
            If None, an empty dictionary is used.
        nineqcv: Integral,
            Number of inequality constraints.
        ineqcv_wt: numpy.ndarray,
            Inequality constraint violation weights.
        ineqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an inequality constraint violation vector.
            The transformation function must be of the form: ``ineqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        ineqcv_trans_kwargs: Optional[dict],
            Keyword arguments for the latent space to inequality constraint violation space transformation function.
            If None, an empty dictionary is used.
        neqcv: Integral
            Number of equality constraints.
        eqcv_wt: numpy.ndarray
            Equality constraint violation weights.
        eqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an equality constraint violation vector.
            The transformation function must be of the form: ``eqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        eqcv_trans_kwargs: dict, None
            Keyword arguments for the latent space to equality constraint violation space transformation function.
            If None, an empty dictionary is used.
        kwargs : dict
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(OptimalPopulationValueSubsetSelectionProblem, self).__init__(
            haplomat = haplomat,
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
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
            **kwargs
        )
        # assignments
        self.haplomat = haplomat

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Optimal Population Value
        Selection.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An OPV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # get max haplotype value
        # (m,n,h,t)[:,(k,),:,:] -> (m,k,h,t)
        # (m,k/2,2,h,t).max((0,1)) -> (h,t)
        # (h,t).sum(0) -> (t,)
        # scalar * (t,) -> (t,)
        out = -self.ploidy * self._haplomat[:,x,:,:].max((0,1)).sum(0)

        return out

# need better interpretation of the Real scenario
# class OptimalPopulationValueRealSelectionProblem(RealSelectionProblem,OPVSProblemProperties):
#     """
#     docstring for OptimalPopulationValueRealSelectionProblem.
#     """
#     pass

# need better interpretation of the Integer scenario
# class OptimalPopulationValueIntegerSelectionProblem(IntegerSelectionProblem,OPVSProblemProperties):
#     """
#     docstring for OptimalPopulationValueIntegerSelectionProblem.
#     """
#     pass

# need better interpretation of the Binary scenario
# class OptimalPopulationValueBinarySelectionProblem(BinarySelectionProblem,OPVSProblemProperties):
#     """
#     docstring for OptimalPopulationValueBinarySelectionProblem.
#     """
#     pass
