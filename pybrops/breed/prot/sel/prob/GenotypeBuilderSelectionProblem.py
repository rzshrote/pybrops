"""
Module implementing Genotype Builder (GB) Selection optimization problems.
"""

__all__ = [
    "GenotypeBuilderSubsetSelectionProblem"
]

from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gt, check_is_lt


class GenotypeBuilderSelectionProblem(SelectionProblem):
    """Helper class to implement properties common to GB."""

    ############################ Object Properties #############################

    ############## Number of latent variables ##############
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in haplotype matrix
        return self._haplomat.shape[3]

    ################### Haplotype matrix ###################
    @property
    def haplomat(self) -> numpy.ndarray:
        """Haplotype effect matrix of shape ``(m,n,b,t)``."""
        return self._haplomat
    @haplomat.setter
    def haplomat(self, value: numpy.ndarray) -> None:
        """Set haplotype effect matrix."""
        check_is_ndarray(value, "gebv")
        check_ndarray_ndim(value, "gebv", 4)
        self._haplomat = value

    ##################### Ploidy level #####################
    @property
    def ploidy(self) -> Integral:
        """ploidy."""
        return self._haplomat.shape[0]

    ########## Number of Best Founders to Select ###########
    @property
    def nbestfndr(self) -> Integral:
        """nbestfndr."""
        return self._nbestfndr
    @nbestfndr.setter
    def nbestfndr(self, value: Integral) -> None:
        """Set nbestfndr."""
        check_is_Integral(value, "nbestfndr")   # must be int
        check_is_gt(value, "nbestfndr", 0)      # int must be >0
        check_is_lt(value, "nbestfndr", self._haplomat.shape[1]) # int must be < number of candidates
        self._nbestfndr = value

class GenotypeBuilderSubsetSelectionProblem(SubsetSelectionProblem,GenotypeBuilderSelectionProblem):
    """
    Class representing Genotype Builder (GB) Selection problems in subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            haplomat: numpy.ndarray,
            nbestfndr: Integral,
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
        Constructor for SubsetGenotypeBuilderSelectionProblem.
        
        Parameters
        ----------
        haplomat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,b,t)``.

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
        super(GenotypeBuilderSubsetSelectionProblem, self).__init__(
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
        self.nbestfndr = nbestfndr

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Genotype Builder
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
            An GB matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # get best haplotype within each individual
        # (m,n,h,t)[:,(k,),:,:] -> (m,k,h,t)
        # (m,k,h,t).max(0) -> (k,h,t)
        bestphase = self._haplomat[:,x,:,:].max(0)

        # sort the bestphase tensor along the individual axis to get best individuals
        # (k,h,t) 
        bestphase.sort(0)

        # get top haplotype index boundaries
        k = len(x)
        st = k - self.nbestfndr
        sp = k

        # get top haplotypes and sum across the individual and haplotype axes
        # scale by nphase / nbestfndr
        # (k,h,t)[(nbestfndr,),:,:] -> (nbestfndr,h,t)
        # (nbestfndr,h,t).sum((0,1)) -> (t,)
        # scalar * (t,) -> (t,)
        out = -(self.ploidy / self.nbestfndr) * bestphase[st:sp,:,:].sum((0,1))

        return out

# need better interpretation for the Real scenario
# class GenotypeBuilderRealSelectionProblem(RealSelectionProblem,GBSProblemProperties):
#     """
#     Class representing Genotype Builder (GB) Selectionproblems in real search spaces.
#     """
#     pass

# need better interpretation for the Integer scenario
# class GenotypeBuilderIntegerSelectionProblem(IntegerSelectionProblem,GBSProblemProperties):
#     """
#     Class representing Genotype Builder (GB) Selectionproblems in integer search spaces.
#     """
#     pass

# need better interpretation for the Binary scenario
# class GenotypeBuilderBinarySelectionProblem(BinarySelectionProblem,GBSProblemProperties):
#     """
#     Class representing Genotype Builder (GB) Selectionproblems in binary search spaces.
#     """
#     pass
