"""
Module implementing weighted genomic selection as a subset optimization problem.
"""

__all__ = [
    "WeightedGenomicBinarySelectionProblem",
    "WeightedGenomicIntegerSelectionProblem",
    "WeightedGenomicRealSelectionProblem",
    "WeightedGenomicSubsetSelectionProblem"
]

from numbers import Integral, Number
from typing import Callable, Optional, Union
import numpy
from pybrops.breed.prot.sel.prob.GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem import GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem, GeneralizedWeightedGenomicEstimatedBreedingValueIntegerSelectionProblem, GeneralizedWeightedGenomicEstimatedBreedingValueRealSelectionProblem, GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem


class WeightedGenomicSubsetSelectionProblem(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem):
    """
    docstring for SubsetWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: numpy.ndarray,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            eqcv_trans_kwargs: Optional[dict],
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call SubsetGeneralizedWeightedGenomicSelectionProblem constructor
        super(WeightedGenomicSubsetSelectionProblem, self).__init__(
            Z_a = Z_a,
            u_a = u_a,
            fafreq = fafreq,
            alpha = 0.5,
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

class WeightedGenomicRealSelectionProblem(GeneralizedWeightedGenomicEstimatedBreedingValueRealSelectionProblem):
    """
    docstring for RealWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: numpy.ndarray,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            eqcv_trans_kwargs: Optional[dict],
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call RealGeneralizedWeightedGenomicSelectionProblem constructor
        super(WeightedGenomicRealSelectionProblem, self).__init__(
            Z_a = Z_a,
            u_a = u_a,
            fafreq = fafreq,
            alpha = 0.5,
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

class WeightedGenomicIntegerSelectionProblem(GeneralizedWeightedGenomicEstimatedBreedingValueIntegerSelectionProblem):
    """
    docstring for IntegerWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: numpy.ndarray,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            eqcv_trans_kwargs: Optional[dict],
            **kwargs: dict
        ) -> None:
        """
        Constructor for IntegerWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call IntegerGeneralizedWeightedGenomicSelectionProblem constructor
        super(WeightedGenomicIntegerSelectionProblem, self).__init__(
            Z_a = Z_a,
            u_a = u_a,
            fafreq = fafreq,
            alpha = 0.5,
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

class WeightedGenomicBinarySelectionProblem(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem):
    """
    docstring for BinaryWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: numpy.ndarray,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]],
            eqcv_trans_kwargs: Optional[dict],
            **kwargs: dict
        ) -> None:
        """
        Constructor for BinaryWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call BinaryGeneralizedWeightedGenomicSelectionProblem constructor
        super(WeightedGenomicBinarySelectionProblem, self).__init__(
            Z_a = Z_a,
            u_a = u_a,
            fafreq = fafreq,
            alpha = 0.5,
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

