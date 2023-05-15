"""
Module implementing generalized weighted genomic selection protocols.
"""

__all__ = [
    "WeightedGenomicBinarySelection",
    "WeightedGenomicIntegerSelection",
    "WeightedGenomicRealSelection",
    "WeightedGenomicSubsetSelection"
]

from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.GeneralizedWeightedGenomicEstimatedBreedingValueSelection import GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelection, GeneralizedWeightedGenomicEstimatedBreedingValueIntegerSelection, GeneralizedWeightedGenomicEstimatedBreedingValueRealSelection, GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelection
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm

class WeightedGenomicSubsetSelection(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelection):
    """
    Weighted Genomic Selection in a subset search space.
    """

    ########################## Special Object Methods ##########################
    # should NOT be an abstract method
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
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
        Constructor for WeightedGenomicSubsetSelection.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(WeightedGenomicSubsetSelection, self).__init__(
            nparent = nparent, 
            ncross = ncross, 
            nprogeny = nprogeny,
            alpha = 0.5,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

class WeightedGenomicRealSelection(GeneralizedWeightedGenomicEstimatedBreedingValueRealSelection):
    """
    Weighted Genomic Selection in a real search space.
    """

    ########################## Special Object Methods ##########################
    # should NOT be an abstract method
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
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
        Constructor for WeightedGenomicRealSelection.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(WeightedGenomicRealSelection, self).__init__(
            nparent = nparent, 
            ncross = ncross, 
            nprogeny = nprogeny,
            alpha = 0.5,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

class WeightedGenomicIntegerSelection(GeneralizedWeightedGenomicEstimatedBreedingValueIntegerSelection):
    """
    Weighted Genomic Selection in an integer search space.
    """

    ########################## Special Object Methods ##########################
    # should NOT be an abstract method
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
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
        Constructor for WeightedGenomicIntegerSelection.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(WeightedGenomicIntegerSelection, self).__init__(
            nparent = nparent, 
            ncross = ncross, 
            nprogeny = nprogeny,
            alpha = 0.5,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

class WeightedGenomicBinarySelection(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelection):
    """
    Weighted Genomic Selection in a subset search space.
    """

    ########################## Special Object Methods ##########################
    # should NOT be an abstract method
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
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
        Constructor for WeightedGenomicBinarySelection.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(WeightedGenomicBinarySelection, self).__init__(
            nparent = nparent, 
            ncross = ncross, 
            nprogeny = nprogeny,
            alpha = 0.5,
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )
