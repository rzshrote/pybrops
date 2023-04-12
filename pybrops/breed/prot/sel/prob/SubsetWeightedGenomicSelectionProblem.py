"""
Module implementing weighted genomic selection as a subset optimization problem.
"""

from numbers import Integral, Number
from typing import Callable, Optional, Union
import numpy
from pybrops.breed.prot.sel.prob.SubsetGeneralizedWeightedGenomicSelectionProblem import SubsetGeneralizedWeightedGenomicSelectionProblem


class SubsetWeightedGenomicSelectionProblem(SubsetGeneralizedWeightedGenomicSelectionProblem):
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
            obj_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
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
        super(SubsetWeightedGenomicSelectionProblem, self).__init__(
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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    # override function, not because it is different, but because it has a different docstring
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Weighted Genomic Selection 
        (WGS). Scoring for WGS is defined as the sum of Weighted Genomic 
        Estimated Breeding Values (wGEBV) for a population.

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
            A wGEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,t)[(ndecn,),:] -> (ndecn,t)    # select individuals
        # Step 2: (ndecn,t).sum(0)  -> (t,)         # sum across all individuals
        return self._wgebv[x,:].sum(0)
