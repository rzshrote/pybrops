"""
Module defining optimization problems for binary optimal constribution selection.
"""

from ast import Tuple
from numbers import Integral
from typing import Callable

import numpy
from pybrops.breed.prot.sel.prob.DenseSubsetSelectionProblem import DenseSubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_is_2d, check_ndarray_is_square

class BinaryOptimalContributionSelectionProblem(DenseSubsetSelectionProblem):
    """
    docstring for BinaryOptimalContributionSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            bv: numpy.ndarray,
            C: numpy.ndarray,
            ndecn: Integral,
            decn_space: numpy.ndarray,
            encode_trans: Callable[[numpy.ndarray,dict],Tuple],
            encode_trans_kwargs: dict,
            nobj: Integral,
            obj_wt: numpy.ndarray,
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for BinaryOptimalContributionSelectionProblem.
        
        Parameters
        ----------
        bv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        C : numpy.ndarray
            An upper triangle matrix of shape ``(n,n)`` resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

            Where:

            - ``n`` is the number of individuals.
        ndecn : Integral
            Number of decision variables.
        decn_space : numpy.ndarray
            A 1d array containing the set of available elements
        encode_trans : Callable[[numpy.ndarray],Tuple]
            Function which transforms outputs from ``latentfn`` to a tuple ``(obj,ineqcv,eqcv)``.
        encode_trans_kwargs : dict
            ``latentfn`` output transformation function keyword arguments.
        nobj : Integral
            Number of objectives.
        obj_wt : numpy.ndarray
            Objective function weights.
        nineqcv : Integral
            Number of inequality constraint violation functions.
        ineqcv_wt : numpy.ndarray
            Inequality constraint violation function weights.
        neqcv : Integral
            Number of equality constraint violation functions.
        eqcv_wt : numpy.ndarray
            Equality constraint violation function weights.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # do not call super constructor
        # order dependent assignments
        self.bv = bv
        self.C = C
        self.ndecn = ndecn
        self.decn_space = decn_space
        self.encode_trans = encode_trans
        self.encode_trans_kwargs = encode_trans_kwargs
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def bv(self) -> numpy.ndarray:
        """Breeding value matrix."""
        return self._bv
    @bv.setter
    def bv(self, value: numpy.ndarray) -> None:
        """Set breeding value matrix."""
        check_is_ndarray(value, "bv")
        check_ndarray_is_2d(value, "bv")
        self._bv = value
    @bv.deleter
    def bv(self) -> None:
        """Delete breeding value matrix."""
        del self._bv

    @property
    def C(self) -> numpy.ndarray:
        """Cholesky decomposition of the kinship matrix."""
        return self._C
    @C.setter
    def C(self, value: numpy.ndarray) -> None:
        """Set Cholesky decomposition of the kinship matrix."""
        check_is_ndarray(value, "C")
        check_ndarray_is_2d(value, "C")
        check_ndarray_is_square(value, "C")
        self._C = value
    @C.deleter
    def C(self) -> None:
        """Delete Cholesky decomposition of the kinship matrix."""
        del self._C

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Encode a candidate solution for the given Problem into an ``l`` 
        dimensional latent evaluation space.
        
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
            A matrix of shape (1+t,).

            The first index in the array is the mean genomic relationship 
            (a minimizing objective):

            .. math::
                MGR = || \\textbf{C} \\textbf{(sel)} ||_2

            The next `t` indices in the array are the sum of breeding values for 
            each of ``t`` traits for the selection (all maximizing objectives).

            Where:

            - ``t`` is the number of traits.
        """
        # calculate MEH
        # (n,n)[:,(k,)] -> (n,k)
        # scalar * (n,k).sum(1) -> (n,)
        Cx = (1.0 / len(x)) * self.C[:,x].sum(1)

        # calculate mean genomic relationship
        # norm2( (n,), keepdims=True ) -> (1,)
        mgr = numpy.linalg.norm(Cx, ord = 2, keepdims = True)

        # calculate breeding value of the selection
        # (n,t)[(k,),:] -> (k,t)
        # (k,t).sum(0) -> (t,)
        gain = self.bv[x,:].sum(0)
        
        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgr,gain])

        # return (1+t,)
        return out

    def evalfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> Tuple:
        """
        Evaluate a candidate solution for the given Problem.
        
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
        out : tuple
            A tuple ``(obj, ineqcv, eqcv)``.
            
            Where:
            
            - ``obj`` is a numpy.ndarray of shape ``(nobj,)`` that contains 
                objective function evaluations.
            - ``ineqcv`` is a numpy.ndarray of shape ``(nineqcv,)`` that contains 
                inequality constraint violation values.
            - ``eqcv`` is a numpy.ndarray of shape ``(neqcv,)`` that contains 
                equality constraint violation values.
        """
        # transform and return tuple 
        return self.encode_trans(
            self.latentfn(x, *args, **kwargs), 
            **self.encode_trans_kwargs
        )

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_BinaryOptimalContributionSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type BinaryOptimalContributionSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BinaryOptimalContributionSelectionProblem):
        raise TypeError("'{0}' must be of type BinaryOptimalContributionSelectionProblem.".format(vname))
