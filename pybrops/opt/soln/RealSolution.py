"""
Partial implementation of the SetSolution interface.
"""

# list of all public imports in the module
__all__ = [
    "RealSolution",
    "check_is_RealSolution",
]

# imports
from numbers import Integral, Real
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_real
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_shape_eq
from pybrops.opt.soln.Solution import Solution

class RealSolution(Solution):
    """
    Class for optimization problem solutions with real decision variables.
    """

    ########################## Special Object Methods ##########################

    # implementation of abstract method
    def __init__(
            self,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            nobj: Integral,
            obj_wt: Union[numpy.ndarray,Real,None],
            nineqcv: Union[Integral,None],
            ineqcv_wt: Union[numpy.ndarray,Real,None],
            neqcv: Union[Integral,None],
            eqcv_wt: Union[numpy.ndarray,Real,None],
            nsoln: Integral,
            soln_decn: numpy.ndarray,
            soln_obj: numpy.ndarray,
            soln_ineqcv: Union[numpy.ndarray,None],
            soln_eqcv: Union[numpy.ndarray,None],
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealSolution.
        
        Parameters
        ----------
        ndecn : Integral
            Number of decision variables.

        decn_space : numpy.ndarray
            A 1d array containing the set of available elements
        
        decn_space_lower : numpy.ndarray
            A 1d array representing the lower bound of the decision variables.
        
        decn_space_upper : numpy.ndarray
            A 1d array representing the upper bound of the decision variables.
        
        nobj : Integral
            Number of objectives for the problem.
        
        obj_wt : numpy.ndarray, Real, None
            Objective function weights. Weights from this vector are applied 
            to objective function values via the Hadamard product. If values 
            are ``1.0`` or ``-1.0``, this can be used to specify minimizing 
            and maximizing objectives, respectively.

            If ``obj_wt`` is ``numpy.ndarray``, then the array must be of shape 
            ``(nobj,)``.

            If ``obj_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``.

            If ``obj_wt`` is ``None``, then the value ``1.0`` is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``. This assumes that all 
            objectives are to be minimized.
        
        nineqcv : Integral, None
            Number of inequality constraint violation functions. This is 
            equivalent to the vector length returned by the ``ineqcv_trans`` 
            function. Must be ``Integral`` greater than or equal to zero.

            If ``nineqcv`` is ``None``, then set to zero.

        ineqcv_wt : numpy.ndarray, None
            Inequality constraint violation function weights. Weights from this 
            vector are applied to inequality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``ineqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(nineqcv,)``. This assumes that 
            all constraints are to be minimized.

        neqcv : Integral, None
            Number of equality constraint violations. This is equivalent to the 
            vector length returned by the ``eqcv_trans`` function. Must be 
            ``Integral`` greater than or equal to zero.
        
            If ``neqcv`` is ``None``, then set to zero.

        eqcv_wt : numpy.ndarray, None
            Equality constraint violation function weights. Weights from this 
            vector are applied to equality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``eqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(neqcv,)``. This assumes that 
            all constraints are to be minimized.

        nsoln : Integral
            Number of solutions.
        
        soln_decn : numpy.ndarray
            Solution decision vectors for each solution.
            This matrix is of shape ``(nsoln,ndecn)``.
        
        soln_obj : numpy.ndarray
            Solution objective vectors for each solution.
            This matrix is of shape ``(nsoln,nobj)``.
        
        soln_ineqcv : numpy.ndarray
            Solution inequality constraint violiation vectors for each solution.
            This matrix is of shape ``(nsoln,nineqcv)``

        soln_eqcv : numpy.ndarray
            Solution equality constraint violiation vectors for each solution.
            This matrix is of shape ``(nsoln,neqcv)``

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        self.ndecn = ndecn
        self.decn_space = decn_space
        self.decn_space_lower = decn_space_lower
        self.decn_space_upper = decn_space_upper
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt
        self.nsoln = nsoln
        self.soln_decn = soln_decn
        self.soln_obj = soln_obj
        self.soln_ineqcv = soln_ineqcv
        self.soln_eqcv = soln_eqcv

    ############################ Object Properties #############################
    
    # override decn_space setter properties
    @Solution.decn_space.setter
    def decn_space(self, value: Union[numpy.ndarray,None]) -> None:
        """Set decision space boundaries."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_real(value, "decn_space")
            check_ndarray_shape_eq(value, "decn_space", (2,self.ndecn))
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space' must be of type numpy.ndarray or None")
        self._decn_space = value

    # override decn_space setter properties
    @Solution.decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_real(value, "decn_space_lower")
            check_ndarray_len_eq(value, "decn_space_lower", self.ndecn)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_lower' must be of type numpy.ndarray, an Integral type, or None")
        self._decn_space_lower = value

    # override decn_space setter properties
    @Solution.decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_real(value, "decn_space_upper")
            check_ndarray_len_eq(value, "decn_space_upper", self.ndecn)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_upper' must be of type numpy.ndarray, an Integral type, or None")
        self._decn_space_upper = value





################################## Utilities ###################################
def check_is_RealSolution(v: object, vname: str) -> None:
    """
    Check if object is of type RealSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealSolution):
        raise TypeError("'{0}' must be of type RealSolution.".format(vname))
