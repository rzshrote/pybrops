"""
Module containing
"""

from numbers import Integral
from numbers import Real
from typing import Union
import numpy
from pybrops.breed.prot.sel.soln.MateSelectionSolution import MateSelectionSolution
from pybrops.breed.prot.sel.soln.IntegerSelectionSolution import IntegerSelectionSolution


class IntegerMateSelectionSolution(IntegerSelectionSolution,MateSelectionSolution):
    """
    Class representing mate selection solutions in subset search spaces.
    """
    
    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,object,None],
            decn_space_upper: Union[numpy.ndarray,object,None],
            decn_space_xmap: numpy.ndarray,
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
        Constructor for IntegerSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(IntegerMateSelectionSolution, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            nsoln = nsoln,
            soln_decn = soln_decn,
            soln_obj = soln_obj,
            soln_ineqcv = soln_ineqcv,
            soln_eqcv = soln_eqcv,
            **kwargs
        )
        # order dependent assignments
        self.decn_space_xmap = decn_space_xmap



################################## Utilities ###################################
def check_is_IntegerMateSelectionSolution(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerMateSelectionSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerMateSelectionSolution):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,IntegerMateSelectionSolution.__name__,type(v).__name__))
