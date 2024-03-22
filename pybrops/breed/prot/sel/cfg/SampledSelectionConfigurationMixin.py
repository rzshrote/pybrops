"""
Mixin class to provide functionality for selection configurations that require sampling.
"""

from abc import ABCMeta
from abc import abstractmethod
from typing import Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.core.random.prng import global_prng
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState, check_is_ndarray


class SampledSelectionConfigurationMixin(
        metaclass = ABCMeta,
    ):
    """
    A mixin class to provide functionality for selection configurations which 
    require random samples to be drawn.
    """

    ########################## Special Object Methods ##########################
    # __init__ cannot be defined since this is a mixin class

    ############################ Object Properties #############################
    @property
    def xconfig_decn(self) -> numpy.ndarray:
        """Decision vector for calculating the cross configuration matrix."""
        return self._xconfig_decn
    @xconfig_decn.setter
    def xconfig_decn(self, value: numpy.ndarray) -> None:
        """Set decision vector for calculating the cross configuration matrix."""
        check_is_ndarray(value, "xconfig_decn")
        self._xconfig_decn = value
    
    @property
    def rng(self) -> Union[Generator,RandomState]:
        """A random number source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState,None]) -> None:
        """Set random number source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value

    ############################## Object Methods ##############################
    @abstractmethod
    def sample_xconfig(self, return_xconfig: bool) -> Union[numpy.ndarray,None]:
        """
        Sample a cross configuration from the decision vector and set it as the
        ``xconfig`` value.

        Parameters
        ----------
        return_xconfig : bool
            Whether to return the sampled ``xconfig`` matrix.

        Returns
        -------
        out : numpy.ndarray, None
            The sampled ``xconfig`` matrix if ``return_xconfig`` is true,
            otherwise return nothing.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
