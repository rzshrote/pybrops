"""
Module defining selection configurations where the decision space is subset in nature.
"""

from numbers import Integral
from typing import Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.breed.prot.sel.cfg.SampledSelectionConfigurationMixin import SampledSelectionConfigurationMixin
from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.random.sampling import axis_shuffle, outcross_shuffle, tiled_choice
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class SubsetSelectionConfiguration(SampledSelectionConfigurationMixin,SelectionConfiguration):
    """
    docstring for SubsetSelectionConfiguration.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            pgmat: PhasedGenotypeMatrix,
            xconfig_decn: numpy.ndarray,
            rng: Optional[Union[Generator,RandomState]] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSelectionConfiguration.
        
        Parameters
        ----------
        ncross : Integral
            Number of cross configurations to consider. Example: ``ncross = 10, nparent = 2``
            specifies 10 two-way crosses.
        nparent : Integral
            The number of parents in a given cross configuration. Example: ``ncross = 10, nparent = 2``
            specifies 10 two-way crosses.
        nmating : Integral, numpy.ndarray
            The number of times an individual cross configuration is executed.
            This becomes important in four-way crosses with heterozygous parents, where
            initial F1 hybrids are unique and can affect the dihybrid composition.
        nprogeny : Integral, numpy.ndarray
            The number of progeny to derive from a mating event.
        pgmat : PhasedGenotypeMatrix
            A genome matrix containing parental candidates
        xconfig_decn : numpy.ndarray
            A decision vector of shape ``(ndecn,)`` containing indices corresponding to individuals in ``pgmat``.
        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments!
        # set shape parameters first
        self.ncross = ncross
        self.nparent = nparent
        # mating parameters second
        self.nmating = nmating
        self.nprogeny = nprogeny
        # set genotypes and cross configuration third
        self.pgmat = pgmat
        self.xconfig_decn = xconfig_decn
        self.rng = rng
        # sample cross configuration
        self.sample_xconfig(return_xconfig=False)

    ############################ Object Properties #############################
    @SampledSelectionConfigurationMixin.xconfig_decn.setter
    def xconfig_decn(self, value: numpy.ndarray) -> None:
        """Set decision vector for calculating the cross configuration matrix."""
        check_is_ndarray(value, "xconfig_decn")
        check_ndarray_ndim(value, "xconfig_decn", 1)
        check_ndarray_dtype_is_integer(value, "xconfig_decn")
        self._xconfig_decn = value

    ############################## Object Methods ##############################
    def sample_xconfig(
            self, 
            return_xconfig: bool = False
        ) -> Union[numpy.ndarray,None]:
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
        # create sample
        out = tiled_choice(
            self.xconfig_decn,
            size = (self.ncross, self.nparent),
            replace = False,
            rng = self.rng
        )

        # at least locally ensure outcrossing
        outcross_shuffle(out, rng = self.rng)

        # shuffle within mating configurations just for good measure
        axis_shuffle(out, 0, rng = self.rng)
        
        # set cross configuration
        self.xconfig = out

        # if we are returning xconfig, then return it, else do not return anything
        if return_xconfig:
            return out

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
