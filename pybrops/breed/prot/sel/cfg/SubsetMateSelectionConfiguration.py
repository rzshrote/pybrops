"""
Module containing representations of subset selection configurations 
where the subset originates from a cross map.
"""

from numbers import Integral
from typing import Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.cfg.MateSelectionConfiguration import MateSelectionConfiguration
from pybrops.breed.prot.sel.cfg.SubsetSelectionConfiguration import SubsetSelectionConfiguration
from pybrops.core.random.sampling import tiled_choice
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class SubsetMateSelectionConfiguration(SubsetSelectionConfiguration,MateSelectionConfiguration):
    """
    Class representing a subset selection configuration where the subset 
    originates from a cross map.
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
            xconfig_xmap: numpy.ndarray,
            rng: Optional[Union[Generator,RandomState]],
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
        self.xconfig_xmap = xconfig_xmap
        self.rng = rng
        # sample cross configuration
        self.sample_xconfig(return_xconfig=False)

    ############################ Object Properties #############################

    ############################## Object Methods ##############################
    def sample_xconfig(
            self, 
            return_xconfig: bool = True
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
        # (ncross,)
        out = tiled_choice(
            self.xconfig_decn,
            size = (self.ncross,),
            replace = False,
            rng = self.rng
        )

        # shuffle within mating configurations just for good measure
        # (ncross,)
        self.rng.shuffle(out)

        # extract mating configurations from cross map
        # (s,d)[(ncross,),:] -> (ncross,d)
        out = self.xconfig_xmap[out,:]

        # set cross configuration
        self.xconfig = out

        # if we are returning xconfig, then return it, else do not return anything
        if return_xconfig:
            return out

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
