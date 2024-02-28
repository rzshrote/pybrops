"""
Module containing abstract class definitions for selection configurations
"""

from numbers import Integral
from typing import Union

import numpy
from pybrops.breed.prot.sel.cfg.MateSelectionConfiguration import MateSelectionConfiguration
from pybrops.breed.prot.sel.cfg.SimpleSelectionConfiguration import SimpleSelectionConfiguration
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class SimpleMateSelectionConfiguration(SimpleSelectionConfiguration,MateSelectionConfiguration):
    """
    A simple selection configuration class containing the basic necessities for
    a SelectionConfiguration object.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            pgmat: PhasedGenotypeMatrix,
            xconfig: numpy.ndarray,
            xconfig_xmap: numpy.ndarray,
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
        xconfig : numpy.ndarray
            A mating configuration matrix of shape ``(ncross,nparent)``.
            This matrix contains indices corresponding to parents in ``pgmat``
            and specifies the manner in which individuals are to be mated.
        xconfig_xmap : numpy.ndarray
            A cross map corresponding to the decision space.
        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments!
        super(SimpleMateSelectionConfiguration, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
            pgmat = pgmat,
            xconfig = xconfig,
            **kwargs
        )
        # make assignments to cross map
        self.xconfig_xmap = xconfig_xmap

    ############################ Object Properties #############################
    
    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
