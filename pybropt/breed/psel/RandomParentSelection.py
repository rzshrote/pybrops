import numpy
# import math
import types

import pybropt.core.random

from . import ParentSelectionOperator

from pybropt.core.error import check_is_bool
from pybropt.core.error import check_is_int
from pybropt.core.error import cond_check_is_Generator

class RandomParentSelection(ParentSelectionOperator):
    """Perform random parent selection"""

    def __init__(self, nparent, ncross, nprogeny, replace = False, rng = None, **kwargs):
        """
        Construct a RandomParentSelection operator.

        Parameters
        ----------
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        replace : bool
            Whether to sample parents with replacement.
        rng : numpy.random.Generator or None
            A random number generator source. Used for optimization algorithms.
            If 'rng' is None, use pybropt.core.random module (NOT THREAD SAFE!).
        """
        super(RandomParentSelection, self).__init__(**kwargs)

        # check inputs
        check_is_int(nparent, "nparent")
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        check_is_bool(replace, "replace")
        cond_check_is_Generator(rng, "rng")

        # assign variables
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.replace = replace
        self.rng = pybropt.core.random if rng is None else rng

    def pselect(self, t_cur, t_max, geno, bval, gmod, nparent = None, **kwargs):
        """
        Select parents for breeding.

        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        geno : dict
            A dict containing genotypic data for all breeding populations.
            Must have the following fields:
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
        bval : dict
            A dict containing breeding value data.
            Must have the following fields:
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                true  | GenomicModel         | True genomic model for trait(s)
        nparent : int
            Number of parents to select
        **kwargs
            Additional keyword arguments to be passed to either
            algorithm.optimize (method = "single") or self.ppareto (method =
            "pareto"), depending on method used.
        """
        # if any optional parameters are None, set to defaults.
        if nparent is None:
            nparent = self.nparent

        # get number of taxa
        ntaxa = geno["cand"].ntaxa

        # make selection decision
        sel = self.rng.choice(
            geno["cand"].ntaxa,
            size = nparent,
            replace = self.replace
        )

        # empty misc dictionary
        misc = {}

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, trans, trans_kwargs, **kwargs):
        raise RuntimeError("Random parent selection has no objective function.")

    def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        raise RuntimeError("Random parent selection has no objective function.")

    def ppareto(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        raise RuntimeError("Random parent selection has no Pareto frontier since it has no objective function.")
