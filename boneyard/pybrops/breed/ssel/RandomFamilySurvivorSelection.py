import numpy

from . import SurvivorSelectionOperator

import pybrops.core.random

from pybrops.core.error import check_is_int
from pybrops.core.error import cond_check_is_Generator

class RandomFamilySurvivorSelection(SurvivorSelectionOperator):
    """Select random individuals with families"""

    def __init__(self, k_f, rng = None, **kwargs: dict):
        super(RandomFamilySurvivorSelection, self).__init__(**kwargs)

        # check inputs
        check_is_int(k_f, "k_f")
        cond_check_is_Generator(rng, "rng")

        self.k_f = k_f
        self.rng = pybrops.core.random if rng is None else rng

    def sselect(self, t_cur, t_max, geno, bval, gmod, k = None, **kwargs: dict):
        """
        Select survivors to serve as potential parents for breeding.

        Parameters
        ----------
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

        Returns
        -------
        out : tuple
            A tuple containing four objects: (geno_new, bval_new, gmod_new, misc)
            geno_new : dict
                A dict containing genotypic data for all breeding populations.
                Must have the following fields:
                    Field | Type                         | Description
                    ------+------------------------------+----------------------
                    cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                    main  | PhasedGenotypeMatrix         | Main breeding population
            bval_new : dict
                A dict containing breeding value data.
                Must have the following fields:
                    Field      | Type                        | Description
                    -----------+-----------------------------+------------------
                    cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                    cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                    main       | BreedingValueMatrix         | Main breeding population breeding values
                    main_true  | BreedingValueMatrix         | Main breeding population true breeding values
            gmod_new : dict
                A dict containing genomic models.
                Must have the following fields:
                    Field | Type                 | Description
                    ------+----------------------+------------------------------
                    cand  | GenomicModel         | Parental candidate breeding population genomic model
                    main  | GenomicModel         | Main breeding population genomic model
                    true  | GenomicModel         | True genomic model for trait(s)
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if k is None:
            k = self.k_f

        taxa_grp = geno["main"].taxa_grp        # get taxa groups
        sel = []                            # construct empty list
        for grp in numpy.unique(taxa_grp):  # for each family group
            mask = (taxa_grp == grp)        # mask for each family
            ix = numpy.flatnonzero(mask)    # get indices for each family member
            s = min(mask.sum(), k)          # min(# in family, k_f)
            sel.append(                     # add indices to list
                self.rng.choice(            # randomly sample
                    ix,                     # array of available indices
                    size = s,               # number to sample from each family
                    replace = False         # never double sample survivors
                )
            )
        sel = numpy.concatenate(sel)        # concatenate to numpy.ndarray
        sel.sort()                          # sort indices ascending

        # shallow copy dict
        geno_new = dict(geno)
        bval_new = dict(bval)
        gmod_new = dict(gmod)

        # update cand fields
        geno_new["cand"] = geno["main"].select_taxa(sel)
        bval_new["cand"] = bval["main"].select_taxa(sel)
        bval_new["cand_true"] = bval["main_true"].select_taxa(sel)
        gmod_new["cand"] = gmod["main"]

        misc = {}   # empty dictionary

        return geno_new, bval_new, gmod_new, misc


    def sobjfn(self, t_cur, t_max, geno, bval, gmod, **kwargs: dict):
        """
        Return a survivor selection objective function.
        """
        raise RuntimeError("Random Survivor Selection has no objective function")

    def sobjfn_vec(self, t_cur, t_max, geno, bval, gmod, **kwargs: dict):
        """
        Return a vectorized survivor objective function.
        """
        raise RuntimeError("Random Survivor Selection has no objective function")
