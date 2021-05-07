import random
import numpy

from . import ParentSelectionOperator
from . import SurvivorSelectionOperator

class RandomSelectionOperator(ParentSelectionOperator,SurvivorSelectionOperator):
    """docstring for RandomSelectionOperator."""

    def __init__(self, k_p, ncross, nprogeny, k_s, **kwargs):
        """
        k_p : int
            Number of individuals to select.
        """
        super(RandomSelectionOperator, self).__init__(**kwargs)
        self.k_p = k_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.k_s = k_s

    def pselect(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        """
        Select parents individuals for breeding.

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
        bval : dict
            A dict containing breeding value data.
            Must have the following fields:
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
                queue      | List of BreedingValueMatrix | Breeding values for populations on queue
                queue_true | List of BreedingValueMatrix | True breeding values for populations on queue
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                queue | List of GenomicModel | Genomic models for populations on queue
                true  | GenomicModel         | True genomic model for trait(s)
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing five objects: (pgvmat, sel, ncross, nprogeny, misc)
            pgvmat : PhasedGenotypeVariantMatrix
                A PhasedGenotypeVariantMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgvmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
        """
        # calculate maximum index
        maxix = geno["cand"].ntaxa - 1
        # construct selection array
        sel = numpy.array([random.randint(0,maxix) for _ in range(self.k_p)])

        misc = {}

        return geno["cand"], sel, ncross, nprogeny, misc

    def sselect(self, t_cur, t_max, geno, bval, gmod, **kwargs):
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
                queue      | List of BreedingValueMatrix | Breeding values for populations on queue
                queue_true | List of BreedingValueMatrix | True breeding values for populations on queue
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                queue | List of GenomicModel | Genomic models for populations on queue
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
                    queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                    ""
            bval_new : dict
                A dict containing breeding value data.
                Must have the following fields:
                    Field      | Type                        | Description
                    -----------+-----------------------------+------------------
                    cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                    cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                    main       | BreedingValueMatrix         | Main breeding population breeding values
                    main_true  | BreedingValueMatrix         | Main breeding population true breeding values
                    queue      | List of BreedingValueMatrix | Breeding values for populations on queue
                    queue_true | List of BreedingValueMatrix | True breeding values for populations on queue
            gmod_new : dict
                A dict containing genomic models.
                Must have the following fields:
                    Field | Type                 | Description
                    ------+----------------------+------------------------------
                    cand  | GenomicModel         | Parental candidate breeding population genomic model
                    main  | GenomicModel         | Main breeding population genomic model
                    queue | List of GenomicModel | Genomic models for populations on queue
                    true  | GenomicModel         | True genomic model for trait(s)
            misc : dict
                Miscellaneous output (user defined).
        """
        # get number of parents in main breeding pool
        ntaxa = geno["main"].ntaxa

        # sample indices
        sel = numpy.random.choice(geno["main"].ntaxa, self.k_s, replace = False)

        # sort selections
        sel.sort()

        # extract candidates
        sel = numpy.array(random.sample([i for i in range(ntaxa)], self.k_s))

        # copy/make dictionaries
        geno_new = dict(geno)
        bval_new = dict(bval)
        gmod_new = dict(gmod)
        misc = {}

        # assign new candidates
        geno_new["cand"] = geno["main"].select(sel, axis = 1)
        bval_new["cand"] = bval["main"].select(sel, axis = 1)
        gmod_new["cand"] = gmod["main"]

        return geno_new, bval_new, gmod_new, misc
