import numpy

from . import InitializationOperator

from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix

class RandomInitializationOperator(InitializationOperator):
    """docstring for RandomInitializationOperator."""

    def __init__(self, ntaxa, nloci, nchr, ntrait, nburn, rng, **kwargs):
        """
        ntaxa : int
            Number of taxa.
        nloci : int
            Number of loci.
        ntrait : int
            Number of traits.
        nburn : int
            Number of burnin generations.
        """
        super(RandomInitializationOperator, self).__init__(**kwargs)
        self.ntaxa = ntaxa
        self.nloci = nloci
        self.nchr = nchr
        self.ntrait = ntrait
        self.nburn = nburn
        self.rng = rng

    def initialize(self, **kwargs):
        """
        Initialize a breeding program.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        out : tuple
            A tuple containing three dict objects: (geno, bval, gmod)
            Elements of this tuple are as follows:
            Element | Description
            --------+-----------------------------------
            geno    | A dict of genotypic data.
            bval    | A dict of breeding value data.
            gmod    | A dict of genomic models.

            Dictionaries must have the following fields:
            ============================================
            geno : dict
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
            bval : dict
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
                queue      | List of BreedingValueMatrix | Breeding values for populations on queue
                queue_true | List of BreedingValueMatrix | True breeding values for populations on queue
            gmod : dict
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                queue | List of GenomicModel | Genomic models for populations on queue
                true  | GenomicModel         | True genomic model for trait(s)
        """
        # create empty containters
        geno = {
            "cand" : None,
            "main" : None,
            "queue" : [],
        }

        bval = {
            "cand" : None,
            "cand_true" : None,
            "main" : None,
            "main_true" : None,
            "queue" : [],
            "queue_true" : [],
        }

        gmod = {
            "cand" : None,
            "main" : None,
            "queue" : [],
            "true" : None,
        }

        # create the true genomic model
        tgmod = GenericLinearGenomicModel(
            mu = self.rng.uniform(0,100,(self.ntrait,1)),
            beta = self.rng.normal(0,1,(self.nloci,self.ntrait)),
        )

        # create starting genotypes

        raise NotImplementedError("method is abstract")
