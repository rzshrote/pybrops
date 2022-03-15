"""
Module implementing random initializatio of a breeding program.
"""

import numpy

from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.core import random as pbo_rng
from pybrops.core.error import cond_check_is_Generator_or_RandomState
from pybrops.core.error import check_is_int
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

class RandomInitializationOperator(InitializationOperator):
    """
    Class implementing random initializatio of a breeding program.
    """

    def __init__(self, ntaxa, nloci, nchr, ntrait, nburn, rng = None, **kwargs):
        """
        Constructor for the concrete class RandomInitializationOperator.

        Parameters
        ----------
        ntaxa : dict
            Number of taxa.
        nloci : int
            Number of loci.
        ntrait : int
            Number of traits.
        nburn : int
            Number of burnin generations.
        """
        super(RandomInitializationOperator, self).__init__(**kwargs)
        # check argument types
        check_is_int(gqlen, "gqlen")
        check_is_int(ntaxa, "ntaxa")
        check_is_int(nloci, "nloci")
        check_is_int(nchr, "nchr")
        check_is_int(ntrait, "ntrait")
        check_is_int(nburn, "nburn")
        cond_check_is_Generator_or_RandomState(rng, "rng")

        # make assignments
        self.gqlen = gqlen
        self.ntaxa = ntaxa
        self.nloci = nloci
        self.nchr = nchr
        self.ntrait = ntrait
        self.nburn = nburn
        self.rng = pbo_rng if rng is None else rng

    def initialize(self, **kwargs):
        """
        Initialize a breeding program.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple of length 5: ``(genome, geno, pheno, bval, gmod)``

            Where:

            - ``genome`` is a ``dict`` of genomes for the breeding program.
            - ``geno`` is a ``dict`` of genotypes for the breeding program.
            - ``pheno`` is a ``dict`` of phenotypes for the breeding program.
            - ``bval`` is a ``dict`` of breeding values for the breeding program.
            - ``gmod`` is a ``dict`` of genomic models for the breeding program.
        """
        # TODO: implement me
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
        tgmod = DenseAdditiveLinearGenomicModel(
            mu = self.rng.uniform(0,100,(self.ntrait,1)),
            beta = self.rng.normal(0,1,(self.nloci,self.ntrait)),
        )

        # create starting genotypes

        raise NotImplementedError("method is abstract")
