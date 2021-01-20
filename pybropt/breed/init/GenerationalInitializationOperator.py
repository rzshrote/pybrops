from . import InitializationOperator

from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.breed.eval import NoGxEEvaluationOperator

class GenerationalInitializationOperator(InitializationOperator):
    """docstring for GenerationalInitializationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, seed_geno, seed_bval, seed_gmod, evalop, **kwargs):
        super(GenerationalInitializationOperator, self).__init__(**kwargs)

        self.seed_geno = seed_geno
        self.seed_bval = seed_bval
        self.seed_gmod = seed_gmod

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
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
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @staticmethod
    def from_vcf(fname, size, rng, gmult, evalop, gmod_true, replace = True):
        """
        Create a GenerationalInitializationOperator from a VCF file.

        Initializes a "main" population with genotypes, and queue populations
        of various lengths. Uses the individuals in "main" to select a set of
        parental candidates using a provided SurvivorSelectionOperator. Then,
        a provided ParentSelectionOperator, and Mating operator is used to
        select and mate parents to create one additional generation that is
        added to the queue.

        Parameters
        ----------
        fname : str
            VCF file name.
        size : dict
            Field | Type        | Description
            ------+-------------+-----------------------------------------------
            cand  | None        | None
            main  | int         | Number of taxa in main breeding population.
            queue | list of int | Number of taxa in breeding populations on queue.
        rng : numpy.random.Generator
            A random number generator object.
        replace : bool, default = True
            Whether genotype sampling is with or without replacement.
            If replace == False:
                If the number of genotypes in the provided file is less than
                the sum of the required initial number of genotypes
                (main + sum(queue)), then sample all genotypes and fill
                remaining genotypes with genotypes sampled without replacement.
        initprogeny : int, default = 0
            Number of times progeny should be selected from the main population.
            If
        """
        # perform error checks
        check_is_dict(size, "size")
        check_keys_in_dict(size, "size", "main", "queue")

        # step 1: load genotype matrix
        dpgvmat = DensePhasedGenotypeVariantMatrix.from_vcf(fname)

        # step 2: count individuals to sample and number of available taxa
        nindiv = size["main"] + sum(size["queue"])
        ntaxa = dpgvmat.ntaxa

        # step 3: sample individuals
        if replace:
            ix = rng.choice(ntaxa, nindiv, replace = True)
        else:
            ix = numpy.concatenate([
                numpy.repeat(numpy.arange(ntaxa), nindiv//ntaxa),   # sample all
                rng.choice(ntaxa, nindiv, replace = False)          # sample remaining
            ])
            rng.shuffle(ix) # shuffle everything

        # define seed_geno
        seed_geno = {
            "cand" : None,
            "main" : None,
            "queue" : []
        }

        # step 4: randomly partition into cohorts
        stix = 0            # start pointer
        spix = size["main"] # stop pointer

        # populate fields with genotypes
        tmp = dpgvmat.select(ix[stix:spix], axix = 1)
        tmp.taxa_grp = -(numpy.arange(stix,spix) + gmult)
        seed_geno["main"] = dpgvmat.select(ix[stix:spix], axis = 1)
        seed_geno["main"].taxa_grp = -(numpy.arange(stix, spix) + gmult)
        for l in size["queue"]:         # for each length in queue
            stix = spix                 # advance start pointer
            spix += l                   # advance stop pointer
            seed_geno["queue"].append(  # subset and append values
                dpgvmat.select(
                    ix[stix:spix],
                    axis = 1
                )
            )

        # add taxa groups
        seed_geno["main"] =

        # define seed_bval
        seed_bval = {
            "cand" : None,
            "cand_true" : None,
            "main" : None,
            "main_true" : None,
        }

        # populate fields with breeding values
        seed_bval["main"]

        # step 3: assign partitions into their corresponding groups
        # assign cand




################################################################################
################################## Utilities ###################################
################################################################################
def is_GenerationalInitializationOperator(v):
    return isinstance(v, GenerationalInitializationOperator)

def check_is_GenerationalInitializationOperator(v, varname):
    if not isinstance(v, GenerationalInitializationOperator):
        raise TypeError("'%s' must be a GenerationalInitializationOperator." % varname)

def cond_check_is_GenerationalInitializationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenerationalInitializationOperator(v, varname)
