class SurvivorSelectionOperator:
    """docstring for SurvivorSelectionOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(SurvivorSelectionOperator, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
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
        raise NotImplementedError("method is abstract")
