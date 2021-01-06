class SelectionOperator:
    """docstring for SelectionOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(SelectionOperator, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, t_cur, t_max, pop, bval, gmod, es, a_max, **kwargs):
        """
        Select individuals for breeding.

        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        pop : dict
            A dict containing genotypic data for all breeding populations.
            Must have the following fields:
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
        bval : dict
            A dict containing breeding value data.
            Must have the following fields:
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
                queue      | List of BreedingValueMatrix | Breeding values for populations on queue
                queue_true | List of BreedingValueMatrix | True breeding values for populations on queue
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------
                main  | GenomicModel         | Main breeding population genomic model
                true  | GenomicModel         | True genomic model for trait(s)
                queue | List of GenomicModel | Genomic models for populations on queue
        es : str
            Selection evolutionary strategy.
            The following options are available:
            Option | Description
            -------+------------
            "m"    | μ strategy: selection considering parents in "main" only
            "m+l"  | μ+λ strategy: selection considering parents in "main" and first on "queue" (queue[0])
            "m,l"  | λ strategy: selection considering parents first on "queue" (queue[0]) only
        a_max : None, int
            Maximum age a line can be to contribute to be selected.
            If None, no age limit is applied.

        Returns
        -------

        """
        raise NotImplementedError("method is abstract")
