from . import PhenotypingProtocol
from pybropt.popgen.ptdf import PandasPhenotypeDataFrame

class TruePhenotyping(PhenotypingProtocol):
    """docstring for TruePhenotyping."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(TruePhenotyping, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def phenotype(self, pgmat, gpmod, **kwargs):
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        gpmod : GenomicModel
            Genomic prediction model to use to determine phenotypes.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            DataFrame containing phenotypes.
        """
        # gather true breeding values
        bvmat = gpmod.predict(pgmat)

        
        raise NotImplementedError("method is abstract")
