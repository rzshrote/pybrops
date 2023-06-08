"""
Module implementing unphased genotyping for dense genotype matrices.
"""

from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import check_is_PhasedGenotypeMatrix

class DenseUnphasedGenotyping(GenotypingProtocol):
    """
    Class implementing unphased genotyping for dense genotype matrices. This
    converts a DensePhasedGenotypeMatrix to a DenseGenotypeMatrix containing
    genotype values.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseUnphasedGenotyping.

        Perform unphased genotyping on a genome.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseUnphasedGenotyping, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def genotype(self, pgmat, miscout = None, **kwargs: dict):
        """
        Genotype a genome. Returned matrix is in {0,1,2} format.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A phased genotype matrix representing the whole simulated genome.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : GenotypeMatrix
            A DenseGenotypeMatrix of genotyped individuals.
        """
        # check for correct data types
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")

        # create genotype matrix
        out = DenseGenotypeMatrix(
            mat = pgmat.mat_asformat("{0,1,2}"),
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp,
            vrnt_chrgrp = pgmat.vrnt_chrgrp,
            vrnt_phypos = pgmat.vrnt_phypos,
            vrnt_name = pgmat.vrnt_name,
            vrnt_genpos = pgmat.vrnt_genpos,
            vrnt_xoprob = pgmat.vrnt_xoprob,
            vrnt_hapgrp = pgmat.vrnt_hapgrp,
            vrnt_hapalt = pgmat.vrnt_hapalt,
            vrnt_hapref = pgmat.vrnt_hapref,
            vrnt_mask = pgmat.vrnt_mask,
            ploidy = pgmat.ploidy,
            **kwargs
        )

        return out
