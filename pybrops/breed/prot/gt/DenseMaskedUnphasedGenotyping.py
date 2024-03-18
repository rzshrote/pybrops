"""
Module implementing unphased genotyping for dense genotype matrices.
"""

from typing import Optional

import numpy
from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol
from pybrops.core.error.error_type_python import check_is_bool
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_PhasedGenotypeMatrix_has_vrnt_mask, check_is_PhasedGenotypeMatrix

class DenseMaskedUnphasedGenotyping(
        GenotypingProtocol
    ):
    """
    Class implementing unphased genotyping for dense genotype matrices. This
    converts a DensePhasedGenotypeMatrix to a DenseGenotypeMatrix containing
    genotype values.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            invert: bool = False,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseMaskedUnphasedGenotyping.

        Perform unphased genotyping on a genome.

        Parameters
        ----------
        invert : bool
            Whether to invert the variant mask when genotyping.
            If ``False``, then variants with mask values == ``True`` are genotyped.
            If ``True``, then variants with mask values == ``False`` are genotyped.
        kwargs : dict
            Additional keyword arguments.
        """
        self.invert = invert

    ############################ Object Properties #############################
    @property
    def invert(self) -> bool:
        """Whether to invert the variant mask when genotyping."""
        return self._invert
    @invert.setter
    def invert(self, value: bool) -> None:
        """Set whether to invert the variant mask when genotyping."""
        check_is_bool(value, "invert")
        self._invert = value

    ############################## Object Methods ##############################
    def genotype(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> GenotypeMatrix:
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

        # get variant mask
        # True means information is extracted, False means information is not extracted
        mask = pgmat.vrnt_mask

        # extract matrix information
        mat         = pgmat.mat_asformat("{0,1,2}")
        taxa        = pgmat.taxa
        taxa_grp    = pgmat.taxa_grp
        vrnt_chrgrp = pgmat.vrnt_chrgrp
        vrnt_phypos = pgmat.vrnt_phypos
        vrnt_name   = pgmat.vrnt_name
        vrnt_genpos = pgmat.vrnt_genpos
        vrnt_xoprob = pgmat.vrnt_xoprob
        vrnt_hapgrp = pgmat.vrnt_hapgrp
        vrnt_hapalt = pgmat.vrnt_hapalt
        vrnt_hapref = pgmat.vrnt_hapref
        vrnt_mask   = pgmat.vrnt_mask
        ploidy      = pgmat.ploidy

        # apply the mask to relevant data 
        if mask is not None:
            # invert the mask if needeed
            if self.invert:
                mask = ~mask
            # apply mask to fields
            if mat is not None:
                mat = mat[:,mask]
            if vrnt_chrgrp is not None:
                vrnt_chrgrp = vrnt_chrgrp[mask]
            if vrnt_phypos is not None:
                vrnt_phypos = vrnt_phypos[mask]
            if vrnt_name is not None:
                vrnt_name = vrnt_name[mask]
            if vrnt_genpos is not None:
                vrnt_genpos = vrnt_genpos[mask]
            if vrnt_xoprob is not None:
                vrnt_xoprob = vrnt_xoprob[mask]
            if vrnt_hapgrp is not None:
                vrnt_hapgrp = vrnt_hapgrp[mask]
            if vrnt_hapalt is not None:
                vrnt_hapalt = vrnt_hapalt[mask]
            if vrnt_hapref is not None:
                vrnt_hapref = vrnt_hapref[mask]
            if vrnt_mask is not None:
                vrnt_mask = vrnt_mask[mask]

        # create genotype matrix
        out = DenseGenotypeMatrix(
            mat         = mat,
            taxa        = taxa,
            taxa_grp    = taxa_grp,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name   = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask   = vrnt_mask,
            ploidy      = ploidy,
            **kwargs
        )

        # get pointers to taxa group metadata
        taxa_grp_name = pgmat.taxa_grp_name
        taxa_grp_stix = pgmat.taxa_grp_stix
        taxa_grp_spix = pgmat.taxa_grp_spix
        taxa_grp_len  = pgmat.taxa_grp_len

        vrnt_chrgrp_name = None
        vrnt_chrgrp_stix = None
        vrnt_chrgrp_spix = None
        vrnt_chrgrp_len  = None

        # calculate new variant metadata
        if pgmat.is_grouped_vrnt():
            # copy arrays
            vrnt_chrgrp_name = numpy.copy(pgmat.taxa_grp_name)
            vrnt_chrgrp_stix = numpy.copy(pgmat.taxa_grp_stix)
            vrnt_chrgrp_spix = numpy.copy(pgmat.taxa_grp_spix)
            vrnt_chrgrp_len  = numpy.copy(pgmat.taxa_grp_len)

            # convert mask to index list
            masknz = numpy.flatnonzero(mask)

            # calculate new group lengths
            for i,(stix,spix) in enumerate(zip(vrnt_chrgrp_stix,vrnt_chrgrp_spix)):
                vrnt_chrgrp_len[i] = numpy.sum((masknz >= stix) & (masknz < spix))
            
            # overwrite start and stop indices
            vrnt_chrgrp_spix[:] = numpy.cumsum(vrnt_chrgrp_len)
            vrnt_chrgrp_stix[0] = 0
            vrnt_chrgrp_stix[1:] = vrnt_chrgrp_spix[:(len(vrnt_chrgrp_spix)-1)]

        # copy metadata
        out.taxa_grp_name    = taxa_grp_name
        out.taxa_grp_stix    = taxa_grp_stix
        out.taxa_grp_spix    = taxa_grp_spix
        out.taxa_grp_len     = taxa_grp_len
        out.vrnt_chrgrp_name = vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = vrnt_chrgrp_spix
        out.vrnt_chrgrp_len  = vrnt_chrgrp_len

        return out
