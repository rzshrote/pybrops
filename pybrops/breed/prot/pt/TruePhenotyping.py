"""
Module implementing phenotyping protocols for extracting true breeding values.
"""

import math
from numbers import Real
from typing import Optional, Union
import numpy
import pandas

from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class TruePhenotyping(PhenotypingProtocol):
    """
    Class implementing phenotyping protocols for extracting true breeding values.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            gpmod: GenomicModel, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class TruePhenotyping.

        Parameters
        ----------
        gpmod : GenomicModel
            Genomic prediction model to use to determine phenotypes.
        kwargs : dict
            Additional keyword arguments
        """
        super(TruePhenotyping, self).__init__(**kwargs)
        self.gpmod = gpmod

    ############################ Object Properties #############################

    ############### Genomic Model Properties ###############
    @property
    def gpmod(self) -> GenomicModel:
        """Genomic prediction model."""
        return self._gpmod
    @gpmod.setter
    def gpmod(self, value: GenomicModel) -> None:
        """Set genomic prediction model"""
        check_is_GenomicModel(value, "gpmod")
        self._gpmod = value

    ################ Stochastic Parameters #################
    @property
    def var_err(self) -> numpy.ndarray:
        """Error variance for each trait."""
        return numpy.repeat(0.0, self.gpmod.ntrait)
    @var_err.setter
    def var_err(self, value: numpy.ndarray) -> None:
        """Set error variance"""
        error_readonly("var_err")

    ############################## Object Methods ##############################
    def phenotype(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        gpmod : GenomicModel, None
            Genomic prediction model to use to determine phenotypes.
            If None, use default genomic prediction model.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : pandas.DataFrame
            A pandas.DataFrame containing phenotypes for individuals.
        """
        # check argument data types
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        if miscout is not None:
            check_is_dict(miscout, "miscout")

        # calculate true breeding values
        gvmat = self.gpmod.gegv(pgmat)

        # construct dataframe labels
        labels_dict = {}
        taxazfill = math.ceil(math.log10(gvmat.ntaxa))+1
        labels_dict["taxa"] = ["Taxon"+str(i+1).zfill(taxazfill) for i in range(gvmat.ntaxa)] if gvmat.taxa is None else gvmat.taxa
        if gvmat.taxa_grp is not None:
            labels_dict["taxa_grp"] = gvmat.taxa_grp
        labels_df = pandas.DataFrame(labels_dict)

        # construct dataframe data values
        mat = gvmat.unscale()   # calculate unscaled breeding values
        traitzfill = math.ceil(math.log10(gvmat.ntrait))+1
        cols = ["Trait"+str(i+1).zfill(traitzfill) for i in range(gvmat.ntrait)] if gvmat.trait is None else gvmat.trait
        values_df = pandas.DataFrame(
            data = mat,
            columns = cols
        )

        # combine data labels and values
        out_df = pandas.concat([labels_df, values_df], axis = 1)

        return out_df

    def set_h2(
            self, 
            h2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise AttributeError("unsupported operation: heritability always set at 1.0")

    def set_H2(
            self, 
            H2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : float, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise AttributeError("unsupported operation: heritability always set at 1.0")
