"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genic variance estimates calculated using three-way DH
formulae.
"""

from typing import Optional
import numpy
import pandas
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.vmat.DenseAdditiveGenicVarianceMatrix import DenseAdditiveGenicVarianceMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class DenseThreeWayDHAdditiveGenicVarianceMatrix(DenseAdditiveGenicVarianceMatrix):
    """
    A concrete class for dense additive genic variance matrices calculated
    for two-way DH progenies.

    The purpose of this concrete class is to implement functionality for:
        1) Genic variance estimation for two-way DH progenies.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ):
        """
        Constructor for the concrete class DenseThreeWayDHAdditiveGenicVarianceMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseThreeWayDHAdditiveGenicVarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @DenseAdditiveGenicVarianceMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 4) # (ntaxa,ntaxa,ntaxa,ntrait)
        self._mat = value

    ############## Square Metadata Properties ##############
    @DenseAdditiveGenicVarianceMatrix.square_axes.getter
    def square_axes(self) -> tuple:
        """Get axis indices for axes that are square"""
        return (0,1,2) # (recurrent, female, male)

    #################### Trait metadata ####################
    @DenseAdditiveGenicVarianceMatrix.trait_axis.getter
    def trait_axis(self) -> int:
        return 3

    ######## Expected parental genome contributions ########
    @DenseAdditiveGenicVarianceMatrix.epgc.getter
    def epgc(self) -> tuple:
        """Get a tuple of the expected parental genome contributions."""
        return (0.5, 0.25, 0.25)

    ############################## Object Methods ##############################
    def to_csv(
            self, 
            fname: str
        ) -> None:
        """
        Write a dense three-way additive genic variance matrix to a csv.

        Parameters
        ----------
        fname : str
            Filename to which to write.
        """
        # get names for taxa and traits
        taxa = [str(e) for e in range(self.ntaxa)] if self.taxa is None else self.taxa
        trait = [str(e) for e in range(self.mat_shape[3])]

        # make dictionary to store output columns
        out_dict = {
            "Recurrent": [],
            "Female": [],
            "Male": [],
            "Trait": [],
            "Variance": []
        }

        # construct columns element by element
        for recurrix in range(self.mat_shape[0]):
            for femaleix in range(self.mat_shape[1]):
                for maleix in range(self.mat_shape[2]):
                    for traitix in range(self.mat_shape[3]):
                        out_dict["Recurrent"].append(taxa[recurrix])
                        out_dict["Female"].append(taxa[femaleix])
                        out_dict["Male"].append(taxa[maleix])
                        out_dict["Trait"].append(trait[traitix])
                        out_dict["Variance"].append(self[recurrix,femaleix,maleix,traitix])

        # create a pandas DataFrame from the data
        out_df = pandas.DataFrame(out_dict)

        # write DataFrame to file
        out_df.to_csv(fname, index = False)

    ############################## Class Methods ###############################
    @classmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int,
            **kwargs: dict
        ) -> 'DenseThreeWayDHAdditiveGenicVarianceMatrix':
        """
        Estimate genic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genic variances.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genic
            variance.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : GenicVarianceMatrix
            A matrix of genic variance estimations.
        """
        # type checks
        check_is_GenomicModel(gmod, "gmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(nprogeny, "nprogeny")

        # if genomic model is an additive linear genomic model, then use specialized routine
        if isinstance(gmod, AdditiveLinearGenomicModel):
            return cls.from_algmod(
                algmod = gmod, 
                pgmat = pgmat, 
                nprogeny = nprogeny, 
                **kwargs
            )
        # otherwise raise error since non-linear support hasn't been implemented yet
        else:
            raise NotImplementedError("support for non-linear models not implemented yet")

    @classmethod
    def from_algmod(
            cls, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            mem: int
        ) -> 'DenseThreeWayDHAdditiveGenicVarianceMatrix':
        """
        Estimate genic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genic variances.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genic
            variance.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : GenicVarianceMatrix
            A matrix of additive genic variance estimations.
        """
        # type checks
        check_is_AdditiveLinearGenomicModel(algmod, "algmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")

        # gather shapes of data input
        ntrait = algmod.ntrait                  # number of traits (t)
        ntaxa = pgmat.ntaxa                     # number of individuals (n)
        ploidy = pgmat.ploidy                   # ploidy level (scalar)
        epgc = (0.5,0.25,0.25)                  # parental contributions

        # gather allele frequencies within each taxon
        tafreq = pgmat.tafreq()                 # (n,p) allele frequencies within taxon
        u = algmod.u_a                          # (p,t) marker effect coefficients
        
        # calculate individual locus variance coeffients for binomial distributions
        # (p,t) -> (p,t)
        varcoef = (ploidy * u)**2

        # allocate a square matrix for each pairwise variance
        var_a = numpy.empty(
            (ntaxa,ntaxa,ntrait),               # (n,n,t) variance matrix
            dtype = float
        )

        # for each mate pair (including selfs)
        for recurr in range(0,ntaxa):
            for female in range(0,ntaxa):
                for male in range(0,female):
                    # calculate the cross allele frequency
                    # (n,p)[(3,),:] -> (3,p)
                    # (3,) . (3,p) -> (p,)
                    p = numpy.dot(epgc, tafreq[(recurr,female,male),:])

                    # calculate the variance
                    # scalar - (p,1) -> (p,1)
                    # (p,t) * (p,1) * (p,1) -> (p,t)
                    # (p,t).sum(0) -> (t,)
                    v = (varcoef * p[:,None] * (1.0 - p[:,None])).sum(0)

                    # store in matrix and copy to lower since matrix is symmetrical
                    var_a[recurr,female,male,:] = v
                    var_a[recurr,male,female,:] = v

        # construct output
        out = cls(
            mat = var_a,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp,
            trait = algmod.trait
        )

        return out