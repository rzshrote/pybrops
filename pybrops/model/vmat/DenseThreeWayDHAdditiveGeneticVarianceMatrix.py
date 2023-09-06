"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genetic variance estimates calculated using three-way DH
formulae.
"""

from numbers import Integral, Real
from typing import Optional, Union
import numpy
import pandas
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral, check_is_Integral_or_None, check_is_Integral_or_inf
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.util.subroutines import srange
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.vmat.util import cov_D1s
from pybrops.model.vmat.util import cov_D2s
from pybrops.model.vmat.DenseAdditiveGeneticVarianceMatrix import DenseAdditiveGeneticVarianceMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction, check_is_GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class DenseThreeWayDHAdditiveGeneticVarianceMatrix(DenseAdditiveGeneticVarianceMatrix):
    """
    A concrete class for dense additive genetic variance matrices calculated
    for three-way DH progenies.

    The purpose of this concrete class is to implement functionality for:
        1) Genetic variance estimation for three-way DH progenies.
        2) I/O for three-way DH progeny variance matrices.
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
        Constructor for the concrete class DenseThreeWayDHAdditiveGeneticVarianceMatrix.

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
        super(DenseThreeWayDHAdditiveGeneticVarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @DenseAdditiveGeneticVarianceMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 4) # (ntaxa,ntaxa,ntaxa,ntrait)
        self._mat = value

    ############## Square Metadata Properties ##############
    @DenseAdditiveGeneticVarianceMatrix.square_axes.getter
    def square_axes(self) -> tuple:
        """Get axis indices for axes that are square"""
        return (0,1,2) # (recurrent, female, male)

    #################### Trait metadata ####################
    @DenseAdditiveGeneticVarianceMatrix.trait_axis.getter
    def trait_axis(self) -> int:
        return 3

    ######## Expected parental genome contributions ########
    @DenseAdditiveGeneticVarianceMatrix.epgc.getter
    def epgc(self) -> tuple:
        """Get a tuple of the expected parental genome contributions."""
        return (0.5, 0.25, 0.25)

    ############################## Object Methods ##############################
    def to_csv(
            self, 
            fname: str
        ) -> None:
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
    # TODO: implement me
    @classmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Union[Integral,Real],
            gmapfn: GeneticMapFunction,
            **kwargs: dict
        ) -> 'DenseThreeWayDHAdditiveGeneticVarianceMatrix':
        """
        Calculate a symmetrical tensor of progeny variances for each possible
        3-way cross between *inbred* individuals.

        Parameters
        ----------
        gmod : GenomicModel
            Genomic Model with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            variance.
        nself : Integral, Real
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-----------------+-------------------------+
            | Example         | Description             |
            +=================+=========================+
            | ``nself = 0``   | Derive gametes from F1  |
            +-----------------+-------------------------+
            | ``nself = 1``   | Derive gametes from F2  |
            +-----------------+-------------------------+
            | ``nself = 2``   | Derive gametes from F3  |
            +-----------------+-------------------------+
            | ``...``         | etc.                    |
            +-----------------+-------------------------+
            | ``nself = inf`` | Derive gametes from SSD |
            +-----------------+-------------------------+
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTwoWayDHAdditiveGeneticVarianceMatrix
            A matrix of additive genetic variance estimations.
        """
        # type checks
        check_is_GenomicModel(gmod, "gmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral_or_inf(nself, "nself")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")

        # if genomic model is an additive linear genomic model, then use specialized routine
        if isinstance(gmod, AdditiveLinearGenomicModel):
            return cls.from_algmod(
                algmod = gmod, 
                pgmat = pgmat, 
                ncross = ncross, 
                nprogeny = nprogeny, 
                nself = nself, 
                gmapfn = gmapfn, 
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
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Union[Integral,Real], 
            gmapfn: GeneticMapFunction, 
            mem: Union[Integral,None] = 1024
        ) -> 'DenseThreeWayDHAdditiveGeneticVarianceMatrix':
        """
        Calculate a symmetrical tensor of progeny variances for each possible
        3-way cross between *inbred* individuals.
        Calculations are derived from Lehermeier et al. (2019).

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            variance.
        nself : Integral, Real
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-----------------+-------------------------+
            | Example         | Description             |
            +=================+=========================+
            | ``nself = 0``   | Derive gametes from F1  |
            +-----------------+-------------------------+
            | ``nself = 1``   | Derive gametes from F2  |
            +-----------------+-------------------------+
            | ``nself = 2``   | Derive gametes from F3  |
            +-----------------+-------------------------+
            | ``...``         | etc.                    |
            +-----------------+-------------------------+
            | ``nself = inf`` | Derive gametes from SSD |
            +-----------------+-------------------------+
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : Integral, default = 1024
            Memory chunk size to use during matrix operations. If ``None``,
            then memory chunk size is not limited.

            WARNING: Setting ``mem = None`` might result in memory allocation
            errors! For reference, ``mem = 1024`` refers to a matrix of size
            1024x1024, which needs about 8.5 MB of storage. Matrices of course
            need a quadratic amount of memory: :math:`O(n^2)`.

        Returns
        -------
        out : GeneticVarianceMatrix
            A matrix of additive genetic variance estimations.
        """
        # type checks
        check_is_AdditiveLinearGenomicModel(algmod, "algmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral_or_inf(nself, "nself")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_Integral_or_None(mem, "mem")
        
        # check for chromosome grouping
        if not pgmat.is_grouped_vrnt():
            raise ValueError("pgmat must be grouped along the vrnt axis")

        # check for genetic positions
        if pgmat.vrnt_genpos is None:
            raise ValueError("pgmat must have genetic positions")

        # gather shapes of data input
        ntrait = algmod.ntrait                  # number ot traits (t)
        ntaxa = pgmat.ntaxa                     # number of individuals (n)

        # gather data pointers
        geno = pgmat.mat                        # (m,n,p) genotype matrix
        chrgrp_stix = pgmat.vrnt_chrgrp_stix    # (g,) chromosome group start indices
        chrgrp_spix = pgmat.vrnt_chrgrp_spix    # (g,) chromosome group stop indices
        genpos = pgmat.vrnt_genpos              # (p,) marker genetic positions
        u = algmod.u_a                          # (p,t) marker effect coefficients

        # allocate a cube matrix for each 3-way variance
        varA_mat = numpy.zeros(
            (ntaxa, ntaxa, ntaxa, ntrait),      # (n,n,n,t) variance matrix
            dtype='float64'
        )

        # for each linkage group
        for lst, lsp in zip(chrgrp_stix, chrgrp_spix):
            # determine memory chunk step (for iterators)
            step = (lsp - lst) if mem is None else mem
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                for cst,csp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                    # create sparse meshgrid indicating where genetic positions are
                    gi, gj = numpy.meshgrid(
                        genpos[rst:rsp],        # row block genetic positions
                        genpos[cst:csp],        # column block genetic positions
                        indexing = 'ij',        # use ij indexing
                        sparse = True           # generate a spare array tuple for speed
                    )

                    # calculate recombination probability matrix for chunk
                    r = gmapfn.mapfn(numpy.abs(gi - gj)) # (rb,cb) covariance matrix

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1 = cov_D1s(r, nself)  # (rb,cb) covariance matrix

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2 = cov_D2s(r, nself)  # (rb,cb) covariance matrix

                    # get marker coefficients for rows and columns
                    ru = u[rst:rsp].T # (rb,t)' -> (t,rb)
                    cu = u[cst:csp].T # (cb,t)' -> (t,cb)

                    # for each 3-way cross (excluding selfs)
                    # subscript codes:
                    #   1 = recurrent parent
                    #   2 = female
                    #   3 = male
                    for recurr in range(0,ntaxa):           # varA slice index
                        for female in range(0,ntaxa):       # varA row index
                            # calculate genotype differences for row, col
                            # phased genotype matrix must be coded in {0,1}
                            # resulting difference matrix has values {-1,0,1}
                            rdgeno21 = geno[0,female,rst:rsp] - geno[0,recurr,rst:rsp] # (rb,)
                            cdgeno21 = geno[0,female,cst:csp] - geno[0,recurr,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect21 = rdgeno21 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect21 = cdgeno21 * cu # (cb,)*(t,cb) -> (t,cb)

                            # calculate the dot product for each trait to get a
                            # parital variance sum for the female-recurrent
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb).sum(1) -> (t,)
                            varA_part21 = (reffect21 @ D1 * ceffect21).sum(1)

                            # only do lower triangle since symmetrical within each slice
                            for male in range(0,female):    # varA col index
                                # calculate genotype differences for row, col
                                rdgeno23 = geno[0,female,rst:rsp] - geno[0,male,rst:rsp] # (rb,)
                                cdgeno23 = geno[0,female,cst:csp] - geno[0,male,cst:csp] # (cb,)

                                rdgeno31 = geno[0,male,rst:rsp] - geno[0,recurr,rst:rsp] # (rb,)
                                cdgeno31 = geno[0,male,cst:csp] - geno[0,recurr,cst:csp] # (cb,)

                                # calculate effect differences
                                reffect23 = rdgeno23 * ru # (rb,)*(t,rb) -> (t,rb)
                                ceffect23 = cdgeno23 * cu # (cb,)*(t,cb) -> (t,cb)

                                reffect31 = rdgeno31 * ru # (rb,)*(t,rb) -> (t,rb)
                                ceffect31 = cdgeno31 * cu # (cb,)*(t,cb) -> (t,cb)

                                # calculate varA parts for crosses with male
                                # compute dot producte for each trait to get
                                # partial variance sums
                                # (t,rb)x(rb,cb) -> (t,cb)
                                # (t,cb)*(t,cb) -> (t,cb)
                                # (t,cb).sum(1) -> (t,)
                                varA_part23 = (reffect23 @ D2 * ceffect23).sum(1)
                                varA_part31 = (reffect31 @ D1 * ceffect31).sum(1)

                                # calculate varA part for this matrix chunk
                                # (scalar * ((t,) + (t,))) + (t,) -> (t,)
                                varA_part = (2.0 * (varA_part21 + varA_part31)) + varA_part23

                                # add this partial variance to the lower triangle
                                # (1,1,1,t) + (t,) -> (1,1,1,t)
                                varA_mat[recurr,female,male,:] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        # multiplication is faster computationally
        varA_mat *= 0.25

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female in range(1, ntaxa):
            for male in range(0, female):
                varA_mat[:,male,female,:] = varA_mat[:,female,male,:]

        # construct output
        out = cls(
            mat = varA_mat,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp
        )

        return out
