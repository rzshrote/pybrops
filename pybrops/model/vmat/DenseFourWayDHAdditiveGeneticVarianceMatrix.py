"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genetic variance estimates calculated using four-way DH
formulae.
"""

import numbers
from typing import Any, Optional, Union
import numpy
import pandas
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.util.subroutines import srange
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.vmat.util import cov_D1s
from pybrops.model.vmat.util import cov_D2s
from pybrops.model.vmat.DenseAdditiveGeneticVarianceMatrix import DenseAdditiveGeneticVarianceMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction, check_is_GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class DenseFourWayDHAdditiveGeneticVarianceMatrix(DenseAdditiveGeneticVarianceMatrix):
    """
    A concrete class for dense additive genetic variance matrices calculated
    for four-way DH progenies.

    The purpose of this concrete class is to implement functionality for:
        1) Genetic variance estimation for four-way DH progenies.
        2) I/O for four-way DH progeny variance matrices.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseFourWayDHAdditiveGeneticVarianceMatrix.

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
        super(DenseFourWayDHAdditiveGeneticVarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Square Metadata Properties ##############
    def square_axes():
        doc = "Axis indices for axes that are square"
        def fget(self):
            """Get axis indices for axes that are square"""
            return (0,1,2,3) # (female2, male2, female1, male1)
        def fset(self, value):
            """Set axis indices for axes that are square"""
            error_readonly("square_axes")
        def fdel(self):
            """Delete axis indices for axes that are square"""
            error_readonly("square_axes")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    square_axes = property(**square_axes())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def to_csv(
            self, 
            fname: str
        ) -> None:
        # get names for taxa and traits
        taxa = [str(e) for e in range(self.ntaxa)] if self.taxa is None else self.taxa
        trait = [str(e) for e in range(self.mat_shape[4])]

        # make dictionary to store output columns
        out_dict = {
            "Female2": [],
            "Male2": [],
            "Female1": [],
            "Male1": [],
            "Trait": [],
            "Variance": []
        }

        # construct columns element by element
        for female2ix in range(self.mat_shape[0]):
            for male2ix in range(self.mat_shape[1]):
                for female1ix in range(self.mat_shape[2]):
                    for male1ix in range(self.mat_shape[3]):
                        for traitix in range(self.mat_shape[4]):
                            out_dict["Female2"].append(taxa[female2ix])
                            out_dict["Male2"].append(taxa[male2ix])
                            out_dict["Female1"].append(taxa[female1ix])
                            out_dict["Male1"].append(taxa[male1ix])
                            out_dict["Trait"].append(trait[traitix])
                            out_dict["Variance"].append(self[female2ix,male2ix,female1ix,male1ix,traitix])

        # create a pandas DataFrame from the data
        out_df = pandas.DataFrame(out_dict)

        # write DataFrame to file
        out_df.to_csv(fname, index = False)

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    # TODO: implement me
    @classmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
            nprogeny: int, 
            nself: Union[int,numbers.Number],
            gmapfn: GeneticMapFunction,
            **kwargs: dict
        ) -> 'DenseFourWayDHAdditiveGeneticVarianceMatrix':
        """
        Calculate a symmetrical tensor of progeny variances for each possible
        3-way cross between *inbred* individuals.

        Parameters
        ----------
        gmod : GenomicModel
            Genomic Model with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : int
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            variance.
        nself : int, numbers.Number
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
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        check_is_int(nself, "nself")
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
            ncross: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            mem: int = 1024
        ) -> 'DenseFourWayDHAdditiveGeneticVarianceMatrix':
        """
        Calculate a symmetrical tensor of progeny variances for each possible
        4-way cross between *inbred* individuals.
        Calculations are derived from Lehermeier et al. (2019).

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : int
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            variance.
        nself : int
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
        mem : int, default = 1024
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
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        check_is_int(nself, "nself")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_int(mem, "mem")
        
        # check for chromosome grouping
        if not pgmat.is_grouped_vrnt():
            raise RuntimeError("pgmat must be grouped along the vrnt axis")

        # gather shapes of data input
        ntrait = algmod.ntrait                  # number ot traits (t)
        ntaxa = pgmat.ntaxa                     # number of individuals (n)

        # gather data pointers
        geno = pgmat.mat                        # (m,n,p) genotype matrix
        chrgrp_stix = pgmat.vrnt_chrgrp_stix    # (g,) chromosome group start indices
        chrgrp_spix = pgmat.vrnt_chrgrp_spix    # (g,) chromosome group stop indices
        genpos = pgmat.vrnt_genpos              # (p,) marker genetic positions
        u = algmod.u_a                          # (p,t) marker effect coefficients

        # allocate a square matrix for each pairwise variance
        varA_mat = numpy.zeros(
            (ntaxa, ntaxa, ntaxa, ntaxa, ntrait),   # (n,n,n,n,t) variance matrix
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

                    # for each 4-way cross (excluding selfs)
                    # subscript codes:
                    #   1 = female 2
                    #   2 = male 2
                    #   3 = female 1
                    #   4 = male 1
                    # TODO: make this more computationally efficient.
                    #       operations can be simplified if you mirror sections
                    #       of the 4d matrix between blocks.
                    #       how to do this is difficult to engineer
                    for female2 in range(0,ntaxa):              # varA block index (change me for efficiency?)
                        for male2 in range(0,ntaxa):            # varA slice index (change me for efficiency?)
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,male2,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                            cdgeno21 = geno[0,male2,cst:csp] - geno[0,female2,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect21 = rdgeno21 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect21 = cdgeno21 * cu # (cb,)*(t,cb) -> (t,cb)

                            # calculate the dot product for each trait to get a
                            # partial variance sum for male2-female2
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb).sum(1) -> (t,)
                            varA_part21 = (reffect21 @ D2 * ceffect21).sum(1)

                            for female1 in range(0,ntaxa):      # varA row index
                                # calculate genotype differences for row, col
                                rdgeno31 = geno[0,female1,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                                cdgeno31 = geno[0,female1,cst:csp] - geno[0,female2,cst:csp] # (cb,)
                                rdgeno32 = geno[0,female1,rst:rsp] - geno[0,male2,rst:rsp] # (rb,)
                                cdgeno32 = geno[0,female1,cst:csp] - geno[0,male2,cst:csp] # (cb,)

                                # calculate effect differences
                                reffect31 = rdgeno31 * ru # (rb,)*(t,rb) -> (t,rb)
                                ceffect31 = cdgeno31 * cu # (cb,)*(t,cb) -> (t,cb)
                                reffect32 = rdgeno32 * ru # (rb,)*(t,rb) -> (t,rb)
                                ceffect32 = cdgeno32 * cu # (cb,)*(t,cb) -> (t,cb)

                                # calculate varA part for female1
                                # (t,rb)x(rb,cb) -> (t,cb)
                                # (t,cb)*(t,cb) -> (t,cb)
                                # (t,cb).sum(1) -> (t,)
                                varA_part31 = (reffect31 @ D1 * ceffect31).sum(1)
                                varA_part32 = (reffect32 @ D1 * ceffect31).sum(1)

                                # only do lower triangle since symmetrical within each slice
                                for male1 in range(0,female1):  # varA col index
                                    # calculate genotype differences for row, col
                                    rdgeno41 = geno[0,male1,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                                    cdgeno41 = geno[0,male1,cst:csp] - geno[0,female2,cst:csp] # (cb,)
                                    rdgeno42 = geno[0,male1,rst:rsp] - geno[0,male2,rst:rsp] # (rb,)
                                    cdgeno42 = geno[0,male1,cst:csp] - geno[0,male2,cst:csp] # (cb,)
                                    rdgeno43 = geno[0,male1,rst:rsp] - geno[0,female1,rst:rsp] # (rb,)
                                    cdgeno43 = geno[0,male1,cst:csp] - geno[0,female1,cst:csp] # (cb,)

                                    # calculate effect differences
                                    reffect41 = rdgeno41 * ru # (rb,)*(t,rb) -> (t,rb)
                                    ceffect41 = rdgeno41 * cu # (cb,)*(t,cb) -> (t,cb)
                                    reffect42 = rdgeno42 * ru # (rb,)*(t,rb) -> (t,rb)
                                    ceffect42 = rdgeno42 * cu # (cb,)*(t,cb) -> (t,cb)
                                    reffect43 = rdgeno43 * ru # (rb,)*(t,rb) -> (t,rb)
                                    ceffect43 = rdgeno43 * cu # (cb,)*(t,cb) -> (t,cb)

                                    # calculate varA parts for crosses with male
                                    # (t,rb)x(rb,cb) -> (t,cb)
                                    # (t,cb)*(t,cb) -> (t,cb)
                                    # (t,cb).sum(1) -> (t,)
                                    varA_part41 = (reffect41 @ D1 * ceffect41).sum(1)
                                    varA_part42 = (reffect42 @ D1 * ceffect42).sum(1)
                                    varA_part43 = (reffect43 @ D2 * ceffect43).sum(1)

                                    # calculate varA part for this matrix chunk
                                    # (t,) + (t,) + (t,) + (t,) + (t,) + (t,) -> (t,)
                                    varA_part = varA_part21 + varA_part31 + varA_part32 + varA_part41 + varA_part42 + varA_part43

                                    # add this partial variance to the lower triangle
                                    # (1,1,1,1,t) + (t,) -> (1,1,1,1,t)
                                    varA_mat[female2,male2,female1,male1,:] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        # multiplication is faster computationally
        varA_mat *= 0.25

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female1 in range(1, ntaxa):
            for male1 in range(0, female1):
                varA_mat[:,:,male1,female1,:] = varA_mat[:,:,female1,male1,:]

        # construct output
        out = cls(
            mat = varA_mat,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp
        )

        return out



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseFourWayDHAdditiveGeneticVarianceMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``DenseFourWayDHAdditiveGeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseFourWayDHAdditiveGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a DenseFourWayDHAdditiveGeneticVarianceMatrix".format(vname))
