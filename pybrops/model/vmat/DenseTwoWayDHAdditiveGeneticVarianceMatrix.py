import numpy

from pybrops.core.util.subroutines import srange
from pybrops.model.vmat.util import cov_D1s
from pybrops.model.vmat.DenseAdditiveGeneticVarianceMatrix import DenseAdditiveGeneticVarianceMatrix

# TODO: implement me
class DenseTwoWayDHAdditiveGeneticVarianceMatrix(DenseAdditiveGeneticVarianceMatrix):
    """docstring for DenseTwoWayDHAdditiveGeneticVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, **kwargs):
        """
        Constructor for the concrete class DenseTwoWayDHAdditiveGeneticVarianceMatrix.

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
        super(DenseTwoWayDHAdditiveGeneticVarianceMatrix, self).__init__(
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
            return (0,1) # (female, male)
        def fset(self, value):
            """Set axis indices for axes that are square"""
            error_readonly("square_axes")
        def fdel(self):
            """Delete axis indices for axes that are square"""
            error_readonly("square_axes")
        return locals()
    square_axes = property(**square_axes())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    # TODO: implement me
    # from_gmod

    @classmethod
    def from_algmod(cls, algmod, pgmat, ncross, nprogeny, s, gmapfn, mem = 1024):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        2-way cross between *inbred* individuals.
        Calculations are derived from Osthushenrich et al. (2017).

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
        s : int
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-------------+-------------------------+
            | Example     | Description             |
            +=============+=========================+
            | ``s = 0``   | Derive gametes from F1  |
            +-------------+-------------------------+
            | ``s = 1``   | Derive gametes from F2  |
            +-------------+-------------------------+
            | ``s = 2``   | Derive gametes from F3  |
            +-------------+-------------------------+
            | ``...``     | etc.                    |
            +-------------+-------------------------+
            | ``s = inf`` | Derive gametes from SSD |
            +-------------+-------------------------+
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
        # check for chromosome grouping
        if not pgmat.is_grouped_vrnt():
            raise RuntimeError("pgmat must be grouped along the vrnt axis")

        # gather shapes of data input
        ntrait = algmod.ntrait                  # number of traits (t)
        ntaxa = pgmat.ntaxa                     # number of individuals (n)

        # gather data pointers from genotypes and linear model
        geno = pgmat.mat                        # (m,n,p) genotype matrix
        chrgrp_stix = pgmat.vrnt_chrgrp_stix    # (g,) chromosome group start indices
        chrgrp_spix = pgmat.vrnt_chrgrp_spix    # (g,) chromosome group stop indices
        genpos = pgmat.vrnt_genpos              # (p,) marker genetic positions
        u = algmod.u                            # (p,t) marker effect coefficients

        # allocate a square matrix for each pairwise variance
        var_A = numpy.zeros(
            (ntaxa, ntaxa, ntrait),             # (n,n,t) variance matrix
            dtype = 'float64'
        )

        # for each linkage group
        for lst, lsp in zip(chrgrp_stix, chrgrp_spix):
            # determine memory chunk step (for iterators)
            step = (lsp - lst) if mem is None else mem
            # for each computational chunk: rb = row block, cb = column block
            for rst,rsp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                for cst,csp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                    # create sparse meshgrid indicating where genetic positions are
                    gi, gj = numpy.meshgrid(
                        genpos[rst:rsp],        # row block genetic positions
                        genpos[cst:csp],        # column block genetic positions
                        indexing='ij',          # use ij indexing
                        sparse=True             # generate a spare array tuple for speed
                    )

                    # calculate recombination probability matrix for chunk
                    r = gmapfn.mapfn(numpy.abs(gi - gj)) # (rb,cb) covariance matrix

                    # calculate a D1 recombination covariance matrix; this is specific to mating scheme
                    D1 = cov_D1s(r, s)  # (rb,cb) covariance matrix

                    # get marker coefficients for rows and columns
                    ru = u[rst:rsp].T # (rb,t)' -> (t,rb)
                    cu = u[cst:csp].T # (cb,t)' -> (t,cb)

                    # for each mate pair (excluding selfs)
                    for female in range(1,ntaxa): # varA row index
                        for male in range(0,female): # varA col index
                            # calculate genotype differences for row, col
                            # phased genotype matrix must be coded in {0,1}
                            # resulting difference matrix has values {-1,0,1}
                            rdgeno = geno[0,female,rst:rsp] - geno[0,male,rst:rsp] # (rb,)
                            cdgeno = geno[0,female,cst:csp] - geno[0,male,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect = rdgeno * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect = cdgeno * cu # (cb,)*(t,cb) -> (t,cb)

                            # compute dot product for each trait to get partial variance sum
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb)[1] -> (t,)
                            var_A_partial = (reffect @ D1 * ceffect).sum(1)

                            # add this partial variance to the lower triangle
                            var_A[female,male,:] += var_A_partial

        # since var_A matrix is symmetrical, copy lower triangle to the upper
        for female in range(1, ntaxa):
            for male in range(0, female):
                var_A[male,female,:] = var_A[female,male,:]

        # construct output
        out = cls(
            mat = var_A,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp
        )

        return out