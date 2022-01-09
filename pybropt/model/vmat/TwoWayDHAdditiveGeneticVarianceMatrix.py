from pybropt.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix
from pybropt.core.util import srange
from pybropt.model.vmat.util import cov_D1s
from pybropt.core.error import check_is_ndarray

# non-sparse matrix
class TwoWayDHAdditiveGeneticVarianceMatrix(AdditiveGeneticVarianceMatrix):
    """docstring for TwoWayDHAdditiveGeneticVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, var_A, **kwargs):
        super(TwoWayDHAdditiveGeneticVarianceMatrix, self).__init__(
            var_A,
            **kwargs
        )
        self.var_A = var_A

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Variance Matrix Data #################
    def vmat():
        doc = "The vmat property."
        def fget(self):
            return self._var_A
        def fset(self, value):
            self._var_A = value
        def fdel(self):
            del self._var_A
        return locals()
    vmat = property(**vmat())

    def vmat_G():
        doc = "The vmat_G property."
        def fget(self):
            return self._var_A
        def fset(self, value):
            self._var_A = value
        def fdel(self):
            del self._var_A
        return locals()
    vmat_G = property(**vmat_G())

    def var_A():
        doc = "The var_A property."
        def fget(self):
            return self._var_A
        def fset(self, value):
            check_is_ndarray(value, "var_A")
            self._var_A = value
        def fdel(self):
            del self._var_A
        return locals()
    var_A = property(**var_A())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @staticmethod
    def from_pgvmat(pgvmat, lgmod, gmapfn, s, mem = None, gmap = None):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        2-way cross between *inbred* individuals.
        Calculations are derived from Osthushenrich et al. (2017).

        Parameters
        ----------
        pgvmat : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix object containing genetic data.
            It is assumed that the genetic map positions (vrnt_genpos) have
            been interpolated in this object. It is also assumed that this
            object has been grouped into chromosome groups.
        lgmod : LinearGenomicModel
            A LinearGenomicModel object.
        gmapfn : GeneticMapFunction
            A GeneticMapFunction to calculate recombination probabilities.
        s : int
            Selfing generation number to derive gametes from.
            Example:
                s = 0    ->  Derive gametes from F1
                s = 1    ->  Derive gametes from F2
                s = 2    ->  Derive gametes from F3
                ...
                k = inf  ->  Derive gametes from SSD
        mem : int, default = None
            Max square matrix dimensions to compute. If not None, will compute
            covariance matrix segments in chunks of shape (mem, mem) or less.
            If None, full matrix chunks will be computed, which may exceed
            physical memory size! Use this option to save random access memory!
        gmap : GeneticMap, default = None
            A GeneticMap object containing genetic map information. If a
            GeneticMap is provided, genetic map positions will be interpolated
            from the provided GeneticMap instead of using native genetic map
            positions in 'pgvmat'.

        Returns
        -------
        agvmat : TwoWayDHAdditiveGeneticVarianceMatrix
            A TwoWayDHAdditiveGeneticVarianceMatrix.
        """
        # check for grouping
        if not pgvmat.is_grouped():
            raise RuntimeError("pgvmat must be grouped")

        # gather data pointers from pgvmat
        geno = pgvmat.geno                      # (m, n, p) matrix
        chrgrp_stix = pgvmat.vrnt_chrgrp_stix   # (g, ) matrix
        chrgrp_spix = pgvmat.vrnt_chrgrp_spix   # (g, ) matrix
        genpos = pgvmat.vrnt_genpos             # (p, ) matrix
        if gmap is not None:                    # if alternative gmap provided, interpolate genpos.
            genpos = gmap.interp_genpos(pgvmat.vrnt_chrgrp, pgvmat.vrnt_phypos)

        # gather data pointers from lgmod
        beta = lgmod.beta                       # (p, t)

        # gather shapes of data input
        ntrait = beta.shape[1]      # number of traits (t)
        ntaxa = geno.shape[1]       # number of individuals (n)

        # allocate a square matrix for each pairwise variance
        var_A = numpy.zeros((ntrait, ntaxa, ntaxa), dtype='float64')

        # for each linkage group
        for lst, lsp in zip(chrgrp_stix, chrgrp_spix):
            # determine memory chunk step (for iterators)
            step = (lsp - lst) if mem is None else mem
            # for each computational chunk: rb = row block, cb = column block
            for rst,rsp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                for cst,csp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                    # create sparse meshgrid indicating where genetic positions are
                    gi, gj = numpy.meshgrid(genpos[rst:rsp], genpos[cst:csp], indexing='ij', sparse=True)

                    # calculate recombination probability matrix for chunk
                    r = gmapfn.mapfn(numpy.abs(gi - gj)) # (rb,cb)

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1_mat = cov_D1s(r, s) # (rb,cb)

                    # get marker coefficients for rows and columns
                    rbeta = beta[rst:rsp].T # (rb,t)' -> (t,rb)
                    cbeta = beta[cst:csp].T # (cb,t)' -> (t,cb)

                    # for each mate pair (excluding selfs)
                    for female in range(1,ntaxa): # varA row index
                        for male in range(0,female): # varA col index
                            # calculate genotype differences for row, col {-1,0,1}
                            rdgeno = geno[0,female,rst:rsp] - geno[0,male,rst:rsp] # (rb,)
                            cdgeno = geno[0,female,cst:csp] - geno[0,male,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect = rdgeno * rbeta # (rb,)*(t,rb) -> (t,rb)
                            ceffect = cdgeno * cbeta # (cb,)*(t,cb) -> (t,cb)

                            # compute dot product for each trait to get partial variance
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb)[1] -> (t,)
                            varA_part = (reffect @ D1_mat * ceffect).sum(1)

                            # add this partial variance to the lower triangle
                            var_A[:,female,male] += varA_part

        # since varA matrix is symmetrical, copy lower triangle to the upper
        for female in range(1, ntaxa):
            for male in range(0, female):
                var_A[:,male,female] = var_A[:,female,male]

        return TwoWayDHAdditiveGeneticVarianceMatrix(
            var_A = var_A
        )
