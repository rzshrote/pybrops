import numpy
from sklearn.cluster import KMeans

from . import ParentSelectionOperator

import pybropt.core.random
from pybropt.core.error import check_is_int
from pybropt.core.error import cond_check_is_Generator

class OptimalPopulationValueParentSelection(ParentSelectionOperator):
    """docstring for OptimalPopulationValueParentSelection."""

    def __init__(self, k_p, traitwt_p, b_p, ncross, nprogeny, algorithm, rng = None, **kwargs):
        """
        k_p : int
            Number of individuals to select (1/2 number of parents).
        b_p : int
            Number of haplotype blocks.
        """
        super(OptimalPopulationValueParentSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(k_p, "k_p")
        # TODO: check traitwt_p
        check_is_int(b_p, "b_p")
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        # TODO: check algorithm
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.k_p = k_p
        self.traitwt_p = traitwt_p
        self.b_p = b_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.algorithm = algorithm
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def calc_nblk(self, genpos, chrgrp_stix, chrgrp_spix):
        """
        Determine the number of markers to give per chromosome
        """
        # calculate genetic lengths of each chromosome
        genlen = genpos[chrgrp_spix-1] - genpos[chrgrp_stix]

        # calculate ideal number of markers per chromosome
        nblk_ideal = (self.b_p / genlen.sum()) * genlen

        # get number of chromosomes
        nchr = len(chrgrp_stix)

        # calculate number of chromosome markers, assuming at least one per chromosome
        nblk_int = numpy.ones(nchr, dtype = "int64")    # start with min of one

        for i in range(self.b_p - nchr):    # forces conformance to self.b_p
            diff = nblk_int - nblk_ideal    # take actual - ideal
            ix = diff.argmin()              # get index of lowest difference
            nblk_int[ix] += 1               # increment at lowest index

        return nblk_int

    def calc_hbin(self, nblk, genpos, chrgrp_stix, chrgrp_spix):
        # initialize
        hbin = numpy.empty(     # create empty array
            len(genpos),        # same length as number of markers
            dtype = "int64"     # force integer
        )

        # for each chromosome
        k = 0                                   # group counter
        for ixst,ixsp,m in zip(chrgrp_stix,chrgrp_spix,nblk):
            if m == 1:                          # if only one group in chromosome
                hbin[ixst:ixsp] = k             # set to k
                k += m                          # increment k
                continue                        # skip to next iteration
            km = KMeans(                        # fit kmeans model
                n_clusters = m                  # m clusters
            ).fit(genpos[ixst:ixsp,None])       # reshape to (n,1)
            hbin[ixst:ixsp] = k + km.labels_    # add k to groups and copy to haplotype bins
            k += m                              # increment k

        # sanitize output
        k = 0                           # group counter
        prev = hbin[0]                  # get previous group label
        hbin[0] = k                     # set element to k
        for i in range(1,len(hbin)):    # for each element
            if hbin[i] != prev:         # if there is a label switch
                k += 1                  # increment k
                prev = hbin[i]          # overwrite the previous group label
            hbin[i] = k                 # set element to k

        return hbin

    def calc_hbin_bounds(self, hbin):
        hstix = [0]                     # starting indices
        hspix = []                      # stopping indices
        prev = hbin[0]                  # get first group label
        for i in range(1,len(hbin)):    # for each group label
            if hbin[i] != prev:         # if the label is different
                hspix.append(i)         # add the stop index for prev
                prev = hbin[i]          # get next label
                hstix.append(i)         # add the start index for next prev
        hspix.append(len(hbin))         # append last stop index
        hstix = numpy.int64(hstix)      # convert to ndarray
        hspix = numpy.int64(hspix)      # convert to ndarray
        hlen = hspix - hstix            # get lengths of each haplotype group
        return hstix, hspix, hlen

    def calc_hmat(self, gmat, mod):
        mat         = gmat.mat              # get genotypes
        genpos      = gmat.vrnt_genpos      # get genetic positions
        chrgrp_stix = gmat.vrnt_chrgrp_stix # get chromosome start indices
        chrgrp_spix = gmat.vrnt_chrgrp_spix # get chromosome stop indices
        chrgrp_len  = gmat.vrnt_chrgrp_len  # get chromosome marker lengths
        beta        = mod.beta              # get regression coefficients

        if (chrgrp_stix is None) or (chrgrp_spix is None):
            raise RuntimeError("markers are not sorted by chromosome position")

        # get number of chromosomes
        nchr = len(chrgrp_stix)

        if self.b_p < nchr:
            raise RuntimeError("number of haplotype blocks is less than the number of chromosomes")

        # calculate number of marker blocks to assign to each chromosome
        nblk = self.calc_nblk(genpos, chrgrp_stix, chrgrp_spix)

        # ensure there are enough markers per chromosome
        if numpy.any(nblk > chrgrp_len):
            raise RuntimeError(
                "number of haplotype blocks assigned to a chromosome greater than number of available markers"
            )

        # calculate haplotype bins
        hbin = self.calc_hbin(nblk, genpos, chrgrp_stix, chrgrp_spix)

        # define shape
        s = (mat.shape[0], mat.shape[1], self.b_p, beta.shape[1]) # (m,n,b,t)

        # allocate haplotype matrix
        hmat = numpy.empty(s, dtype = beta.dtype)   # (m,n,b,t)

        # get boundary indices
        hstix, hspix, hlen = self.calc_hbin_bounds(hbin)

        # OPTIMIZE: perhaps eliminate one loop using dot function
        # fill haplotype matrix
        for i in range(hmat.shape[3]):                              # for each trait
            for j,(st,sp) in enumerate(zip(hstix,hspix)):           # for each haplotype block
                hmat[:,:,j,i] = mat[:,:,st:sp].dot(beta[st:sp,i])   # take dot product and fill

        return hmat

    def pselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitwt = None, **kwargs):
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        geno : dict
            A dict containing genotypic data for all breeding populations.
            Must have the following fields:
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
        bval : dict
            A dict containing breeding value data.
            Must have the following fields:
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                true  | GenomicModel         | True genomic model for trait(s)
        k : int
        traitwt : numpy.ndarray
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing five objects: (pgvmat, sel, ncross, nprogeny, misc)
            pgvmat : PhasedGenotypeVariantMatrix
                A PhasedGenotypeVariantMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgvmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if k is None:
            k = self.k_p
        if traitwt is None:
            traitwt = self.traitwt_p

        # get objective function
        objfn = self.pobjfn(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            traitwt = traitwt
        )

        # optimize solution using algorithm
        soln_dict = self.algorithm.optimize(
            objfn,
            k = k,
            setspace = numpy.arange(geno["cand"].ntaxa),
            objwt = 1.0
        )

        # extract solution
        sel = soln_dict["soln"]

        misc = {}

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        # get haplotype matrix
        mat = self.calc_hmat(geno["cand"], gmod["cand"])    # (t,m,n,h)

        def objfn(sel, mat = mat, traitwt = traitwt):
            """
            Parameters
            ----------
            sel : numpy.ndarray, None
                A selection indices matrix of shape (k,)
                Where:
                    'k' is the number of individuals to select.
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
            mat : numpy.ndarray
                A haplotype effect matrix of shape (m, n, b, t).
                Where:
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'b' is the number of haplotype blocks.
                    't' is the number of traits.
            traitwt : numpy.ndarray, None
                A trait objective coefficients matrix of shape (t,).
                Where:
                    't' is the number of trait objectives.
                These are used to weigh objectives in the weight sum method.
                If None, do not multiply GEBVs by a weight sum vector.

            Returns
            -------
            cgs : numpy.ndarray
                A GEBV matrix of shape (k, t) if objwt is None.
                A GEBV matrix of shape (k,) if objwt shape is (t,)
                Where:
                    'k' is the number of individuals selected.
                    't' is the number of traits.
            """
            # get max haplotype value
            # (m,n,h,t)[:,(k,),:,:] -> (m,k,h,t)
            # (m,k/2,2,h,t).max((0,1)) -> (h,t)
            # (h,t).sum(0) -> (t,)
            opv = mat[:,sel,:,:].max((0,1)).sum(0)

            # apply objective weights
            # (t,) dot (t,) -> scalar
            if traitwt is not None:
                opv = opv.dot(traitwt)

            return opv

        return objfn

    def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        # get haplotype matrix
        mat = self.calc_hmat(geno["cand"], gmod["cand"])    # (t,m,n,h)

        def objfn(sel, mat = mat, traitwt = traitwt):
            """
            Parameters
            ----------
            sel : numpy.ndarray, None
                A selection indices matrix of shape (j,k)
                Where:
                    'j' is the number of configurations to score.
                    'k' is the number of individuals to select.
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
            mat : numpy.ndarray
                A haplotype effect matrix of shape (t, m, n, b).
                Where:
                    't' is the number of traits.
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'b' is the number of haplotype blocks.
            traitwt : numpy.ndarray, None
                A trait objective coefficients matrix of shape (t,).
                Where:
                    't' is the number of trait objectives.
                These are used to weigh objectives in the weight sum method.
                If None, do not multiply GEBVs by a weight sum vector.

            Returns
            -------
            cgs : numpy.ndarray
                A GEBV matrix of shape (j, t) if traitwt is None.
                A GEBV matrix of shape (j,) if traitwt shape is (t,)
                Where:
                    'j' is the number of configurations to score.
                    't' is the number of traits.
            """
            # get max haplotype value
            # (m,n,h,t)[:,(j,k),:,:] -> (m,j,k,h,t)
            # (m,j,k,h,t).max((0,2)) -> (j,h,t)
            # (j,h,t).sum(1) -> (j,t)
            opv = mat[:,sel,:,:].max((0,2)).sum(1)

            # apply objective weights
            # (j,t) dot (t,) -> scalar
            if traitwt is not None:
                opv = opv.dot(traitwt)

            return opv

        return objfn
