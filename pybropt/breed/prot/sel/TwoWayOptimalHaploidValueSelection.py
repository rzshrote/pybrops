import numpy
from sklearn.cluster import KMeans

from . import SelectionProtocol

import pybropt.core.random
from pybropt.core.error import check_is_int
from pybropt.core.error import cond_check_is_Generator

class TwoWayOptimalHaploidValueParentSelection(SelectionProtocol):
    """docstring for TwoWayOptimalHaploidValueParentSelection."""

    def __init__(self, nparent, ncross, nprogeny, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = 1.0, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = 1.0, rng = None, **kwargs):
        """
        Constructor for two-way Optimal Haploid Value Selection (OHV).

        Parameters
        ----------
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        objfn_trans : function, callable, None
        objfn_trans_kwargs : dict, None
        objfn_wt : float, numpy.ndarray
        ndset_trans : function, callable, None
        ndset_trans_kwargs : dict, None
        ndset_wt : float
        rng : numpy.Generator
        """
        super(TwoWayOptimalHaploidValueParentSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(nparent, "nparent")
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        cond_check_is_callable(objfn_trans, "objfn_trans")
        cond_check_is_dict(objfn_trans_kwargs, "objfn_trans_kwargs")
        # TODO: check objfn_wt
        cond_check_is_callable(ndset_trans, "ndset_trans")
        cond_check_is_dict(ndset_trans_kwargs, "ndset_trans_kwargs")
        # TODO: check ndset_wt
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = {} if objfn_trans_kwargs is None else objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = {} if ndset_trans_kwargs is None else ndset_trans_kwargs
        self.ndset_wt = ndset_wt
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
            pgvmat : PhasedGenotypeMatrix
                A PhasedGenotypeMatrix of parental candidates.
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

        # construct pairs
        ntaxa = geno["cand"].ntaxa
        a = []
        for i in range(ntaxa):
            for j in range(i+1,ntaxa):
                a.append(i)
                a.append(j)
        a = numpy.int64(a)

        # get objective function
        objfn = self.pobjfn(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            traitwt = traitwt
        )

        gebv = objfn(a)                         # get all GEBVs
        if gebv.ndim == 1:                      # if there is one trait objective
            pairs = gebv.argsort()[::-1][:k]    # get indices of top k GEBVs
            self.rng.shuffle(pairs)             # shuffle indices
            sel = []
            for i in pairs:
                sel.append(a[2*i])      # female index
                sel.append(a[2*i+1])    # male index
            sel = numpy.int64(sel)
        elif gebv.ndim == 2:                    # TODO: ND-selection
            raise RuntimeError("non-dominated genomic selection not implemented")

        misc = {}

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, method = "single"):
        pass

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        # get haplotype matrix
        mat = self.calc_hmat(geno["cand"], gmod["cand"])    # (t,m,n,h)

        def objfn(sel, mat = mat, traitwt = traitwt):
            """
            Score a population of individuals based on Conventional Genomic Selection
            (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
            Genomic Estimated Breeding Values (GEBV) for a population.

            CGS selects the 'q' individuals with the largest GEBVs.

            Parameters
            ----------
            sel : numpy.ndarray, None
                A selection indices matrix of shape (k,)
                Where:
                    'k' is the number of individuals to select. (k/2 pairs)
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
                A GEBV matrix of shape (k, t) if objwt is None.
                A GEBV matrix of shape (k,) if objwt shape is (t,)
                Where:
                    'k' is the number of individuals selected.
                    't' is the number of traits.
            """
            # get female and male selections
            fsel = sel[0::2]
            msel = sel[1::2]

            # construct selection tuple
            s = tuple(zip(fsel,msel))

            # get max haplotype value
            # (m,n,h,t)[:,(k/2,2),:,:] -> (m,k/2,2,h,t)
            # (m,k/2,2,h,t).max((0,2)) -> (k/2,h,t)
            # (k/2,h,t).sum(1) -> (k/2,t)
            ohv = mat[:,s,:,:].max((0,2)).sum(1)

            # apply objective weights
            # (k/2,t) dot (t,) -> (k/2,)
            if traitwt is not None:
                ohv = ohv.dot(traitwt)

            return ohv

        return objfn

    # TODO: implementation of this function
    # def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
    #     """
    #     Return a vectorized objective function.
    #     """
    #     mat = geno["cand"].mat      # genotype matrix
    #     mu = gmod["cand"].mu        # trait means
    #     beta = gmod["cand"].beta    # regression coefficients
    #
    #     def objfn_vec(sel, mat = mat, mu = mu, beta = beta, traitwt = traitwt):
    #         """
    #         Score a population of individuals based on Conventional Genomic Selection
    #         (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
    #         Genomic Estimated Breeding Values (GEBV) for a population.
    #
    #         Parameters
    #         ----------
    #         sel : numpy.ndarray
    #             A selection indices matrix of shape (j,k)
    #             Where:
    #                 'j' is the number of selection configurations.
    #                 'k' is the number of individuals to select.
    #             Each index indicates which individuals to select.
    #             Each index in 'sel' represents a single individual's row.
    #             If 'sel' is None, use all individuals.
    #         mat : numpy.ndarray
    #             A int8 binary genotype matrix of shape (m, n, p).
    #             Where:
    #                 'm' is the number of chromosome phases (2 for diploid, etc.).
    #                 'n' is the number of individuals.
    #                 'p' is the number of markers.
    #         mu : numpy.ndarray
    #             A trait mean matrix of shape (t, 1)
    #             Where:
    #                 't' is the number of traits.
    #         beta : numpy.ndarray
    #             A trait prediction coefficients matrix of shape (p, t).
    #             Where:
    #                 'p' is the number of markers.
    #                 't' is the number of traits.
    #         traitwt : numpy.ndarray, None
    #             A trait objective coefficients matrix of shape (t,).
    #             Where:
    #                 't' is the number of objectives.
    #             These are used to weigh objectives in the weight sum method.
    #             If None, do not multiply GEBVs by a weight sum vector.
    #
    #         Returns
    #         -------
    #         cgs : numpy.ndarray
    #             A trait GEBV matrix of shape (j,k,t) if objwt is None.
    #             A trait GEBV matrix of shape (j,k) if objwt shape is (t,)
    #             OR
    #             A weighted GEBV matrix of shape (t,).
    #             Where:
    #                 'k' is the number of individuals selected.
    #                 't' is the number of traits.
    #         """
    #         # (m,n,p)[:,(j,k),:] -> (m,j,k,p)
    #         # (m,j,k,p) -> (j,k,p)
    #         # (j,k,p) . (p,t) -> (j,k,t)
    #         # (j,k,t) + (1,t) -> (j,k,t)
    #         cgs = mat[:,sel,:].sum(0).dot(beta) + mu.T
    #
    #         # (j,k,t) . (t,) -> (j,k)
    #         if traitwt is not None:
    #             cgs = cgs.dot(traitwt)
    #
    #         return cgs
    #
    #     return objfn_vec

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, trans, kwargs):
        """
        Score individuals based on Optimal Haploid Value Selection for two-way
        crosses.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select. (k/2 pairs)
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
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        ohv : numpy.ndarray
            A GEBV matrix of shape (t,) if 'trans' is None.
            Otherwise, of shape specified by 'trans'.
            Where:
                't' is the number of traits.
        """
        # get female and male selections
        fsel = sel[0::2]
        msel = sel[1::2]

        # construct selection tuple
        s = tuple(zip(fsel,msel))

        # get max haplotype value
        # (m,n,h,t)[:,(k/2,2),:,:] -> (m,k/2,2,h,t)
        # (m,k/2,2,h,t).max((0,2)) -> (k/2,h,t)
        # (k/2,h,t).sum(1) -> (k/2,t)
        ohv = mat[:,s,:,:].max((0,2)).sum(1)

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            cgs = trans(cgs, **kwargs)

        return ohv

    def objfn(sel, mat, traitwt):
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

        CGS selects the 'q' individuals with the largest GEBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select. (k/2 pairs)
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
            A GEBV matrix of shape (k, t) if objwt is None.
            A GEBV matrix of shape (k,) if objwt shape is (t,)
            Where:
                'k' is the number of individuals selected.
                't' is the number of traits.
        """
        # get female and male selections
        fsel = sel[0::2]
        msel = sel[1::2]

        # construct selection tuple
        s = tuple(zip(fsel,msel))

        # get max haplotype value
        # (m,n,h,t)[:,(k/2,2),:,:] -> (m,k/2,2,h,t)
        # (m,k/2,2,h,t).max((0,2)) -> (k/2,h,t)
        # (k/2,h,t).sum(1) -> (k/2,t)
        ohv = mat[:,s,:,:].max((0,2)).sum(1)

        # apply objective weights
        # (k/2,t) dot (t,) -> (k/2,)
        if traitwt is not None:
            ohv = ohv.dot(traitwt)

        return ohv
