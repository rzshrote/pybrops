from . import ParentSelectionOperator

from pybropt.core.error import check_is_int

class TwoWayOptimalHaploidValueParentSelection(ParentSelectionOperator):
    """docstring for TwoWayOptimalHaploidValueParentSelection."""

    def __init__(self, k_p, traitwt_p, b_p, ncross, nprogeny, rng, **kwargs):
        """
        b_p : int
            Number of haplotype blocks.
        """
        super(TwoWayOptimalHaploidValueParentSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(k_p, "k_p")
        check_is_int(b_p, "b_p")

        # variable assignment
        self.k_p = k_p
        self.traitwt_p = traitwt_p
        self.b_p = b_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def calc_nmkr(self, genpos, chrgrp_stix, chrgrp_spix):
        """
        Determine the number of markers to give per chromosome
        """
        # calculate genetic lengths of each chromosome
        genlen = genpos[chrgrp_spix-1] - genpos[chrgrp_stix]

        # calculate ideal number of markers per chromosome
        nmkr_ideal = (self.b_p / genlen.sum()) * genlen

        # calculate number of chromosome markers, assuming at least one per chromosome
        nmkr_int = numpy.ones(nchr, dtype = "int64")    # start with min of one

        for i in range(self.b_p - nchr):    # forces conformance to self.b_p
            diff = nmkr_int - nmkr_ideal    # take actual - ideal
            ix = diff.argmin()              # get index of lowest difference
            nmkr_int[ix] += 1               # increment at lowest index

        return nmkr_int

    def calc_hbin(self, nmkr, genpos, chrgrp_stix, chrgrp_spix):
        # calculate genetic lengths of each chromosome
        genpos_st = genpos[chrgrp_stix]
        genpos_sp = genpos[chrgrp_spix-1]
        genlen = genpos_sp - genpos_st

        # create empty array
        hbin = numpy.empty(len(genpos), dtype = "int64")

        # TODO: use k-means algorithm to solve marker dispersion
        for ixst,ixsp,gest,gesp,l in zip(chrgrp_stix,chrgrp_spix,genpos_st,genpos_sp,genlen):
            # create linspace
            pos_ideal = numpy.linspace(gest, gesp, num=l+1)


        pass

    def calc_hmat(self, geno, gmod):
        geno        = geno["cand"].mat              # get genotypes
        genpos      = geno["cand"].vrnt_genpos      # get genetic positions
        chrgrp_stix = geno["cand"].vrnt_chrgrp_stix # get chromosome start indices
        chrgrp_spix = geno["cand"].vrnt_chrgrp_spix # get chromosome stop indices
        chrgrp_len  = geno["cand"].vrnt_chrgrp_len  # get chromosome marker lengths

        if (vrnt_chrgrp_stix is None) or (vrnt_chrgrp_spix is None):
            raise RuntimeError("markers are not sorted by chromosome position")

        # get number of chromosomes
        nchr = len(vrnt_chrgrp_stix)

        if self.b_p < nchr:
            raise RuntimeError("number of haplotype blocks is less than the number of chromosomes")

        # calculate number of markers to assign to each chromosome
        nmkr = self.calc_nmkr(genpos, chrgrp_stix, chrgrp_spix)

        if numpy.any(nmkr > chrgrp_len):
            raise RuntimeError("number of markers assigned to chromosome greater than number of available markers")


        pass

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

        gebv = objfn(None)                  # get all GEBVs
        if gebv.ndim == 1:                  # if there is one trait objective
            sel = gebv.argsort()[::-1][:k]  # get indices of top k GEBVs
            self.rng.shuffle(sel)           # shuffle indices
        elif gebv.ndim == 2:                # TODO: ND-selection
            raise RuntimeError("non-dominated genomic selection not implemented")

        misc = {
            "gebv" : gebv
        }

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        mat = geno["cand"].mat      # genotype matrix
        mu = gmod["cand"].mu        # trait means
        beta = gmod["cand"].beta    # regression coefficients

        def objfn(sel, mat = mat, mu = mu, beta = beta, traitwt = traitwt):
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
                    'k' is the number of individuals to select.
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
            mat : numpy.ndarray
                A int8 binary genotype matrix of shape (m, n, p).
                Where:
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'p' is the number of markers.
            mu : numpy.ndarray
                A trait mean matrix of shape (t, 1)
                Where:
                    't' is the number of traits.
            beta : numpy.ndarray
                A trait prediction coefficients matrix of shape (p, t).
                Where:
                    'p' is the number of markers.
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
            # if sel is None, slice all individuals
            if sel is None:
                sel = slice(None)

            # CGS calculation explanation
            # Step 1: (m,k,p) -> (k,p)
            # Step 2: (k,p) . (p,t) -> (k,t)
            # Step 3: (k,t) + (1,t) -> (k,t)
            cgs = mat[:,sel,:].sum(0).dot(beta) + mu.T

            # apply objective weights
            # (k,t) . (t,) -> (k,)
            if traitwt is not None:
                cgs = cgs.dot(traitwt)

            return cgs

        return objfn

    def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a vectorized objective function.
        """
        mat = geno["cand"].mat      # genotype matrix
        mu = gmod["cand"].mu        # trait means
        beta = gmod["cand"].beta    # regression coefficients

        def objfn_vec(sel, mat = mat, mu = mu, beta = beta, traitwt = traitwt):
            """
            Score a population of individuals based on Conventional Genomic Selection
            (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
            Genomic Estimated Breeding Values (GEBV) for a population.

            Parameters
            ----------
            sel : numpy.ndarray
                A selection indices matrix of shape (j,k)
                Where:
                    'j' is the number of selection configurations.
                    'k' is the number of individuals to select.
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
            mat : numpy.ndarray
                A int8 binary genotype matrix of shape (m, n, p).
                Where:
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'p' is the number of markers.
            mu : numpy.ndarray
                A trait mean matrix of shape (t, 1)
                Where:
                    't' is the number of traits.
            beta : numpy.ndarray
                A trait prediction coefficients matrix of shape (p, t).
                Where:
                    'p' is the number of markers.
                    't' is the number of traits.
            traitwt : numpy.ndarray, None
                A trait objective coefficients matrix of shape (t,).
                Where:
                    't' is the number of objectives.
                These are used to weigh objectives in the weight sum method.
                If None, do not multiply GEBVs by a weight sum vector.

            Returns
            -------
            cgs : numpy.ndarray
                A trait GEBV matrix of shape (j,k,t) if objwt is None.
                A trait GEBV matrix of shape (j,k) if objwt shape is (t,)
                OR
                A weighted GEBV matrix of shape (t,).
                Where:
                    'k' is the number of individuals selected.
                    't' is the number of traits.
            """
            # (m,n,p)[:,(j,k),:] -> (m,j,k,p)
            # (m,j,k,p) -> (j,k,p)
            # (j,k,p) . (p,t) -> (j,k,t)
            # (j,k,t) + (1,t) -> (j,k,t)
            cgs = mat[:,sel,:].sum(0).dot(beta) + mu.T

            # (j,k,t) . (t,) -> (j,k)
            if traitwt is not None:
                cgs = cgs.dot(traitwt)

            return cgs

        return objfn_vec

numpy.random.seed(42)
genlen = numpy.random.uniform(0,2,5)
mkrlen = 20*genlen/genlen.sum()
out = numpy.ones(5, dtype="int")
for i in range(20 - 5):
    ix = (out - mkrlen).argmin()
    out[ix] += 1
