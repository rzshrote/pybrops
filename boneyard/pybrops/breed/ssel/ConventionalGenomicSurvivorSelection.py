from . import SurvivorSelectionOperator

import pybrops.core.random
from pybrops.core.error import check_is_int
from pybrops.core.error import cond_check_is_Generator

class ConventionalGenomicSurvivorSelection(SurvivorSelectionOperator):
    """docstring for ConventionalGenomicSurvivorSelection."""

    def __init__(self, k_s, traitwt_s, rng = None, **kwargs: dict):
        super(ConventionalGenomicSurvivorSelection, self).__init__(**kwargs)

        # check data types
        check_is_int(k_s, "k_s")
        # TODO: check traitwt_s
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.k_s = k_s
        self.traitwt_s = traitwt_s
        self.rng = pybrops.core.random if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def sselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitwt = None, **kwargs: dict):
        """
        Select survivors to serve as potential parents for breeding.

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
        trait : numpy.ndarray, None

        Returns
        -------
        out : tuple
            A tuple containing four objects: (geno_new, bval_new, gmod_new, misc)
            geno_new : dict
                A dict containing genotypic data for all breeding populations.
                Must have the following fields:
                    Field | Type                         | Description
                    ------+------------------------------+----------------------
                    cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                    main  | PhasedGenotypeMatrix         | Main breeding population
            bval_new : dict
                A dict containing breeding value data.
                Must have the following fields:
                    Field      | Type                        | Description
                    -----------+-----------------------------+------------------
                    cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                    cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                    main       | BreedingValueMatrix         | Main breeding population breeding values
                    main_true  | BreedingValueMatrix         | Main breeding population true breeding values
            gmod_new : dict
                A dict containing genomic models.
                Must have the following fields:
                    Field | Type                 | Description
                    ------+----------------------+------------------------------
                    cand  | GenomicModel         | Parental candidate breeding population genomic model
                    main  | GenomicModel         | Main breeding population genomic model
                    true  | GenomicModel         | True genomic model for trait(s)
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if k is None:
            k = self.k_s
        if traitwt is None:
            traitwt = self.traitwt_s

        # get objective function
        objfn = self.sobjfn(
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

        # shallow copy dict
        geno_new = dict(geno)
        bval_new = dict(bval)
        gmod_new = dict(gmod)

        # update cand fields
        geno_new["cand"] = geno["main"].select(sel, axis = 1)
        bval_new["cand"] = bval["main"].select(sel, axis = 0)
        bval_new["cand_true"] = bval["main_true"].select(sel, axis = 0)
        gmod_new["cand"] = gmod["main"]

        misc = {}   # empty dictionary

        return geno_new, bval_new, gmod_new, misc

    def sobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs: dict):
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

    def sobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs: dict):
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
