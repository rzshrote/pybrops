from . import SelectionProtocol

import pybropt.core.random
from pybropt.core.error import cond_check_is_Generator
from pybropt.core.error import check_is_int

class ConventionalPhenotypicSelection(SelectionProtocol):
    """docstring for ConventionalPhenotypicSelection."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, k_p, traitwt_p, ncross, nprogeny, rng = None, **kwargs):
        super(ConventionalPhenotypicSelection, self).__init__(**kwargs)

        # check data types
        check_is_int(k_p, "k_p")
        # TODO: check traitwt_p
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.k_p = k_p
        self.traitwt_p = traitwt_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
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

        ebv = objfn(None)                   # get all EBVs
        if ebv.ndim == 1:                   # if there is one trait objective
            sel = ebv.argsort()[::-1][:k]   # get indices of top k GEBVs
            self.rng.shuffle(sel)           # shuffle indices
        elif ebv.ndim == 2:                 # TODO: ND-selection
            raise RuntimeError("non-dominated phenotypic selection not implemented")

        misc = {}   # empty dictionary

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        mat = bval["cand"].mat      # breeding value matrix

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
                A breeding value matrix of shape (n, t).
                Where:
                    'n' is the number of individuals.
                    't' is the number of traits.
            traitwt : numpy.ndarray, None
                A trait objective coefficients matrix of shape (t,).
                Where:
                    't' is the number of trait objectives.
                These are used to weigh objectives in the weight sum method.
                If None, do not multiply GEBVs by a weight sum vector.

            Returns
            -------
            cps : numpy.ndarray
                A matrix of shape (k, t) if objwt is None.
                A matrix of shape (k,) if objwt shape is (t,)
                Where:
                    'k' is the number of individuals selected.
                    't' is the number of traits.
            """
            # if sel is None, slice all individuals
            if sel is None:
                sel = slice(None)

            # (n,t) -> (k,t)
            cps = mat[sel,:]

            # apply objective weights
            # (k,t) . (t,) -> (k,)
            if traitwt is not None:
                cps = cps.dot(traitwt)

            return cps

        return objfn

    def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a vectorized objective function.
        """
        mat = bval["cand"].mat      # breeding value matrix

        def objfn_vec(sel, mat = mat, traitwt = traitwt):
            """
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
            # (n,t) -> (k,t)
            # (n,t)[(j,k),:] -> (j,k,t)
            cps = mat[sel,:]

            # (j,k,t) . (t,) -> (j,k)
            if traitwt is not None:
                cps = cps.dot(traitwt)

            return cps

        return objfn_vec

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, trans, kwargs):
        """
        Score a selection configuration based on its breeding values
        (Conventional Phenotype Selection).

        CPS selects the 'q' individuals with the largest EBVs.

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
            A breeding value matrix of shape (n,t).
            Where:
                'n' is the number of individuals.
                't' is the number of traits.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        cps : numpy.ndarray
            A EBV matrix of shape (t,) if 'trans' is None.
            Otherwise, of shape specified by 'trans'.
            Where:
                't' is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # CPS calculation explanation
        # Step 1: (n,t) -> (k,t)        # select individuals
        # Step 2: (k,t).sum(0) -> (t,)  # sum across all individuals
        cps = mat[sel,:].sum(0)

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            cps = trans(cps, **kwargs)

        return cps

    @staticmethod
    def objfn_vec_static(sel, mat, trans, kwargs):
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape (j,k)
            Where:
                'j' is the number of selection configurations.
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, score each individual separately: (n,1)
        mat : numpy.ndarray
            A genotype matrix of shape (n, p).
            Where:
                'n' is the number of individuals.
                'p' is the number of markers.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        cps : numpy.ndarray
            A GEBV matrix of shape (j,t) if 'trans' is None.
            Otherwise, of shape specified by 'trans'.
            Where:
                'j' is the number of selection configurations.
                't' is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            n = mat.shape[0]
            sel = numpy.arange(n).reshape(n,1)

        # CPS calculation explanation
        # Step 1: (n,t)[(j,k),:] -> (j,k,t) # select individuals
        # Step 2: (j,k,t).sum(1) -> (j,t)   # sum across all individuals
        cps = mat[sel,:].sum(1)

        # apply transformations
        # (j,t) ---trans---> (?,?)
        if trans:
            cps = trans(cps, **kwargs)

        return cps
