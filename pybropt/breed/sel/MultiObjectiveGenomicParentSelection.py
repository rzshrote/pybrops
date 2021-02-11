import numpy

from . import ParentSelectionOperator

from pybropt.core.error import check_is_int

class MultiObjectiveGenomicParentSelection(ParentSelectionOperator):
    """docstring for MultiObjectiveGenomicParentSelection."""

    def __init__(self, k_p, traitwt_p, ncross, nprogeny, target, weight, rng, **kwargs):
        """
        k_p : int
            Number of parents to select.
        target : str or numpy.ndarray
            If target is a string, check value and follow these rules:
                Value         | Description
                --------------+-------------------------------------------------
                "positive"    | Select alleles with the most positive effect.
                "negative"    | Select alleles with the most negate effect.
                "stabilizing" | Set target allele frequency to 0.5.
            If target is a numpy.ndarray, use values as is.
        weight : str or numpy.ndarray
            If weight is a string, check value and follow these rules:
                Value       | Description
                ------------+---------------------------------------------------
                "magnitude" | Assign weights using the magnitudes of regression coefficients.
                "equal"     | Assign weights equally.
        """
        super(MultiObjectiveGenomicParentSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(k_p, "k_p")
        check_is_int(b_p, "b_p")

        # variable assignment
        self.k_p = k_p
        self.traitwt_p = traitwt_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.target = target
        self.weight = weight
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def calc_mkrwt(self, weight, beta):
        if isinstance(weight, str):
            if weight == "magnitude":           # return abs(beta)
                return numpy.absolute(beta)
            elif weight == "equal":             # return 1s matrix
                return numpy.full(beta.shape, 1.0, dtype='float64')
            else:
                raise ValueError("string value for 'weight' not recognized")
        elif isinstance(weight, numpy.ndarray):
            return weight
        else:
            raise TypeError("variable 'weight' must be a string or numpy.ndarray")

    def calc_tfreq(self, target, beta):
        if isinstance(target, str):
            if target == "positive":
                return numpy.float64(beta >= 0.0)   # positive alleles are desired
            elif target == "negative":
                return numpy.float64(beta <= 0.0)   # negative alleles are desired
            elif target == "stabilizing":
                return 0.5                          # both alleles desired
                # return numpy.full(coeff.shape, 0.5, dtype = 'float64')
            else:
                raise ValueError("string value for 'target' not recognized")
        elif isinstance(target, numpy.ndarray):
            return target
        else:
            raise TypeError("variable 'target' must be a string or numpy.ndarray")

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

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitobjwt, traitsum, objsum, **kwargs):
        """
        Return a parent selection objective function.
        """
        mat = geno["cand"].mat                      # (m,n,p) get genotype matrix
        beta = gmod["cand"].beta                    # (p,t) get regression coefficients
        mkrwt = self.calc_mkrwt(self.weight, beta)  # (p,t) get marker weights
        tfreq = self.calc_tfreq(self.target, beta)  # (p,t) get target allele frequencies

        def objfn(sel, mat = mat, tfreq = tfreq, mkrwt = mkrwt, traitobjwt = traitobjwt, traitsum = traitsum, objsum = objsum):
            """
            Multi-objective genomic selection objective function.
                The goal is to minimize this function. Lower is better.
                This is a bare bones function. Minimal error checking is done.

            Given a 2D weight vector 'dcoeff', calculate the Euclidian distance from the
            origin according to:
                dist = dot( dcoeff, F(x) )
                Where:
                    F(x) is a vector of objective functions:
                        F(x) = < f_PAU(x), f_PAFD(x) >

            f_PAU(x):

            Given the provided genotype matrix 'geno' and row selections from it 'sel',
            calculate the selection allele freq. From the selection allele frequencies
            and the target allele frequencies, determine if the target frequencies
            cannot be attained after unlimited generations and selection rounds.
            Multiply this vector by a weight coefficients vector 'wcoeff'.

            f_PAFD(x):

            Given a genotype matrix, a target allele frequency vector, and a vector of
            weights, calculate the distance between the selection frequency and the
            target frequency.

            Parameters
            ==========
            sel : numpy.ndarray, None
                A selection indices matrix of shape (k,)
                Where:
                    'k' is the number of individuals to select.
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
            mat : numpy.ndarray, None
                A int8 binary genotype matrix of shape (m, n, p).
                Where:
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'p' is the number of markers.
                Remarks:
                    Shape of the matrix is most critical. Underlying matrix
                    operations will support other numeric data types.
            tfreq : floating, numpy.ndarray
                A target allele frequency matrix of shape (p, t).
                Where:
                    'p' is the number of markers.
                    't' is the number of traits.
                Example:
                    tfreq = numpy.array([0.2, 0.6, 0.7])
            mkrwt : numpy.ndarray
                A marker weight coefficients matrix of shape (p, t).
                Where:
                    'p' is the number of markers.
                    't' is the number of traits.
                Remarks: Values in 'wcoeff' have an assumption:
                    All values must be non-negative.
            traitobjwt : numpy.ndarray, None
                Combined trait weights and objective weights matrix.
                This matrix must be compatible with shape (2,t).
                Specifying:
                    Only trait weights: (1,t)
                    Only objective weights: (2,1)
                    Trait and objective weights: (2,t)
            traitsum : bool
                Sum across traits.
                If True:
                    (2,t) -> (2,)
                Else:
                    (2,t)
            objsum : bool
                Sum across objectives.
                If True:
                    (2,t) -> (t,)
                Else:
                    (2,t)

            Returns
            =======
            mogs : numpy.ndarray
                A MOGS score matrix of shape (2,t) or other.
            """
            # if no selection, select all
            if sel is None:
                sel = slice(None)

            # generate a view of the geno matrix that only contains 'sel' rows.
            # (m,(k,),p) -> (m,k,p)
            sgeno = geno[:,sel,:]

            # calculate reciprocal number of phases
            rphase = 1.0 / (sgeno.shape[0] * sgeno.shape[1])

            # calculate population frequencies; add axis for correct broadcast
            # (m,k,p).sum[0,1] -> (p,)
            # (p,) * scalar -> (p,)
            # (p,None) -> (p,1) # need (p,1) for broadcasting with (p,t) arrays
            pfreq = (sgeno.sum((0,1)) * rphase)[:,None]

            # calculate some inequalities for use multiple times
            pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
            pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

            # calculate allele unavailability
            allele_unavail = numpy.where(
                tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
                pfreq_lteq_0,           # then set True if sel has allele freq == 0
                numpy.where(            # else
                    tfreq > 0.0,        # if 0.0 < target freq < 1.0
                    numpy.logical_or(   # then set True if pop freq is outside (0.0,1.0)
                        pfreq_lteq_0,
                        pfreq_gteq_1
                    ),
                    pfreq_gteq_1        # else set True if pop freq is >= 1.0
                )
            )

            # calculate distance between target and population
            # (p,t)-(p,1) -> (p,t)
            dist = numpy.absolute(tfreq - pfreq)

            # compute f_PAU(x)
            # (p,t) * (p,t) -> (p,t)
            # (p,t).sum[0] -> (t,)
            pau = (wcoeff * allele_unavail).sum(0)

            # compute f_PAFD(x)
            # (p,t) * (p,t) -> (p,t)
            # (p,t).sum[0] -> (t,)
            pafd = (wcoeff * dist).sum(0)

            # stack to make MOGS matrix
            # (2,t)
            mogs = numpy.stack([pau, pafd])

            # weight traits
            if traitobjwt is not None:
                mogs *= traitobjwt

            # sum across any axes as needed.
            if objsum or traitsum:
                axistuple = tuple()
                if objsum:
                    axistuple += (0,)
                if traitsum:
                    axistuple += (1,)
                mogs = mogs.sum(axistuple)

            return mogs

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
