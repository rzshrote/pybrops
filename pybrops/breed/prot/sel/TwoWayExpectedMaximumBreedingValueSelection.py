"""
UNDER CONSTRUCTION!
Module implementing selection protocols for two-way expected maximum breeding
value selection.
"""

import numpy

from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error import check_is_bool
from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
from pybrops.core.random.prng import global_prng

class TwoWayExpectedMaximumBreedingValueSelection(SelectionProtocol):
    """
    UNDER CONSTRUCTION!
    Class implementing selection protocols for two-way expected maximum breeding
    value selection.

    # TODO: add formulae for methodology.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, k_p, traitwt_p, ncross, nprogeny, nrep, selfing = False, rng = None, **kwargs: dict):
        """
        k_p : int
            Number of crosses to select (1/2 number of parents).
        nrep : int
            Number of simulation replicates
        """
        super(TwoWayExpectedMaximumBreedingValueSelection, self).__init__(**kwargs)

        # variable assignment
        self.k_p = k_p
        self.traitwt_p = traitwt_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.nrep = nrep
        self.selfing = selfing
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def k_p():
        doc = "The k_p property."
        def fget(self):
            """Get value for k_p."""
            return self._k_p
        def fset(self, value):
            """Set value for k_p."""
            check_is_int(value, "k_p")
            self._k_p = value
        def fdel(self):
            """Delete value for k_p."""
            del self._k_p
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    k_p = property(**k_p())

    # TODO: error checks
    def traitwt_p():
        doc = "The traitwt_p property."
        def fget(self):
            """Get value for traitwt_p."""
            return self._traitwt_p
        def fset(self, value):
            """Set value for traitwt_p."""
            self._traitwt_p = value
        def fdel(self):
            """Delete value for traitwt_p."""
            del self._traitwt_p
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    traitwt_p = property(**traitwt_p())

    def ncross():
        doc = "The ncross property."
        def fget(self):
            """Get value for ncross."""
            return self._ncross
        def fset(self, value):
            """Set value for ncross."""
            check_is_int(value, "ncross")
            self._ncross = value
        def fdel(self):
            """Delete value for ncross."""
            del self._ncross
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ncross = property(**ncross())

    def nprogeny():
        doc = "The nprogeny property."
        def fget(self):
            """Get value for nprogeny."""
            return self._nprogeny
        def fset(self, value):
            """Set value for nprogeny."""
            check_is_int(value, "nprogeny")
            self._nprogeny = value
        def fdel(self):
            """Delete value for nprogeny."""
            del self._nprogeny
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    nprogeny = property(**nprogeny())

    def nrep():
        doc = "The nrep property."
        def fget(self):
            """Get value for nrep."""
            return self._nrep
        def fset(self, value):
            """Set value for nrep."""
            check_is_int(value, "nrep")
            self._nrep = value
        def fdel(self):
            """Delete value for nrep."""
            del self._nrep
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    nrep = property(**nrep())

    def selfing():
        doc = "The selfing property."
        def fget(self):
            """Get value for selfing."""
            return self._selfing
        def fset(self, value):
            """Set value for selfing."""
            check_is_bool(value, "selfing")
            self._selfing = value
        def fdel(self):
            """Delete value for selfing."""
            del self._selfing
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    selfing = property(**selfing())

    def rng():
        doc = "The rng property."
        def fget(self):
            """Get value for rng."""
            return self._rng
        def fset(self, value):
            """Set value for rng."""
            if value is not None:
                check_is_Generator_or_RandomState(value, "rng")
            else:
                value = global_prng
            self._rng = value
        def fdel(self):
            """Delete value for rng."""
            del self._rng
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    rng = property(**rng())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def pselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitwt = None, **kwargs: dict):
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
        kwargs : dict
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

        # construct parent pairs
        ntaxa = geno["cand"].ntaxa
        a = None
        if self.selfing:
            a = numpy.int64([[i,j] for i in range(ntaxa) for j in range(i,ntaxa)])
        else:
            a = numpy.int64([[i,j] for i in range(ntaxa) for j in range(i+1,ntaxa)])

        # get objective function
        objfn_vec = self.pobjfn_vec(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            traitwt = traitwt
        )

        gebv = objfn_vec(a)                     # get all GEBVs
        if gebv.ndim == 1:                      # if there is one trait objective
            pairs = gebv.argsort()[::-1][:k]    # get indices of top k GEBVs
            self.rng.shuffle(pairs)             # shuffle indices
            sel = []
            for i in pairs:
                for j in a[i]:
                    sel.append(j)       # append each parent
            sel = numpy.int64(sel)
        elif gebv.ndim == 2:                    # TODO: ND-selection
            raise RuntimeError("non-dominated genomic selection not implemented")

        misc = {}

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs: dict):
        """
        Return a parent selection objective function.
        """
        # create mating operator
        mateop = TwoWayDHCross(rng = self.rng)

        # extract matrices and models
        pgvmat = geno["cand"]   # get phased genotype variant matrix
        model = gmod["cand"]    # get genomic prediction model

        # get mating parameters
        ncross = self.ncross
        nprogeny = self.nprogeny
        nrep = self.nrep

        def objfn(sel, pgvmat = gmat, model = model, mateop = mateop, t_cur = t_cur, t_max = t_max, ncross = ncross, nprogeny = nprogeny, nrep = nrep, traitwt = traitwt):
            """
            Parameters
            ----------
            sel : numpy.ndarray, None
                A selection indices matrix of shape (2,)
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
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
            # variable for tracking the average
            avg = 0

            # run progeny simulations
            for i in range(nrep):
                # create progeny
                progeny = mateop.mate(
                    t_cur = t_cur,
                    t_max = t_max,
                    pgvmat = pgvmat,
                    sel = sel,
                    ncross = ncross,
                    nprogeny = nprogeny
                )

                # predict progeny breeding values
                bvmat = model.predict(progeny)

                # extract phenotype matrix
                # (nprogeny,t)
                embv = bvmat.mat

                # take dot product if we have trait weights
                # (nprogeny,t) dot (t,) -> (nprogeny,)
                if traitwt is not None:
                    embv = embv.dot(traitwt)

                # find max value in the array and add to avg
                # (nprogeny,t).max(0) -> (t,)
                # (nprogeny,).max(0) -> scalar
                avg += embv.max(0)

            # divide by the number of replicates
            avg = avg / nrep

            return avg

        return objfn

    def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs: dict):
        """
        Return a parent selection objective function.
        """
        # create mating operator
        mateop = TwoWayDHCross(rng = self.rng)

        # extract matrices and models
        pgvmat = geno["cand"]   # get phased genotype variant matrix
        model = gmod["cand"]    # get genomic prediction model

        # get mating parameters
        ncross = self.ncross
        nprogeny = self.nprogeny
        nrep = self.nrep

        def objfn_vec(sel, pgvmat = pgvmat, model = model, mateop = mateop, t_cur = t_cur, t_max = t_max, ncross = ncross, nprogeny = nprogeny, nrep = nrep, traitwt = traitwt):
            """
            Parameters
            ----------
            sel : numpy.ndarray, None
                A selection indices matrix of shape (j,2)
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
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
            # create a list for temporarily storing things
            out = []

            # for each cross
            for i in range(len(sel)):
                # variable for tracking the average
                avg = 0

                # run progeny simulations
                for j in range(nrep):
                    # create progeny
                    progeny, misc = mateop.mate(
                        t_cur = t_cur,
                        t_max = t_max,
                        pgvmat = pgvmat,
                        sel = sel[i],
                        ncross = ncross,
                        nprogeny = nprogeny
                    )

                    # predict progeny breeding values
                    bvmat = model.predict(progeny)

                    # extract phenotype matrix
                    # (nprogeny,t)
                    embv = bvmat.mat

                    # take dot product if we have trait weights
                    # (nprogeny,t) dot (t,) -> (nprogeny,)
                    if traitwt is not None:
                        embv = embv.dot(traitwt)

                    # find max value in the array and add to avg
                    # (nprogeny,t).max(0) -> (t,)
                    # (nprogeny,).max(0) -> scalar
                    avg += embv.max(0)

                # divide by the number of replicates
                avg = avg / nrep

                # append the average to the list
                out.append(avg)

            # stack the output
            # j*scalar -> (j,)
            # j*(t,) -> (j,t)
            out = numpy.stack(out)

            return out

        return objfn_vec
