import numpy

from . import SurvivorSelectionOperator

class FamilyPhenotypicSurvivorSelection(SurvivorSelectionOperator):
    """docstring for FamilyPhenotypicSurvivorSelection."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, k_f, traitwt_f, rng, **kwargs):
        super(FamilyPhenotypicSurvivorSelection, self).__init__(**kwargs)
        self.k_f = k_f
        self.traitwt_f = traitwt_f
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def sselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitwt = None, **kwargs):
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
            k = self.k_f
        if traitwt is None:
            traitwt = self.traitwt_f

        # get objective function
        objfn = self.sobjfn(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            traitwt = traitwt
        )

        taxa_grp = bval["main"].taxa_grp        # get taxa groups
        ebv = objfn(None)                       # get all EBVs
        if ebv.ndim == 1:                       # if there is one trait objective
            sel = []                            # construct empty list
            ord = ebv.argsort()[::-1]           # get order of EBVs
            for taxa in numpy.unique(taxa_grp): # for each family
                mask = (taxa_grp[ord] == taxa)  # mask for each family
                s = min(mask.sum(), k)          # min(# in family, k_f)
                sel.append(ord[mask][:s])       # add indices to list
            sel = numpy.concatenate(sel)        # concatenate to numpy.ndarray
            #self.rng.shuffle(sel)               # shuffle indices
        elif ebv.ndim == 2:                     # TODO: ND-selection
            raise RuntimeError("non-dominated phenotypic selection not implemented")

        # shallow copy dict
        geno_new = dict(geno)
        bval_new = dict(bval)
        gmod_new = dict(gmod)

        # update cand fields
        # print(geno["main"].mat.shape)
        # print(bval["main"].mat.shape)
        # print("sel:", sel)
        geno_new["cand"] = geno["main"].select(sel, axis = 1)
        bval_new["cand"] = bval["main"].select(sel, axis = 0)
        bval_new["cand_true"] = bval["main_true"].select(sel, axis = 0)
        gmod_new["cand"] = gmod["main"]

        misc = {}   # empty dictionary

        return geno_new, bval_new, gmod_new, misc

    def sobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        mat = bval["main"].mat      # breeding value matrix

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

    def sobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a vectorized objective function.
        """
        mat = bval["main"].mat      # breeding value matrix

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