import breed

class CGS(breed.GenomicSelection):
    """docstring for CGS."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    @classmethod
    def __init__(self, population, cross, method = "CGS"):
        super(CGS, self).__init__(population, cross, method)

        # check that we have marker coefficients
        check_is_ParametricGenomicModel(self._population.genomic_model)

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def objfn(self, sel, objcoeff = None):
        """
        Breeding method objective function. Implement this in derived classes.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        objcoeff : numpy.ndarray, None
            An objective coefficients matrix of shape (t,).
            Where:
                't' is the number of objectives.
            These are used to weigh objectives in the weight sum method.
            If None, do not multiply GEBVs by a weight sum vector.

        Returns
        -------
        cgs : numpy.ndarray
            A GEBV matrix of shape (k,) or (k,t) depending on whether 'objcoeff'
            was specified or not, respectively.
        """
        # calculate GEBVs
        cgs = CGS.objfn_mat(
            sel,
            self._population.geno,
            self._population.genomic_model.coeff
        )

        # if we have objective weights, take dot product for weight sum method
        if objcoeff is not None:
            cgs = cgs.dot(objcoeff)

        return cgs

    @classmethod
    def objfn_vec(self, sel, objcoeff = None):
        """
        Breeding method objective function. Implement this in derived classes.

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
        objcoeff : numpy.ndarray, None
            An objective coefficients matrix of shape (t,) or (t,j).
            Where:
                't' is the number of objectives.
                'j' is the number of selection configurations.
            These are used to weigh objectives in the weight sum method.
            If None, do not multiply GEBVs by a weight sum vector.

        Returns
        -------
        cgs : numpy.ndarray
            A GEBV matrix.
            Shape rules are as follows:
                (j,k)   if objcoeff shape is (t,)
                (j,k,t) if objcoeff shape is (t,j)
                (j,k,t) if objcoeff is None
            Where:
                'j' is the number of selection configurations.
                'k' is the number of individuals to select.
                't' is the number of objectives.
        """
        # calculate GEBVs
        cgs = CGS.objfn_vec_mat(
            sel,
            self._population.geno,
            self._population.genomic_model.coeff
        )

        # if we have objective weights, take dot product for weight sum method
        if objcoeff is not None:
            cgs = cgs.dot(objcoeff)

        return cgs

    @classmethod
    def optimize(self, k, objcoeff, algorithm = None):
        """
        k : int
            Number of individuals to select.
        objcoeff : numpy.ndarray, None
            An objective coefficients matrix of shape (t,).
            Where:
                't' is the number of objectives.
            These are used to weigh objectives in the weight sum method.
            If None, do not multiply GEBVs by a weight sum vector.
        algorithm : None, {'quicksort', 'mergesort', 'heapsort', 'stable'}
            Optional specification for algorithm to use for sorting GEBVs.
            Default is 'quicksort'. Only worry about this for an extremely
            large number of individuals.

        Returns
        -------
        sel : numpy.ndarray
            A selection indices matrix.
             of shape (k,)
            Shape rules are as follows:
                (k,)    if objcoeff shape is (t,)
                (t,k)   if objcoeff is None
            Where:
                'k' is the number of individuals to select.
                't' is the number of objectives.
        """
        # calculate GEBVs
        cgs = self.objfn(None, objcoeff)

        sel = None
        if cgs.ndim == 1:
            cgs_argsort = cgs.argsort(axis = 0)     # get sorted indices
            sel = cgs_argsort[-k:]                  # get last k indices
        elif cgs.ndim == 2:
            cgs_argsort = cgs.argsort(axis = 1)     # get sorted indices
            sel = cgs_argsort[:,-k:]                # get last k indices

        return sel

    @classmethod
    def simulate(self, algorithm):
        """
        """
        raise NotImplementedError("This method is not implemented.")

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def objfn_mat(sel, geno, coeff):
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
        geno : numpy.ndarray, None
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
            Remarks:
                Shape of the matrix is most critical. Underlying matrix
                operations will support other numeric data types.
        coeff : numpy.ndarray, None
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.

        Returns
        -------
        cgs : numpy.ndarray
            A trait prediction (GEBV) matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)
            #return geno.sum(0).dot(coeff)

        # CGS calculation explanation
        # Step 1: get sum of alleles for each individual: shape=(len(sel),L)
        # Step 2: take dot product to get CGS values
        cgs = geno[:,sel].sum(0).dot(coeff)

        return cgs

    @staticmethod
    def objfn_vec_mat(sel, geno, coeff):
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

        CGS selects the 'q' individuals with the largest GEBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape (j,k)
            Where:
                'j' is the number of selection configurations.
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        geno : numpy.ndarray
            An array of allele states. The dtype of 'geno' should be 'uint8'. Array
            shape should be (depth, row, column) = (M, N, L) where 'M' represents
            number of chromosome phases, 'N' represents number of individuals, 'L'
            represents number of markers. Array format should be the 'C' format.
        coeff : numpy.ndarray
            An array of coefficients for allele effects. The dtype of 'coeff'
            should be either 'float32' or 'float64'. This array should be single
            dimensional (column,) = (N,) where 'N' represents number of markers.

        Returns
        -------
        cgs : numpy.ndarray
            Returns an array of floating point number representing GEBVs for each
            individual.
        """
        if sel is None:
            return geno.sum(0).dot(coeff)
        return geno[:,sel].sum(0).dot(coeff)
