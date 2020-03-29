from . import GenomicSelection

class WGS(GenomicSelection):
    """docstring for CGS."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    @classmethod
    def __init__(self, population, cross, method = "WGS"):
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
        # calculate weighted GEBVs
        wgs = WGS.objfn_mat(
            sel,
            self._population.geno,
            self._population.genomic_model.coeff
        )

        # if we have objective weights, take dot product for weight sum method
        if objcoeff is not None:
            wgs = wgs.dot(objcoeff)

        return wgs

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
        wgs = WGS.objfn_vec_mat(
            sel,
            self._population.geno,
            self._population.genomic_model.coeff
        )

        # if we have objective weights, take dot product for weight sum method
        if objcoeff is not None:
            wgs = wgs.dot(objcoeff)

        return wgs

    @classmethod
    def optimize(self, k, objcoeff = None, algorithm = None):
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
        wgs = self.objfn(None, objcoeff)

        sel = None
        if wgs.ndim == 1:
            wgs_argsort = wgs.argsort(axis = 0)     # get sorted indices
            sel = wgs_argsort[-k:]                  # get last k indices
        elif wgs.ndim == 2:
            wgs_argsort = wgs.argsort(axis = 1)     # get sorted indices
            sel = wgs_argsort[:,-k:]                # get last k indices

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
    def objfn_mat(sel, geno, coeff = None, weight = None, wcoeff = None):
        """
        Score a population of individuals based on Weighted Genomic Selection (WGS)
        (Goddard, 2009; Jannink, 2010). Scoring for WGS is defined as the adjusted
        Genomic Estimated Breeding Values (GEBV) for an individual. Marker
        coefficients are adjusted according to the frequency of the most beneficial
        allele.

        WGS selects the 'q' individuals with the largest weighted GEBVs.

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
        weight : numpy.ndarray, None
            A allele effect weight coefficients matrix of shape (p,t)
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            This is the frequency of the most beneficial allele. If this
            frequency is 0, set it to 1 instead.
            Ignored if 'wcoeff' is provided.
        wcoeff : numpy.ndarray, None
            A allele weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            This is the Hadamard product of 'coeff' and 1/sqrt('weight')
            If None, will be calculated.

        Returns
        -------
        wgs : numpy.ndarray
            A weighted GEBV matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # if no wcoeff matrix, calculate it
        if wcoeff is None:
            # if no weight matrix is given, calculate it
            if weight is None:
                weight = WGS.weight(coeff, geno)
            # calculate wcoeff
            wcoeff = WGS.wcoeff(coeff, weight)

        # CGS calculation explanation
        # Step 1: get sum of alleles for each individual: shape=(len(sel),L)
        # Step 2: take dot product to get CGS values
        wgs = geno[:,sel].sum(0).dot(wcoeff)

        return wgs

    @staticmethod
    def objfn_vec_mat(sel, geno, coeff = None, weight = None, wcoeff = None):
        """
        Score a population of individuals based on Weighted Genomic Selection (WGS)
        (Goddard, 2009; Jannink, 2010). Scoring for WGS is defined as the adjusted
        Genomic Estimated Breeding Values (GEBV) for an individual. Marker
        coefficients are adjusted according to the frequency of the most beneficial
        allele.

        WGS selects the 'q' individuals with the largest weighted GEBVs.

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
        weight : numpy.ndarray, None
            A allele effect weight coefficients matrix of shape (p,t)
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        wcoeff : numpy.ndarray, None
            A allele weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            This is the Hadamard product of 'coeff' and 'wcoeff'

        Returns
        -------
        wgs : numpy.ndarray
            Returns an array of floating point number representing GEBVs for each
            individual.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # if no wcoeff matrix, calculate it
        if wcoeff is None:
            wcoeff = WGS.wcoeff(coeff, weight)

        # CGS calculation explanation
        # Step 1: get sum of alleles for each individual: shape=(len(sel),L)
        # Step 2: take dot product to get CGS values
        wgs = geno[:,sel].sum(0).dot(wcoeff)

        return wgs

    @staticmethod
    def weight(coeff, geno):
        """
        Parameters
        ----------
        coeff : numpy.ndarray, None
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        """
        # get the allele frequencies within the population.
        wght = geno.sum((0,1)) / (geno.shape[0] * geno.shape[1])

        # duplicate matrix to a (p,t) matrix
        wght = numpy.repeat(wght[:,None], coeff.shape[1], 1)

        # get mask of where coeff is deleterious
        mask = coeff < 0

        # invert wght where coeff is negative
        wght[mask] = 1.0 - wght[mask]

        # Where the weight is 0.0, set to 1.0 (we cannot have division by zero)
        wght[wght == 0.0] = 1.0

        return wght

    @staticmethod
    def wcoeff(coeff, weight):
        """
        Parameters
        ----------
        coeff : numpy.ndarray, None
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        weight : numpy.ndarray, None
            A allele effect weight coefficients matrix of shape (p,t)
            Where:
                'p' is the number of markers.
                't' is the number of traits.

        Returns
        -------
        wcoeff : numpy.ndarray
            A allele weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            This is the Hadamard product of 'coeff' and 'wcoeff'
        """
        # divide coeff by sqrt(weight)
        wcoeff = coeff / numpy.sqrt(weight)

        return weight
