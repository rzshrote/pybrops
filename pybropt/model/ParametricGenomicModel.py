# 3rd party
import numpy

# our libraries
from . import GenomicModel
# import pybropt.model
import pybropt.util

class ParametricGenomicModel(GenomicModel):
    """docstring for ParametricGenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, trait, coeff, mkr_name = None, model_name = "Parametric"):
        """
        Constructor for ParametricGenomicModel class.

        Parameters
        ----------
        trait : numpy.ndarray
            A vector of trait names.
            dtype of this vector must be 'object'.
        coeff : numpy.ndarray
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        model_name : str, None
            Name of the prediction method of this model. Optional.
        """
        # call super constructor
        super(ParametricGenomicModel, self).__init__(trait, mkr_name, model_name)

        # error check the input
        pybropt.util.check_is_matrix(coeff, "coeff")
        pybropt.util.check_matrix_dtype_is_numeric(coeff, "coeff")
        if coeff.ndim == 1:         # if coeff is a vector,
            coeff = coeff[:,None]   # add new axis so shape=(len(coeff),1)
        pybropt.util.check_matrix_ndim(coeff, "coeff", 2)
        pybropt.util.check_matrix_axis_len(coeff, "coeff", 1, self.ntrait)

        # set private variable
        self._coeff = coeff

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def coeff():
        doc = "The coeff property."
        def fget(self):
            return self._coeff
        def fset(self, value):
            self._coeff = value
        def fdel(self):
            del self._coeff
        return locals()
    coeff = property(**coeff())

    def nloci():
        doc = "Number of markers this ParametricGenomicModel has."
        def fget(self):
            return self._coeff.shape[0]
        def fset(self, value):
            error_readonly("nloci")
        def fdel(self):
            error_readonly("nloci")
        return locals()
    nloci = property(**nloci())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def predict(self, geno):
        """
        Predict trait values given a genotype matrix.

        Parameters
        ----------
        geno : numpy.ndarray
            A genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
            Matrix assumptions:
                * Binary (but underlying matrix operations permit any type)
                * 'int8' dtype (underlying matrix operations permit any type)

        Returns
        -------
        pred : numpy.ndarray
            A prediction matrix of shape (n, t)
            Where:
                'n' is the number of individuals.
                't' is the number of traits.
        """
        # take the sum across the matrix stacks O(n^2)
        # NOTE: this reduces time complexity of the dot product
        unphased_geno = geno.sum(0)

        # calculate predictions through matrix operations
        pred = unphased_geno.dot(self._coeff)

        # return predictions
        return pred

    def reorder(self, indices):
        """
        Reorder the internals of this model (for matrix operations).

        Parameters
        ----------
        a : numpy.ndarray
            An array of indices specifying reordering of self.coeff.
        """
        self._coeff = self._coeff[indices]

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def from_csv(fpath, sep = ',', header = 0, mkr_name_ix = 0, coeff_ix = None):
        """
        Extract data from a csv. Needs a header

        coeff_ix : int, list
            Columns to extract for coefficients.
        """
        # check args
        if coeff_ix is None:
            raise ValueError("coeff_ix cannot be None")

        # read the file using pandas.read_csv
        df = pandas.read_csv(
            fpath,
            sep = sep,
            header = header
        )

        #extract marker names
        mkr_name = numpy.string_(df.iloc[:,mkr_name_ix].values)

        # extract coefficients
        coeff = df.iloc[:,coeff_ix].values

        # extract trait names
        trait = None
        if df.columns.dtype != 'object':
            trait = numpy.string_([str(e) for e in df.columns])
        else:
            trait = numpy.string(df.columns)

        # make genomic model
        genomic_model = ParametricGenomicModel(trait, coeff, mkr_name)

        return genomic_model
