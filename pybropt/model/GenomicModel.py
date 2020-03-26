# append paths
import sys
import os
model_dir = os.path.dirname(os.path.realpath(__file__))     # get pybropt/model
pybropt_dir = os.path.dirname(model_dir)                    # get pybropt
sys.path.append(pybropt_dir)                                # append pybropt

# 3rd party

# our libraries
import util


class GenomicModel:
    """docstring for GenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, trait, model_name = None):
        """
        Constructor for the GenomicModel class.

        Parameters
        ----------
        trait : numpy.ndarray
            A vector of strings. Must be of type numpy.string_
        model_name : str
            Name of the model type.
        """
        # error check the input
        util.check_is_matrix(trait, "trait")
        util.check_matrix_ndim(trait, "trait", 1)
        util.check_matrix_dtype_is_string_(trait, "trait")
        util.cond_check_is_string(model_name, "model_name")

        # set private variables
        self._trait = trait
        self._model_name = model_name

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def trait():
        doc = "The trait property."
        def fget(self):
            return self._trait
        def fset(self, value):
            self._trait = value
        def fdel(self):
            del self._trait
        return locals()
    trait = property(**trait())

    def ntrait():
        doc = "The ntrait property."
        def fget(self):
            return len(self._trait)
        def fset(self, value):
            error_readonly("ntrait")
        def fdel(self):
            error_readonly("ntrait")
        return locals()
    ntrait = property(**ntrait())

    def model_name():
        doc = "The model_name property."
        def fget(self):
            return self._model_name
        def fset(self, value):
            self._model_name = value
        def fdel(self):
            del self._model_name
        return locals()
    model_name = property(**model_name())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def train(self, geno, pheno):
        raise NotImplementedError("The method 'train' is abstract.")

    def predict(self, geno):
        raise NotImplementedError("The method 'predict' is abstract.")

    def reorder(self, a):
        raise NotImplementedError("The method 'reorder' is abstract.")
