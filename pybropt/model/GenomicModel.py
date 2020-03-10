class GenomicModel:
    """docstring for GenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################

    def __init__(self, trait, model_name = None):
        # error check the input
        check_is_matrix(trait, "trait")
        check_matrix_dtype_is_object_(trait, "trait")
        check_matrix_ndim(trait, "trait", 1)

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
            return len(trait)
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
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def train(self, geno, pheno):
        raise NotImplementedError("The method 'train' is not implemented.")

    @classmethod
    def predict(self, geno):
        raise NotImplementedError("The method 'predict' is not implemented.")

    @classmethod
    def reorder(self, a):
        raise NotImplementedError("This method is not implemented yet.")
