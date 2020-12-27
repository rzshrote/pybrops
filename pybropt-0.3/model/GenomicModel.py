# 3rd party
import numpy

# our libraries
import pybropt.util

class GenomicModel:
    """docstring for GenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, trait, mkr_name = None, model_name = None):
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
        pybropt.util.check_is_matrix(trait, "trait")
        pybropt.util.check_matrix_ndim(trait, "trait", 1)
        pybropt.util.check_matrix_dtype_is_string_(trait, "trait")

        pybropt.util.cond_check_is_matrix(mkr_name, "mkr_name")
        pybropt.util.cond_check_matrix_ndim(mkr_name, "mkr_name", 1)
        pybropt.util.cond_check_matrix_dtype_is_string_(mkr_name, "mkr_name")

        pybropt.util.cond_check_is_string(model_name, "model_name")

        # set private variables
        self._trait = trait
        self._mkr_name = mkr_name
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
            pybropt.util.error_readonly("ntrait")
        def fdel(self):
            pybropt.util.error_readonly("ntrait")
        return locals()
    ntrait = property(**ntrait())

    def mkr_name():
        doc = "The mkr_name property."
        def fget(self):
            return self._mkr_name
        def fset(self, value):
            self._mkr_name = value
        def fdel(self):
            del self._mkr_name
        return locals()
    mkr_name = property(**mkr_name())

    def nloci():
        doc = "The nloci property."
        def fget(self):
            return len(self._mkr_name)
        def fset(self, value):
            pybropt.util.error_readonly("nloci")
        def fdel(self):
            pybropt.util.error_readonly("nloci")
        return locals()
    nloci = property(**nloci())

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

    def reorder(self, indices):
        raise NotImplementedError("The method 'reorder' is abstract.")

    def remove(self, indices):
        raise NotImplementedError("The method 'remove' is abstract.")

    def mkr_where(self, mkr_name):
        """
        Find indices for where a marker is. Indices match the order of the
        input.
        """
        # if nothing was given return nothing
        if mkr_name is None:
            return None

        if self._mkr_name is None:
            raise RuntimeError("marker names not provided")

        # build marker name location arrays
        finds = [numpy.flatnonzero(self._mkr_name == e) for e in mkr_name]

        # check for multi-mapping + build where array
        where = numpy.empty(len(finds), dtype = 'int64')
        for i,f in enumerate(finds):
            # get length of overlap numpy.ndarray
            l = len(f)

            # if length of array is > 1, we have multi-mapping; raise error
            if l > 1:
                # build error message
                s = ""
                if chr_grp is not None:
                    s += " chr_grp = " + str(chr_grp[i])
                if chr_start is not None:
                    s += " chr_start = " + str(chr_start[i])
                if chr_stop is not None:
                    s += " chr_stop = " + str(chr_stop[i])
                if mkr_name is not None:
                    s += " mkr_name = " + str(mkr_name[i])
                raise ValueError("Multi-mapping detected at:%s" % s)

            # set matrix value
            where[i] = f[0] if l > 0 else -1

        # return where array
        return where
