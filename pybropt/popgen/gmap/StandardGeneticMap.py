from . import GeneticMap
import pybropt.util

class StandardGeneticMap(GeneticMap):
    """docstring for StandardGeneticMap."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, vrnt_chrgrp, vrnt_phypos, vrnt_genpos, **kwargs):
        super(StandardGeneticMap, self).__init__(**kwargs)
        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_genpos = vrnt_genpos
        # TODO: check all lengths equivalent

    def __len__(self):
        return len(self._vrnt_genpos)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def vrnt_chrgrp():
        doc = "The vrnt_chrgrp property."
        def fget(self):
            return self._vrnt_chrgrp
        def fset(self, value):
            pybropt.util.check_is_matrix(value, "vrnt_chrgrp")
            pybropt.util.check_matrix_dtype(value, "vrnt_chrgrp", numpy.int64)
            pybropt.util.check_matrix_ndim(value, "vrnt_chrgrp", 1)
            self._vrnt_chrgrp = value
        def fdel(self):
            del self._vrnt_chrgrp
        return locals()
    vrnt_chrgrp = property(**vrnt_chrgrp())

    def vrnt_phypos():
        doc = "The vrnt_phypos property."
        def fget(self):
            return self._vrnt_phypos
        def fset(self, value):
            pybropt.util.check_is_matrix(value, "vrnt_phypos")
            pybropt.util.check_matrix_dtype(value, "vrnt_phypos", numpy.int64)
            pybropt.util.check_matrix_ndim(value, "vrnt_phypos", 1)
            self._vrnt_phypos = value
        def fdel(self):
            del self._vrnt_phypos
        return locals()
    vrnt_phypos = property(**vrnt_phypos())

    def vrnt_genpos():
        doc = "The vrnt_genpos property."
        def fget(self):
            return self._vrnt_genpos
        def fset(self, value):
            pybropt.util.check_is_matrix(value, "vrnt_genpos")
            pybropt.util.check_matrix_dtype(value, "vrnt_genpos", numpy.float64)
            pybropt.util.check_matrix_ndim(value, "vrnt_genpos", 1)
            self._vrnt_genpos = value
        def fdel(self):
            del self._vrnt_genpos
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################