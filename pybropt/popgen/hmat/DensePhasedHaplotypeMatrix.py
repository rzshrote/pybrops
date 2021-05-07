# TODO: finish writing this class
from . import DenseHaplotypeMatrix
from pybropt.core.error import check_is_ndarray

class DensePhasedHaplotypeMatrix(DenseHaplotypeMatrix):
    """docstring for DensePhasedHaplotypeMatrix ."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        """
        DensePhasedHaplotypeMatrix constructor

        Parameters
        ----------
        haplo : numpy.ndarray
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DensePhasedHaplotypeMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            # TODO: add more checks
            check_is_ndarray(value, "mat")
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())

    def ploidy():
        doc = "The ploidy property."
        def fget(self):
            return self._ploidy
        def fset(self, value):
            self._ploidy = value
        def fdel(self):
            del self._ploidy
        return locals()
    ploidy = property(**ploidy())

    def ntaxa():
        doc = "The ntaxa property."
        def fget(self):
            return self._ntaxa
        def fset(self, value):
            self._ntaxa = value
        def fdel(self):
            del self._ntaxa
        return locals()
    ntaxa = property(**ntaxa())

    def nloci():
        doc = "The nloci property."
        def fget(self):
            return self._nloci
        def fset(self, value):
            self._nloci = value
        def fdel(self):
            del self._nloci
        return locals()
    nloci = property(**nloci())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
