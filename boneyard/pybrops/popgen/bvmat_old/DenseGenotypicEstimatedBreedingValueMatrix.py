import copy
import numpy

from . import DenseEstimatedBreedingValueMatrix

from pybrops.core.error import cond_check_is_ndarray
from pybrops.core.error import cond_check_ndarray_axis_len
from pybrops.core.error import cond_check_ndarray_dtype_is_float64
from pybrops.core.error import cond_check_ndarray_ndim

# WARNING: under construction!!!
class DenseGenotypicEstimatedBreedingValueMatrix(DenseEstimatedBreedingValueMatrix):
    """docstring for GenotypicEstimatedBreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, raw = None, se = None, trait = None, taxa = None, taxa_grp = None, **kwargs: dict):
        super(DenseGenotypicEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            raw = raw,
            trait = trait,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )
        self.se = se

    def __copy__(self):
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        # construct new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            raw = copy.copy(self.raw),
            se = copy.copy(self.se),
            trait = copy.copy(self.trait),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
        )

        # copy metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        return out

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Matrix
        """
        # construct new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat),
            raw = copy.deepcopy(self.raw),
            se = copy.deepcopy(self.se),
            trait = copy.deepcopy(self.trait),
            taxa = copy.deepcopy(self.taxa),
            taxa_grp = copy.deepcopy(self.taxa_grp),
        )

        # copy metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len)

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Breeding Value Data ##################
    def se():
        doc = "The se property."
        def fget(self):
            return self._se
        def fset(self, value):
            cond_check_is_ndarray(value, "se")
            cond_check_ndarray_dtype_is_float64(value, "se")
            cond_check_ndarray_ndim(value, "se", 2)
            cond_check_ndarray_axis_len(value, "se", 0, self._mat.shape[0])
            cond_check_ndarray_axis_len(value, "se", 1, self._mat.shape[1])
            self._se = value
        def fdel(self):
            del self._se
        return locals()
    se = property(**se())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def reorder(self, indices, axis = -1):
        """
        Reorder the VariantMatrix.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        axis : int
            The axis over which to reorder values.

        """
        # transform axis number to an index
        axis = self.axis_index(axis)

        ########################################################################
        if axis == 0:                                           ### TAXA AXIS
            self._mat = self._mat[indices,:]                    # reorder mat array
            if self._raw is not None:
                self._raw = self._raw[:,indices,:]
            if self._taxa is not None:
                self._taxa = self._taxa[indices]                # reorder taxa array
            if self._taxa_grp is not None:
                self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array
            if self._se is not None:
                self._se = self._se[indices,:]                  # reorder standard error
        ########################################################################
        elif axis == 1:                                         ### LOCUS AXIS
            self._mat = self._mat[:,indices]                    # reorder mat array
            if self._raw is not None:
                self._raw = self._raw[:,:,indices]
            if self._trait is not None:
                self._trait = self._trait[indices]              # reorder trait array
            if self._se is not None:
                self._se = self._se[:,indices]                  # reorder standard error



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGenotypicEstimatedBreedingValueMatrix(v):
    return isinstance(v, DenseGenotypicEstimatedBreedingValueMatrix)

def check_is_DenseGenotypicEstimatedBreedingValueMatrix(v, varname):
    if not isinstance(v, DenseGenotypicEstimatedBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseGenotypicEstimatedBreedingValueMatrix." % varname)

def cond_check_is_DenseGenotypicEstimatedBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseGenotypicEstimatedBreedingValueMatrix(v, varname)
