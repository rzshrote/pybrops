from . import DenseEstimatedBreedingValueMatrix
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_ndim

class DenseGenotypicEstimatedBreedingValueMatrix(DenseEstimatedBreedingValueMatrix):
    """docstring for GenotypicEstimatedBreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, trait = None, se = None, **kwargs):
        super(GenotypicEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            se = se,
            **kwargs
        )

        # make error checks and assignments
        # self.mat = mat            # already checked in super constructor
        # self.taxa = taxa          # already checked in super constructor
        # self.taxa_grp = taxa_grp  # already checked in super constructor
        # self.trait = trait        # already checked in super constructor
        self.se = se

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
            pybropt.util.cond_check_matrix_dtype(value, "se", numpy.float64)
            cond_check_ndarray_ndim(value, "se", 2)
            pybropt.util.cond_check_matrix_axis_len(value, "se", 0, self._mat.shape[0])
            pybropt.util.cond_check_matrix_axis_len(value, "se", 1, self._mat.shape[1])
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
            if self._taxa is not None:
                self._taxa = self._taxa[indices]                # reorder taxa array
            if self._taxa_grp is not None:
                self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array
            if self._se is not None:
                self._se = self._se[indices,:]                  # reorder standard error
        ########################################################################
        elif axis == 1:                                         ### LOCUS AXIS
            self._mat = self._mat[:,indices]                    # reorder mat array
            if self._trait is not None:
                self._trait = self._trait[indices]              # reorder trait array
            if self._se is not None:
                self._se = self._se[:,indices]                  # reorder standard error
