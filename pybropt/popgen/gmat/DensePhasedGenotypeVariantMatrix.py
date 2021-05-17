# import 3rd party modules we'll need
import copy
import cyvcf2
import numpy

# import our libraries
from . import GenotypeVariantMatrix
from . import DensePhasedGenotypeMatrix

from pybropt.core.mat import get_axis
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.error import check_ndarray_dtype
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_ndim
from pybropt.core.error import cond_check_ndarray_dtype
from pybropt.core.error import cond_check_ndarray_axis_len
from pybropt.core.error import check_ndarray_axis_len
from pybropt.core.error import cond_check_ndarray_dtype_is_object
from pybropt.core.error import check_is_iterable

from pybropt.popgen.gmap import check_is_GeneticMap
from pybropt.popgen.gmap import check_is_GeneticMapFunction

class DensePhasedGenotypeVariantMatrix(DensePhasedGenotypeMatrix,GenotypeVariantMatrix):
    """docstring for PhasedVariantMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, vrnt_chrgrp, vrnt_phypos, taxa = None, taxa_grp = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        PhasedGenotypeVariantMatrix object constructor.

        Parameters
        ----------
        mat : numpy.int8
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
        vrnt_chrgrp : numpy.ndarray
            A int64 chromosome group array of shape (p).
            Where:
                'p' is the number of markers.
        vrnt_phypos : numpy.ndarray
        vrnt_name : numpy.ndarray
        vrnt_genpos : numpy.ndarray
        vrnt_xoprob : numpy.ndarray
        vrnt_hapgrp : numpy.ndarray
        vrnt_mask : numpy.ndarray
        taxa : numpy.ndarray
        taxa_grp : numpy.ndarray
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        # call all parent constructors (some arguments not used)
        super(DensePhasedGenotypeVariantMatrix, self).__init__(
            mat = mat
        )

        # set variables
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_name = vrnt_name
        self.vrnt_genpos = vrnt_genpos
        self.vrnt_xoprob = vrnt_xoprob
        self.vrnt_hapgrp = vrnt_hapgrp
        self.vrnt_mask = vrnt_mask

        # set sort metadata to None
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len = None

    #################### Matrix copying ####################
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
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            vrnt_chrgrp = copy.copy(self.vrnt_chrgrp),
            vrnt_phypos = copy.copy(self.vrnt_phypos),
            vrnt_name = copy.copy(self.vrnt_name),
            vrnt_genpos = copy.copy(self.vrnt_genpos),
            vrnt_xoprob = copy.copy(self.vrnt_xoprob),
            vrnt_hapgrp = copy.copy(self.vrnt_hapgrp),
            vrnt_mask = copy.copy(self.vrnt_mask),
        )

        # copy metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)
        out.vrnt_chrgrp_name = copy.copy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.copy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.copy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len = copy.copy(self.vrnt_chrgrp_len)

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
            taxa = copy.deepcopy(self.taxa),
            taxa_grp = copy.deepcopy(self.taxa_grp),
            vrnt_chrgrp = copy.deepcopy(self.vrnt_chrgrp),
            vrnt_phypos = copy.deepcopy(self.vrnt_phypos),
            vrnt_name = copy.deepcopy(self.vrnt_name),
            vrnt_genpos = copy.deepcopy(self.vrnt_genpos),
            vrnt_xoprob = copy.deepcopy(self.vrnt_xoprob),
            vrnt_hapgrp = copy.deepcopy(self.vrnt_hapgrp),
            vrnt_mask = copy.deepcopy(self.vrnt_mask),
        )

        # copy metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len)
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len = copy.deepcopy(self.vrnt_chrgrp_len)

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Taxa Data Properites #################
    def taxa():
        doc = "The taxa property."
        def fget(self):
            return self._taxa
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa")
            cond_check_ndarray_dtype_is_object(value, "taxa")
            cond_check_ndarray_ndim(value, "taxa", 1)
            cond_check_ndarray_axis_len(value, "taxa", 0, self._mat.shape[1])
            self._taxa = value
        def fdel(self):
            del self._taxa
        return locals()
    taxa = property(**taxa())

    def taxa_grp():
        doc = "The taxa_grp property."
        def fget(self):
            return self._taxa_grp
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp")
            cond_check_ndarray_dtype(value, "taxa_grp", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp", 1)
            cond_check_ndarray_axis_len(value, "taxa_grp", 0, self._mat.shape[1])
            self._taxa_grp = value
        def fdel(self):
            del self._taxa_grp
        return locals()
    taxa_grp = property(**taxa_grp())

    ############### Taxa Metadata Properites ###############
    def taxa_grp_name():
        doc = "The taxa_grp_name property."
        def fget(self):
            return self._taxa_grp_name
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_name")
            cond_check_ndarray_dtype(value, "taxa_grp_name", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_name", 1)
            self._taxa_grp_name = value
        def fdel(self):
            del self._taxa_grp_name
        return locals()
    taxa_grp_name = property(**taxa_grp_name())

    def taxa_grp_stix():
        doc = "The taxa_grp_stix property."
        def fget(self):
            return self._taxa_grp_stix
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_stix")
            cond_check_ndarray_dtype(value, "taxa_grp_stix", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_stix", 1)
            self._taxa_grp_stix = value
        def fdel(self):
            del self._taxa_grp_stix
        return locals()
    taxa_grp_stix = property(**taxa_grp_stix())

    def taxa_grp_spix():
        doc = "The taxa_grp_spix property."
        def fget(self):
            return self._taxa_grp_spix
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_spix")
            cond_check_ndarray_dtype(value, "taxa_grp_spix", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_spix", 1)
            self._taxa_grp_spix = value
        def fdel(self):
            del self._taxa_grp_spix
        return locals()
    taxa_grp_spix = property(**taxa_grp_spix())

    def taxa_grp_len():
        doc = "The taxa_grp_len property."
        def fget(self):
            return self._taxa_grp_len
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa_grp_len")
            cond_check_ndarray_dtype(value, "taxa_grp_len", numpy.int64)
            cond_check_ndarray_ndim(value, "taxa_grp_len", 1)
            self._taxa_grp_len = value
        def fdel(self):
            del self._taxa_grp_len
        return locals()
    taxa_grp_len = property(**taxa_grp_len())

    ############### Variant Data Properites ################
    def vrnt_chrgrp():
        doc = "The vrnt_chrgrp property."
        def fget(self):
            return self._vrnt_chrgrp
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_chrgrp")
            cond_check_ndarray_dtype(value, "vrnt_chrgrp", numpy.int64)
            cond_check_ndarray_ndim(value, "vrnt_chrgrp", 1)
            cond_check_ndarray_axis_len(value, "vrnt_chrgrp", 0, self._mat.shape[2])
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
            cond_check_is_ndarray(value, "vrnt_phypos")
            cond_check_ndarray_dtype(value, "vrnt_phypos", numpy.int64)
            cond_check_ndarray_ndim(value, "vrnt_phypos", 1)
            cond_check_ndarray_axis_len(value, "vrnt_phypos", 0, self._mat.shape[2])
            self._vrnt_phypos = value
        def fdel(self):
            del self._vrnt_phypos
        return locals()
    vrnt_phypos = property(**vrnt_phypos())

    def vrnt_name():
        doc = "The vrnt_name property."
        def fget(self):
            return self._vrnt_name
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_name")
            cond_check_ndarray_dtype_is_object(value, "vrnt_name")
            cond_check_ndarray_ndim(value, "vrnt_name", 1)
            cond_check_ndarray_axis_len(value, "vrnt_name", 0, self._mat.shape[2])
            self._vrnt_name = value
        def fdel(self):
            del self._vrnt_name
        return locals()
    vrnt_name = property(**vrnt_name())

    def vrnt_genpos():
        doc = "The vrnt_genpos property."
        def fget(self):
            return self._vrnt_genpos
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_genpos")
            cond_check_ndarray_dtype(value, "vrnt_genpos", numpy.float64)
            cond_check_ndarray_ndim(value, "vrnt_genpos", 1)
            cond_check_ndarray_axis_len(value, "vrnt_genpos", 0, self._mat.shape[2])
            self._vrnt_genpos = value
        def fdel(self):
            del self._vrnt_genpos
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    def vrnt_xoprob():
        doc = "The vrnt_xoprob property."
        def fget(self):
            return self._vrnt_xoprob
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_xoprob")
            cond_check_ndarray_dtype(value, "vrnt_xoprob", numpy.float64)
            cond_check_ndarray_ndim(value, "vrnt_xoprob", 1)
            cond_check_ndarray_axis_len(value, "vrnt_xoprob", 0, self._mat.shape[2])
            self._vrnt_xoprob = value
        def fdel(self):
            del self._vrnt_xoprob
        return locals()
    vrnt_xoprob = property(**vrnt_xoprob())

    def vrnt_hapgrp():
        doc = "The vrnt_hapgrp property."
        def fget(self):
            return self._vrnt_hapgrp
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_hapgrp")
            cond_check_ndarray_dtype(value, "vrnt_hapgrp", numpy.int64)
            cond_check_ndarray_ndim(value, "vrnt_hapgrp", 1)
            cond_check_ndarray_axis_len(value, "vrnt_hapgrp", 0, self._mat.shape[2])
            self._vrnt_hapgrp = value
        def fdel(self):
            del self._vrnt_hapgrp
        return locals()
    vrnt_hapgrp = property(**vrnt_hapgrp())

    def vrnt_mask():
        doc = "The vrnt_mask property."
        def fget(self):
            return self._vrnt_mask
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_mask")
            cond_check_ndarray_dtype(value, "vrnt_mask", numpy.bool_)
            cond_check_ndarray_ndim(value, "vrnt_mask", 1)
            cond_check_ndarray_axis_len(value, "vrnt_mask", 0, self._mat.shape[2])
            self._vrnt_mask = value
        def fdel(self):
            del self._vrnt_mask
        return locals()
    vrnt_mask = property(**vrnt_mask())

    ############# Variant Metadata Properites ##############
    def vrnt_chrgrp_name():
        doc = "The vrnt_chrgrp_name property."
        def fget(self):
            return self._vrnt_chrgrp_name
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_chrgrp_name")
            cond_check_ndarray_dtype(value, "vrnt_chrgrp_name", numpy.int64)
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_name", 1)
            self._vrnt_chrgrp_name = value
        def fdel(self):
            del self._vrnt_chrgrp_name
        return locals()
    vrnt_chrgrp_name = property(**vrnt_chrgrp_name())

    def vrnt_chrgrp_stix():
        doc = "The vrnt_chrgrp_stix property."
        def fget(self):
            return self._vrnt_chrgrp_stix
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_chrgrp_stix")
            cond_check_ndarray_dtype(value, "vrnt_chrgrp_stix", numpy.int64)
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_stix", 1)
            self._vrnt_chrgrp_stix = value
        def fdel(self):
            del self._vrnt_chrgrp_stix
        return locals()
    vrnt_chrgrp_stix = property(**vrnt_chrgrp_stix())

    def vrnt_chrgrp_spix():
        doc = "The vrnt_chrgrp_spix property."
        def fget(self):
            return self._vrnt_chrgrp_spix
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_chrgrp_spix")
            cond_check_ndarray_dtype(value, "vrnt_chrgrp_spix", numpy.int64)
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_spix", 1)
            self._vrnt_chrgrp_spix = value
        def fdel(self):
            del self._vrnt_chrgrp_spix
        return locals()
    vrnt_chrgrp_spix = property(**vrnt_chrgrp_spix())

    def vrnt_chrgrp_len():
        doc = "The vrnt_chrgrp_len property."
        def fget(self):
            return self._vrnt_chrgrp_len
        def fset(self, value):
            cond_check_is_ndarray(value, "vrnt_chrgrp_len")
            cond_check_ndarray_dtype(value, "vrnt_chrgrp_len", numpy.int64)
            cond_check_ndarray_ndim(value, "vrnt_chrgrp_len", 1)
            self._vrnt_chrgrp_len = value
        def fdel(self):
            del self._vrnt_chrgrp_len
        return locals()
    vrnt_chrgrp_len = property(**vrnt_chrgrp_len())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DensePhasedGenotypeMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedGenotypeMatrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        # BUG: need to check that values are compatible in the first place
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DensePhasedGenotypeVariantMatrix(values):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DensePhasedGenotypeVariantMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if axis == 0:
            raise ValueError("adjoin along axis 0 not supported") ## TODO: implement me
        elif axis == 1:
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot adjoin: taxa_grp argument is required")
        elif axis == 2:
            if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
                raise TypeError("cannot adjoin: vrnt_chrgrp argument is required")
            if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
                raise TypeError("cannot adjoin: vrnt_phypos argument is required")
            if (self._vrnt_name is not None) and (vrnt_name is None):
                vrnt_name = numpy.object_([None] * values.shape[2])     # fill with None
            if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
                raise TypeError("cannot adjoin: vrnt_genpos argument is required")
            if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
                raise TypeError("cannot adjoin: vrnt_xoprob argument is required")
            if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
                raise TypeError("cannot adjoin: vrnt_hapgrp argument is required")
            if (self._vrnt_mask is not None) and (vrnt_mask is None):
                raise TypeError("cannot adjoin: vrnt_mask argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # adjoin values
        values = numpy.append(self._mat, values, axis = axis)
        if axis == 1:
            if self._taxa is not None:
                taxa = numpy.append(self._taxa, taxa, axis = 0)
            if self._taxa_grp is not None:
                taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)
        elif axis == 2:
            if self._vrnt_chrgrp is not None:
                vrnt_chrgrp = numpy.append(self._vrnt_chrgrp, vrnt_chrgrp, axis = 0)
            if self._vrnt_phypos is not None:
                vrnt_phypos = numpy.append(self._vrnt_phypos, vrnt_phypos, axis = 0)
            if self._vrnt_name is not None:
                vrnt_name = numpy.append(self._vrnt_name, vrnt_name, axis = 0)
            if self._vrnt_genpos is not None:
                vrnt_genpos = numpy.append(self._vrnt_genpos, vrnt_genpos, axis = 0)
            if self._vrnt_xoprob is not None:
                vrnt_xoprob = numpy.append(self._vrnt_xoprob, vrnt_xoprob, axis = 0)
            if self._vrnt_hapgrp is not None:
                vrnt_hapgrp = numpy.append(self._vrnt_hapgrp, vrnt_hapgrp, axis = 0)
            if self._vrnt_mask is not None:
                vrnt_mask = numpy.append(self._vrnt_mask, vrnt_mask, axis = 0)

        out = self.__class__(
            mat = values,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        return out

    def adjoin_taxa(self, values, taxa = None, taxa_grp = None, **kwargs):
        """
        Add additional elements to the end of the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        return self.adjoin(
            values = values,
            axis = 1,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def adjoin_vrnt(self, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Add additional elements to the end of the Matrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        return self.adjoin(
            values = values,
            axis = 2,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

    def delete(self, obj, axis = -1, **kwargs):
        """
        Delete sub-arrays along an axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to delete the subarray defined by obj.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        if axis == 0:
            raise ValueError("delete along axis 0 not supported") ## TODO: implement me

        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp
        vrnt_chrgrp = self._vrnt_chrgrp
        vrnt_phypos = self._vrnt_phypos
        vrnt_name = self._vrnt_name
        vrnt_genpos = self._vrnt_genpos
        vrnt_xoprob = self._vrnt_xoprob
        vrnt_hapgrp = self._vrnt_hapgrp
        vrnt_mask = self._vrnt_mask

        # delete values
        mat = numpy.delete(mat, obj, axis = axis)
        if axis == 1:
            if taxa is not None:
                taxa = numpy.delete(taxa, obj, axis = 0)
            if taxa_grp is not None:
                taxa_grp = numpy.delete(taxa_grp, obj, axis = 0)
        elif axis == 2:
            if vrnt_chrgrp is not None:
                vrnt_chrgrp = numpy.delete(vrnt_chrgrp, obj, axis = 0)
            if vrnt_phypos is not None:
                vrnt_phypos = numpy.delete(vrnt_phypos, obj, axis = 0)
            if vrnt_name is not None:
                vrnt_name = numpy.delete(vrnt_name, obj, axis = 0)
            if vrnt_genpos is not None:
                vrnt_genpos = numpy.delete(vrnt_genpos, obj, axis = 0)
            if vrnt_xoprob is not None:
                vrnt_xoprob = numpy.delete(vrnt_xoprob, obj, axis = 0)
            if vrnt_hapgrp is not None:
                vrnt_hapgrp = numpy.delete(vrnt_hapgrp, obj, axis = 0)
            if vrnt_mask is not None:
                vrnt_mask = numpy.delete(vrnt_mask, obj, axis = 0)

        out = self.__class__(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        return out

    def delete_taxa(self, obj, **kwargs):
        """
        Delete sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.delete(
            obj = obj,
            axis = 1,
            **kwargs
        )

    def delete_vrnt(self, obj, **kwargs):
        """
        Delete sub-arrays along the variant axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.delete(
            obj = obj,
            axis = 2,
            **kwargs
        )

    def insert(self, obj, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : DensePhasedGenotypeVariantMatrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
            If values is a DensePhasedGenotypeVariantMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedGenotypeVariantMatrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # BUG: need to check that values are compatible in the first place
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DensePhasedGenotypeVariantMatrix(values):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DensePhasedGenotypeVariantMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if axis == 0:
            raise ValueError("insert along axis 0 not supported") ## TODO: implement me
        elif axis == 1:
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot insert: taxa_grp argument is required")
        elif axis == 2:
            if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
                raise TypeError("cannot insert: vrnt_chrgrp argument is required")
            if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
                raise TypeError("cannot insert: vrnt_phypos argument is required")
            if (self._vrnt_name is not None) and (vrnt_name is None):
                vrnt_name = numpy.object_([None] * values.shape[2])     # fill with None
            if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
                raise TypeError("cannot insert: vrnt_genpos argument is required")
            if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
                raise TypeError("cannot insert: vrnt_xoprob argument is required")
            if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
                raise TypeError("cannot insert: vrnt_hapgrp argument is required")
            if (self._vrnt_mask is not None) and (vrnt_mask is None):
                raise TypeError("cannot insert: vrnt_mask argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # insert values
        values = numpy.insert(self._mat, obj, values, axis = axis)
        if axis == 1:
            if self._taxa is not None:
                taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
            if self._taxa_grp is not None:
                taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)
        elif axis == 2:
            if self._vrnt_chrgrp is not None:
                vrnt_chrgrp = numpy.insert(self._vrnt_chrgrp, obj, vrnt_chrgrp, axis = 0)
            if self._vrnt_phypos is not None:
                vrnt_phypos = numpy.insert(self._vrnt_phypos, obj, vrnt_phypos, axis = 0)
            if self._vrnt_name is not None:
                vrnt_name = numpy.insert(self._vrnt_name, obj, vrnt_name, axis = 0)
            if self._vrnt_genpos is not None:
                vrnt_genpos = numpy.insert(self._vrnt_genpos, obj, vrnt_genpos, axis = 0)
            if self._vrnt_xoprob is not None:
                vrnt_xoprob = numpy.insert(self._vrnt_xoprob, obj, vrnt_xoprob, axis = 0)
            if self._vrnt_hapgrp is not None:
                vrnt_hapgrp = numpy.insert(self._vrnt_hapgrp, obj, vrnt_hapgrp, axis = 0)
            if self._vrnt_mask is not None:
                vrnt_mask = numpy.insert(self._vrnt_mask, obj, vrnt_mask, axis = 0)

        # create output
        out = self.__class__(
            mat = values,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        return out

    def insert_taxa(self, obj, values, taxa = None, taxa_grp = None, **kwargs):
        """
        Insert values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.insert(
            obj = obj,
            values = values,
            axis = 1,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def insert_vrnt(self, obj, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Insert values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        return self.insert(
            obj = obj,
            values = values,
            axis = 2,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

    def select(self, indices, axis = -1, **kwargs):
        """
        Select certain values from the matrix.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        axis : int
            The axis along which values are selected.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        if axis == 0:
            raise ValueError("select along axis 0 not supported") ## TODO: implement me

        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp
        vrnt_chrgrp = self._vrnt_chrgrp
        vrnt_phypos = self._vrnt_phypos
        vrnt_name = self._vrnt_name
        vrnt_genpos = self._vrnt_genpos
        vrnt_xoprob = self._vrnt_xoprob
        vrnt_hapgrp = self._vrnt_hapgrp
        vrnt_mask = self._vrnt_mask

        # select values
        mat = numpy.take(mat, indices, axis = axis)
        if axis == 1:
            if taxa is not None:
                taxa = numpy.take(taxa, indices, axis = 0)
            if taxa_grp is not None:
                taxa_grp = numpy.take(taxa_grp, indices, axis = 0)
        elif axis == 2:
            if vrnt_chrgrp is not None:
                vrnt_chrgrp = numpy.take(vrnt_chrgrp, indices, axis = 0)
            if vrnt_phypos is not None:
                vrnt_phypos = numpy.take(vrnt_phypos, indices, axis = 0)
            if vrnt_name is not None:
                vrnt_name = numpy.take(vrnt_name, indices, axis = 0)
            if vrnt_genpos is not None:
                vrnt_genpos = numpy.take(vrnt_genpos, indices, axis = 0)
            if vrnt_xoprob is not None:
                vrnt_xoprob = numpy.take(vrnt_xoprob, indices, axis = 0)
            if vrnt_hapgrp is not None:
                vrnt_hapgrp = numpy.take(vrnt_hapgrp, indices, axis = 0)
            if vrnt_mask is not None:
                vrnt_mask = numpy.take(vrnt_mask, indices, axis = 0)

        out = self.__class__(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        return out

    def select_taxa(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the taxa axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        return self.select(
            indices = indices,
            axis = 1,
            **kwargs
        )

    def select_vrnt(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the variant axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        return self.select(
            indices = indices,
            axis = 2,
            **kwargs
        )

    @staticmethod
    def concat(mats, axis = -1, **kwargs):
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : array_like of matrices
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : DensePhasedGenotypeVariantMatrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        # ensure that we have an iterable object
        check_is_iterable(mats, "mats")

        # get length of mats
        mats_len = len(mats)

        # ensure that we have an array_like of length >= 1
        if mats_len <= 0:
            raise ValueError("need at least one Matrix to concatenate")

        # ensure that all items in mats are DensePhasedGenotypeVariantMatrix
        for i in range(mats_len):
            check_is_DensePhasedGenotypeVariantMatrix(mats[i], "mats[{0}]".format(i))

        # get first matrix
        mats0 = mats[0]

        # get axis
        axis = get_axis(axis, mats0.mat.ndim)

        # extract tuples of shape parameters for Matrix
        nphase_t, ntaxa_t, nloci_t = zip(*[m.mat.shape for m in mats])

        # extract first Matrix shape parameters
        nphase, ntaxa, nloci = mats0.mat.shape

        # create matrix lists
        mat_l = [m.mat for m in mats]
        taxa_l = None
        taxa_grp_l = None
        vrnt_chrgrp_l = None
        vrnt_phypos_l = None
        vrnt_name_l = None
        vrnt_genpos_l = None
        vrnt_xoprob_l = None
        vrnt_hapgrp_l = None
        vrnt_mask_l = None

        # check shapes and add to list
        if axis == 0:
            raise ValueError("concat along axis 0 not supported") # TODO: implement me
        elif axis == 1:                                         # concatenate additional taxa
            # check matrix shapes
            if any(e != nphase for e in nphase_t):              # raise error if any have different phase number
                raise ValueError("Matrix shapes do not all align along axis 0 (phase axis)")
            if any(e != nloci for e in nloci_t):                # raise error if any have different loci number
                raise ValueError("Matrix shapes do not all align along axis 2 (loci axis)")
            # add taxa related attributes to lists
            if mats0.taxa is not None:                        # populate taxa_l
                taxa_l = [numpy.object_([None]*m.ntaxa) if m.taxa is None else m.taxa for m in mats]
            if mats0.taxa_grp is not None:
                taxa_grp_l = [m.taxa_grp for m in mats]
                if any(e is None for e in taxa_grp_l):
                    raise ValueError("cannot concat: taxa_grp needed for all Matrix in list")
        elif axis == 2:                                         # concatenate additional loci
            # check matrix shapes
            if any(e != nphase for e in nphase_t):              # raise error if any have different phase number
                raise ValueError("Matrix shapes do not all align along axis 0 (phase axis)")
            if any(e != ntaxa for e in ntaxa_t):                # raise error if any have different taxa number
                raise ValueError("Matrix shapes do not all align along axis 1 (taxa axis)")
            # add loci related attributes to lists
            if mats0.vrnt_chrgrp is not None:
                vrnt_chrgrp_l = [m.vrnt_chrgrp for m in mats]
                if any(e is None for e in vrnt_chrgrp_l):
                    raise ValueError("cannot concat: vrnt_chrgrp needed for all Matrix in list")
            if mats0.vrnt_phypos is not None:
                vrnt_phypos_l = [m.vrnt_phypos for m in mats]
                if any(e is None for e in vrnt_phypos_l):
                    raise ValueError("cannot concat: vrnt_phypos needed for all Matrix in list")
            if mats0.vrnt_name is not None:
                vrnt_name_l = [numpy.object_([None]*m.nloci) if m.vrnt_name is None else m.vrnt_name for m in mats]
            if mats0.vrnt_genpos is not None:
                vrnt_genpos_l = [m.vrnt_genpos for m in mats]
                if any(e is None for e in vrnt_genpos_l):
                    raise ValueError("cannot concat: vrnt_genpos needed for all Matrix in list")
            if mats0.vrnt_xoprob is not None:
                vrnt_xoprob_l = [m.vrnt_xoprob for m in mats]
                if any(e is None for e in vrnt_xoprob_l):
                    raise ValueError("cannot concat: vrnt_xoprob needed for all Matrix in list")
            if mats0.vrnt_hapgrp is not None:
                vrnt_hapgrp_l = [m.vrnt_hapgrp for m in mats]
                if any(e is None for e in vrnt_hapgrp_l):
                    raise ValueError("cannot concat: vrnt_hapgrp needed for all Matrix in list")
            if mats0.vrnt_mask is not None:
                vrnt_mask_l = [m.vrnt_mask for m in mats]
                if any(e is None for e in vrnt_mask_l):
                    raise ValueError("cannot concat: vrnt_mask needed for all Matrix in list")

        # concatenate everything
        mat = numpy.concatenate(mat_l, axis = axis)
        vrnt_chrgrp = mats0.vrnt_chrgrp if vrnt_chrgrp_l is None else numpy.concatenate(vrnt_chrgrp_l, axis = 0)
        vrnt_phypos = mats0.vrnt_phypos if vrnt_phypos_l is None else numpy.concatenate(vrnt_phypos_l, axis = 0)
        taxa = mats0.taxa if taxa_l is None else numpy.concatenate(taxa_l, axis = 0)
        taxa_grp = mats0.taxa_grp if taxa_grp_l is None else numpy.concatenate(taxa_grp_l, axis = 0)
        vrnt_name = mats0.vrnt_name if vrnt_name_l is None else numpy.concatenate(vrnt_name_l, axis = 0)
        vrnt_genpos = mats0.vrnt_genpos if vrnt_genpos_l is None else numpy.concatenate(vrnt_genpos_l, axis = 0)
        vrnt_xoprob = mats0.vrnt_xoprob if vrnt_xoprob_l is None else numpy.concatenate(vrnt_xoprob_l, axis = 0)
        vrnt_hapgrp = mats0.vrnt_hapgrp if vrnt_hapgrp_l is None else numpy.concatenate(vrnt_hapgrp_l, axis = 0)
        vrnt_mask = mats0.vrnt_mask if vrnt_mask_l is None else numpy.concatenate(vrnt_mask_l, axis = 0)

        # concatenate everything and put into new DensePhasedGenotypeVariantMatrix
        out = DensePhasedGenotypeVariantMatrix(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        return out

    @staticmethod
    def concat_taxa(mats, **kwargs):
        """
        Concatenate list of Matrix together along the taxa axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        return DensePhasedGenotypeVariantMatrix.concat(
            mats = mats,
            axis = 1,
            **kwargs
        )

    @staticmethod
    def concat_vrnt(mats, **kwargs):
        """
        Concatenate list of Matrix together along the variant axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        return DensePhasedGenotypeVariantMatrix.concat(
            mats = mats,
            axis = 2,
            **kwargs
        )

    #################### Matrix pruning ####################
    def prune(self, axis = -1, nt = None, M = None, **kwargs):
        """
        Prune markers evenly across all chromosomes.

        Parameters
        ----------
        nt : int
            Target distance between each selected marker in nucleotides.
        M : float
            Target distance between each selected marker in Morgans.
            If this option is specified, selection based on Morgans takes first
            priority. If the physical distance between two markers selected
            based on their genetic distance exceeds 'nt' (if provided), the
            additional markers are sought between those regions.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # raise error if axis != 2 (variant axis)
        if axis != 2:
            raise ValueError("pruning not applicable along axis {0}".format(axis))

        # check if we have acceptible inputs
        if (nt is None) and (M is None):
            raise ValueError("'nt' and 'M' cannot both be None")

        # if not sorted and grouped, sort and group
        if not self.is_grouped_vrnt():
            self.group_vrnt()

        # make empty index list to store selected marker indices
        indices = []

        # make a generic array pointer to a position; this position can be a
        # physical position (self._vrnt_phypos) or a genetic position
        # (self._vrnt_genpos); this is used in initial marker selection below.
        # genetic position takes dominance
        position = self._vrnt_genpos if M is not None else self._vrnt_phypos

        # generic spacing variable
        spacing = M if M is not None else nt

        # for each chromosome
        for st,sp in zip(self._vrnt_chrgrp_stix, self._vrnt_chrgrp_spix):
            # calculate chromosome length given start, end marker positions
            dist = position[sp-1] - position[st]

            # calculate the target distance between each marker (float)
            step = dist / int(math.ceil(dist / spacing))

            # force addition of the first marker on the chromosome
            indices.append(st)

            # target site; we want markers as close to this value (float)
            target = position[st] + step

            # for each locus index in the chromosome
            for i in range(st+1, sp):
                # if position exceeds target, determine which marker to add
                if position[i] >= target:
                    # get distance between target and previous marker
                    downstream = target - position[i-1]

                    # get distance between target and current marker
                    upstream = position[i] - target

                    # determine which index to add
                    ix = i-1 if downstream < upstream else i

                    # if we haven't added this index previously, add it
                    if ix != indices[-1]:
                        indices.append(ix)

                    # increment target site position
                    target += step

            # final check to make sure we've added last marker on chromosome
            if (sp-1) != indices[-1]:
                indices.append(sp-1)

        # secondary marker selection based on 'nt' if both 'M' and 'nt' provided
        if (M is not None) and (nt is not None):
            # make new indices list to store M indices + nt indices
            new_indices = []

            # for each neighbor marker pair
            for up,down in zip(indices[:-1], indices[1:]):
                # append the upstream index
                new_indices.append(up)

                # if they are on the same chromosome
                if self._vrnt_chrgrp[up] == self._vrnt_chrgrp[down]:
                    # calculate physical distance between two selected markers
                    dist = self._vrnt_phypos[down] - self._vrnt_phypos[up]

                    # if we exceed 'nt' distance
                    if dist > nt:
                        # calculate the target distance between each marker (float)
                        step = dist / int(math.ceil(dist / nt))

                        # target site; we want markers as close to this value (float)
                        target = self._vrnt_phypos[up] + step

                        # for each locus between upstream and downstream markers
                        for i in range(up+1, down):
                            # if position exceeds target, determine which marker to add
                            if self._vrnt_phypos[i] >= target:
                                # get distance between target and previous marker
                                downstream = target - self._vrnt_phypos[i-1]

                                # get distance between target and current marker
                                upstream = self._vrnt_phypos[i] - target

                                # determine which index to add
                                ix = i-1 if downstream < upstream else i

                                # if we haven't added this index previously, add it
                                if ix != new_indices[-1]:
                                    new_indices.append(ix)

                                # increment target site position
                                target += step

            # append the last index of 'indices' to 'new_indices' since we skipped it
            new_indices.append(indices[-1])

            # replace indices with new_indices
            indices = new_indices

        # convert indices into an array
        indices = numpy.array(indices)

        # return selected indices
        return indices

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DensePhasedGenotypeVariantMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        """
        # OPTIMIZE: this is a hot mess.
        # BUG: need to check that values are compatible in the first place
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DensePhasedGenotypeVariantMatrix(values):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DensePhasedGenotypeVariantMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if axis == 0:
            raise ValueError("adjoin along axis 0 not supported") ## TODO: implement me
        elif axis == 1:
            if self._mat.shape[0] != values.shape[0]:
                raise ValueError("Matrix shapes do not all align along axis 0 (phase axis)")
            if self._mat.shape[2] != values.shape[2]:
                raise ValueError("Matrix shapes do not all align along axis 2 (loci axis)")
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot append: taxa_grp argument is required")
        elif axis == 2:
            if self._mat.shape[0] != values.shape[0]:
                raise ValueError("Matrix shapes do not all align along axis 0 (phase axis)")
            if self._mat.shape[1] != values.shape[1]:
                raise ValueError("Matrix shapes do not all align along axis 1 (taxa axis)")
            if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
                raise TypeError("cannot append: vrnt_chrgrp argument is required")
            if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
                raise TypeError("cannot append: vrnt_phypos argument is required")
            if (self._vrnt_name is not None) and (vrnt_name is None):
                vrnt_name = numpy.object_([None] * values.shape[2])     # fill with None
            if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
                raise TypeError("cannot append: vrnt_genpos argument is required")
            if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
                raise TypeError("cannot append: vrnt_xoprob argument is required")
            if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
                raise TypeError("cannot append: vrnt_hapgrp argument is required")
            if (self._vrnt_mask is not None) and (vrnt_mask is None):
                raise TypeError("cannot append: vrnt_mask argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # append values
        self._mat = numpy.append(self._mat, values, axis = axis)
        if axis == 1:
            # set fields
            if self._taxa is not None:
                self._taxa = numpy.append(self._taxa, taxa, axis = 0)
            if self._taxa_grp is not None:
                self._taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)
            # reset metadata
            self._taxa_grp_len = None
            self._taxa_grp_name = None
            self._taxa_grp_stix = None
            self._taxa_grp_spix = None
        elif axis == 2:
            # set fields
            if self._vrnt_chrgrp is not None:
                self._vrnt_chrgrp = numpy.append(self._vrnt_chrgrp, vrnt_chrgrp, axis = 0)
            if self._vrnt_phypos is not None:
                self._vrnt_phypos = numpy.append(self._vrnt_phypos, vrnt_phypos, axis = 0)
            if self._vrnt_name is not None:
                self._vrnt_name = numpy.append(self._vrnt_name, vrnt_name, axis = 0)
            if self._vrnt_genpos is not None:
                self._vrnt_genpos = numpy.append(self._vrnt_genpos, vrnt_genpos, axis = 0)
            if self._vrnt_xoprob is not None:
                self._vrnt_xoprob = numpy.append(self._vrnt_xoprob, vrnt_xoprob, axis = 0)
            if self._vrnt_hapgrp is not None:
                self._vrnt_hapgrp = numpy.append(self._vrnt_hapgrp, vrnt_hapgrp, axis = 0)
            if self._vrnt_mask is not None:
                self._vrnt_mask = numpy.append(self._vrnt_mask, vrnt_mask, axis = 0)
            # reset metadata
            self._vrnt_chrgrp_len = None
            self._vrnt_chrgrp_name = None
            self._vrnt_chrgrp_stix = None
            self._vrnt_chrgrp_spix = None

    def append_taxa(self, values, taxa = None, taxa_grp = None, **kwargs):
        """
        Append values to the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.append(
            values = values,
            axis = 1,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def append_vrnt(self, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Append values to the Matrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to append to the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to append to the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to append to the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to append to the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to append to the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to append to the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to append to the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.append(
            values = values,
            axis = 2,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

    def remove(self, obj, axis = -1, **kwargs):
        """
        Remove sub-arrays along an axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to remove the subarray defined by obj.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        if axis == 0:
            raise ValueError("delete along axis 0 not supported") ## TODO: implement me

        # delete values
        self._mat = numpy.delete(self._mat, obj, axis = axis)
        if axis == 1:
            if self._taxa is not None:
                self._taxa = numpy.delete(self._taxa, obj, axis = 0)
            if self._taxa_grp is not None:
                self._taxa_grp = numpy.delete(self._taxa_grp, obj, axis = 0)
            # reset metadata
            self._taxa_grp_len = None
            self._taxa_grp_name = None
            self._taxa_grp_stix = None
            self._taxa_grp_spix = None
        elif axis == 2:
            if self._vrnt_chrgrp is not None:
                self._vrnt_chrgrp = numpy.delete(self._vrnt_chrgrp, obj, axis = 0)
            if self._vrnt_phypos is not None:
                self._vrnt_phypos = numpy.delete(self._vrnt_phypos, obj, axis = 0)
            if self._vrnt_name is not None:
                self._vrnt_name = numpy.delete(self._vrnt_name, obj, axis = 0)
            if self._vrnt_genpos is not None:
                self._vrnt_genpos = numpy.delete(self._vrnt_genpos, obj, axis = 0)
            if self._vrnt_xoprob is not None:
                self._vrnt_xoprob = numpy.delete(self._vrnt_xoprob, obj, axis = 0)
            if self._vrnt_hapgrp is not None:
                self._vrnt_hapgrp = numpy.delete(self._vrnt_hapgrp, obj, axis = 0)
            if self._vrnt_mask is not None:
                self._vrnt_mask = numpy.delete(self._vrnt_mask, obj, axis = 0)

    def remove_taxa(self, obj, **kwargs):
        """
        Remove sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.remove(
            obj = obj,
            axis = 1,
            **kwargs
        )

    def remove_vrnt(self, obj, **kwargs):
        """
        Remove sub-arrays along the variant axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.remove(
            obj = obj,
            axis = 2,
            **kwargs
        )

    def incorp(self, obj, values, axis = -1, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Incorporate values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : array_like
            Values to incorporate into the matrix.
        axis : int
            The axis along which values are incorporated.
        **kwargs
            Additional keyword arguments.
        """
        # BUG: ensuring correct shape is needed & is currently a hot mess.
        # BUG: need to check that values are compatible in the first place
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if is_DensePhasedGenotypeVariantMatrix(values):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            if vrnt_chrgrp is None:
                vrnt_chrgrp = values.vrnt_chrgrp
            if vrnt_phypos is None:
                vrnt_phypos = values.vrnt_phypos
            if vrnt_name is None:
                vrnt_name = values.vrnt_name
            if vrnt_genpos is None:
                vrnt_genpos = values.vrnt_genpos
            if vrnt_xoprob is None:
                vrnt_xoprob = values.vrnt_xoprob
            if vrnt_hapgrp is None:
                vrnt_hapgrp = values.vrnt_hapgrp
            if vrnt_mask is None:
                vrnt_mask = values.vrnt_mask
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DensePhasedGenotypeVariantMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if axis == 0:
            raise ValueError("incorp along axis 0 not supported") ## TODO: implement me
        elif axis == 1:
            if (self._taxa is not None) and (taxa is None):
                taxa = numpy.object_([None] * values.shape[1])          # fill with None
            if (self._taxa_grp is not None) and (taxa_grp is None):
                raise TypeError("cannot incorp: taxa_grp argument is required")
        elif axis == 2:
            if (self._vrnt_chrgrp is not None) and (vrnt_chrgrp is None):
                raise TypeError("cannot incorp: vrnt_chrgrp argument is required")
            if (self._vrnt_phypos is not None) and (vrnt_phypos is None):
                raise TypeError("cannot incorp: vrnt_phypos argument is required")
            if (self._vrnt_name is not None) and (vrnt_name is None):
                vrnt_name = numpy.object_([None] * values.shape[2])     # fill with None
            if (self._vrnt_genpos is not None) and (vrnt_genpos is None):
                raise TypeError("cannot incorp: vrnt_genpos argument is required")
            if (self._vrnt_xoprob is not None) and (vrnt_xoprob is None):
                raise TypeError("cannot incorp: vrnt_xoprob argument is required")
            if (self._vrnt_hapgrp is not None) and (vrnt_hapgrp is None):
                raise TypeError("cannot incorp: vrnt_hapgrp argument is required")
            if (self._vrnt_mask is not None) and (vrnt_mask is None):
                raise TypeError("cannot incorp: vrnt_mask argument is required")

        # Remark:
        # Only test if self.field is not None.
        # Error check above guarantees that field is not None

        # OPTIMIZE: Consider merging the if statements above and below.
        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = axis)
        if axis == 1:
            if self._taxa is not None:
                self._taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
            if self._taxa_grp is not None:
                self._taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)
            # reset metadata
            self._taxa_grp_len = None
            self._taxa_grp_name = None
            self._taxa_grp_stix = None
            self._taxa_grp_spix = None
        elif axis == 2:
            if self._vrnt_chrgrp is not None:
                self._vrnt_chrgrp = numpy.insert(self._vrnt_chrgrp, obj, vrnt_chrgrp, axis = 0)
            if self._vrnt_phypos is not None:
                self._vrnt_phypos = numpy.insert(self._vrnt_phypos, obj, vrnt_phypos, axis = 0)
            if self._vrnt_name is not None:
                self._vrnt_name = numpy.insert(self._vrnt_name, obj, vrnt_name, axis = 0)
            if self._vrnt_genpos is not None:
                self._vrnt_genpos = numpy.insert(self._vrnt_genpos, obj, vrnt_genpos, axis = 0)
            if self._vrnt_xoprob is not None:
                self._vrnt_xoprob = numpy.insert(self._vrnt_xoprob, obj, vrnt_xoprob, axis = 0)
            if self._vrnt_hapgrp is not None:
                self._vrnt_hapgrp = numpy.insert(self._vrnt_hapgrp, obj, vrnt_hapgrp, axis = 0)
            if self._vrnt_mask is not None:
                self._vrnt_mask = numpy.insert(self._vrnt_mask, obj, vrnt_mask, axis = 0)

    def incorp_taxa(self, obj, values, taxa = None, taxa_grp = None, **kwargs):
        """
        Incorporate values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        taxa : numpy.ndarray
            Taxa names to incorporate into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.incorp(
            obj = obj,
            values = values,
            axis = 1,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    def incorp_vrnt(self, obj, values, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, **kwargs):
        """
        Incorporate values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to incorporate into the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to incorporate into the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to incorporate into the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to incorporate into the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to incorporate into the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to incorporate into the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to incorporate into the Matrix.
        **kwargs
            Additional keyword arguments.
        """
        self.incorp(
            obj = obj,
            values = values,
            axis = 2,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

    ################### Sorting Methods ####################
    def lexsort(self, keys = None, axis = -1, **kwargs):
        """
        Perform an indirect stable sort using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using vrnt_chrgrp as primary key, and
            vrnt_phypos as secondary key.
        axis : int
            The axis of the Matrix over which to sort values.

        Returns
        -------
        indices : numpy.ndarray
            Array of indices that sort the keys.
        """
        axis = get_axis(axis, self._mat.ndim)                   # transform axis number to an index
        emess = None                                            # error message

        if keys is None:                                        # if no keys were provided, set a default
            if axis == 0:                                       # phase axis
                keys = (None,)                                  # phase default keys
                emess = "axis unsortable by default"            # phase error message
            elif axis == 1:                                     # taxa axis
                keys = (self._taxa, self._taxa_grp)             # taxa default keys
                emess = "taxa, taxa_grp are None"               # taxa error message
            elif axis == 2:                                     # loci axis
                keys = (self._vrnt_phypos, self._vrnt_chrgrp)   # loci default keys
                emess = "vrnt_phypos, vrnt_chrgrp are None"     # loci error message

            keys = tuple(k for k in keys if k is not None)      # remove None keys
            if len(keys) == 0:                                  # raise error if needed
                raise RuntimeError("cannot lexsort on axis {0}: {1}".format(axis, emess))
        else:
            l = self._mat.shape[axis]
            for i,k in enumerate(keys):
                if len(k) != l:
                    raise ValueError("cannot lexsort on axis {0}: key {1} is incompatible with axis length {2}".format(axis,i,l))

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def lexsort_taxa(self, keys = None, **kwargs):
        """
        Perform an indirect stable sort using a sequence of keys along the taxa
        axis.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        return self.lexsort(
            keys = keys,
            axis = 1,
            **kwargs
        )

    def lexsort_vrnt(self, keys = None, **kwargs):
        """
        Perform an indirect stable sort using a sequence of keys along the
        variant axis.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        return self.lexsort(
            keys = keys,
            axis = 2,
            **kwargs
        )

    def reorder(self, indices, axis = -1, **kwargs):
        """
        Reorder the VariantMatrix.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        axis : int
            The axis over which to reorder values.

        """
        axis = get_axis(axis, self._mat.ndim)                   # transform axis number to an index

        ########################################################################
        if axis == 0:                                           ### PHASE AXIS
            self._mat = self._mat[indices,:,:]                  # reorder geno array
        ########################################################################
        elif axis == 1:                                         ### TAXA AXIS
            self._mat = self._mat[:,indices,:]                  # reorder geno array
            if self._taxa is not None:
                self._taxa = self._taxa[indices]                # reorder taxa array
            if self._taxa_grp is not None:
                self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array
        ########################################################################
        elif axis == 2:                                         ### LOCUS AXIS
            self._mat = self._mat[:,:,indices]                  # reorder geno array
            self._vrnt_chrgrp = self._vrnt_chrgrp[indices]      # reorder chromosome group array
            self._vrnt_phypos = self._vrnt_phypos[indices]      # reorder physical position array
            if self._vrnt_name is not None:
                self._vrnt_name = self._vrnt_name[indices]      # reorder marker name array
            if self._vrnt_genpos is not None:
                self._vrnt_genpos = self._vrnt_genpos[indices]  # reorder map position array
            if self._vrnt_xoprob is not None:
                self._vrnt_xoprob = self._vrnt_xoprob[indices]  # reorder crossover probability array
            if self._vrnt_hapgrp is not None:
                self._vrnt_hapgrp = self._vrnt_hapgrp[indices]  # reorder haplotype group array
            if self._vrnt_mask is not None:
                self._vrnt_mask = self._vrnt_mask[indices]      # reorder variant mask array

    def reorder_taxa(self, indices, **kwargs):
        """
        Reorder elements of the Matrix along the taxa axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.reorder(
            indices = indices,
            axis = 1,
            **kwargs
        )

    def reorder_vrnt(self, indices, **kwargs):
        """
        Reorder elements of the Matrix along the variant axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        self.reorder(
            indices = indices,
            axis = 2,
            **kwargs
        )

    def sort(self, keys = None, axis = -1, **kwargs):
        """
        Reset metadata for corresponding axis: name, stix, spix, len.
        Sort the VariantMatrix using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using vrnt_chrgrp as primary key, and
            vrnt_phypos as secondary key.
        axis : int
            The axis over which to sort values.
        """
        axis = get_axis(axis, self._mat.ndim)                   # transform axis number to an index

        if axis == 0:
            pass
        elif axis == 1:
            # reset taxa group metadata
            self.taxa_grp_name = None
            self.taxa_grp_stix = None
            self.taxa_grp_spix = None
            self.taxa_grp_len = None
        elif axis == 2:
            # reset variant group metadata
            self.vrnt_chrgrp_name = None
            self.vrnt_chrgrp_stix = None
            self.vrnt_chrgrp_spix = None
            self.vrnt_chrgrp_len = None

        # get indices for sort
        indices = self.lexsort(keys, axis)

        # reorder internals
        self.reorder(indices, axis)

    def sort_taxa(self, keys = None, **kwargs):
        """
        Sort slements of the Matrix along the taxa axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.
        """
        self.sort(
            keys = keys,
            axis = 1,
            **kwargs
        )

    def sort_vrnt(self, keys = None, **kwargs):
        """
        Sort slements of the Matrix along the variant axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.
        """
        self.sort(
            keys = keys,
            axis = 2,
            **kwargs
        )

    ################### Grouping Methods ###################
    def group(self, axis = -1, **kwargs):
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        axis = get_axis(axis, self._mat.ndim)                   # transform axis number to an index

        # sort along taxa axis
        self.sort(axis = axis)

        if axis == 0:
            pass    # no phase grouping protocols
        elif axis == 1:
            if self._taxa_grp is not None:
                # get unique taxa group names, starting indices, group lengths
                uniq = numpy.unique(self._taxa_grp, return_index = True, return_counts = True)
                # make assignments to instance data
                self._taxa_grp_name, self._taxa_grp_stix, self._taxa_grp_len = uniq
                # calculate stop indices
                self._taxa_grp_spix = self._taxa_grp_stix + self._taxa_grp_len
        elif axis == 2:
            # get unique chromosome group names, starting indices, group lengths
            uniq = numpy.unique(self._vrnt_chrgrp, return_index = True, return_counts = True)
            # make assignments to instance data
            self._vrnt_chrgrp_name, self._vrnt_chrgrp_stix, self._vrnt_chrgrp_len = uniq
            # calculate stop indices
            self._vrnt_chrgrp_spix = self._vrnt_chrgrp_stix + self._vrnt_chrgrp_len

    def group_taxa(self, **kwargs):
        """
        Sort the Matrix along the taxa axis, then populate grouping indices for
        the taxa axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.
        """
        self.group(
            axis = 1,
            **kwargs
        )

    def group_vrnt(self, **kwargs):
        """
        Sort the Matrix along the variant axis, then populate grouping indices
        for the variant axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.
        """
        self.group(
            axis = 2,
            **kwargs
        )

    def is_grouped(self, axis = -1, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        axis = get_axis(axis, self._mat.ndim)                   # transform axis number to an index

        if axis == 0:
            pass    # not grouped because it cannot be grouped
        elif axis == 1:
            return (
                (self._taxa_grp_name is not None) and
                (self._taxa_grp_stix is not None) and
                (self._taxa_grp_spix is not None) and
                (self._taxa_grp_len is not None)
            )
        elif axis == 2:
            return (
                (self._vrnt_chrgrp_name is not None) and
                (self._vrnt_chrgrp_stix is not None) and
                (self._vrnt_chrgrp_spix is not None) and
                (self._vrnt_chrgrp_len is not None)
            )
        return False

    def is_grouped_taxa(self, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped along the taxa
        axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        return self.is_grouped(
            axis = 1,
            **kwargs
        )

    def is_grouped_vrnt(self, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped along the
        variant axis.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        return self.is_grouped(
            axis = 2,
            **kwargs
        )

    ################# Interpolation Methods ################
    def interp_genpos(self, gmap, **kwargs):
        """
        Interpolate genetic map postions for variants using a GeneticMap

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        """
        # check if gmap is a GeneticMap
        check_is_GeneticMap(gmap, "gmap")

        # interpolate postions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

    def interp_xoprob(self, gmap, gmapfn, **kwargs):
        """
        Interpolate genetic map positions AND crossover probabilities between
        sequential markers using a GeneticMap and a GeneticMapFunction.

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        gmapfn : GeneticMapFunction
            A genetic map function from which to interpolate crossover
            probabilities for loci within the VariantMatrix.
        """
        # check data types
        check_is_GeneticMap(gmap, "gmap")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")

        # check if self has been sorted and grouped
        if not self.is_grouped():
            raise RuntimeError("must be grouped first before interpolation of crossover probabilities")

        # interpolate genetic positions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

        # interpolate crossover probabilities
        self.vrnt_xoprob = gmapfn.rprob1g(gmap, self._vrnt_chrgrp, self._vrnt_genpos)

    ################## Clustering Methods ##################
    def assign_hapgrp(self, k, **kwargs):
        """
        Assign haplotype groups using k-means clustering.

        Parameters
        ----------
        k : int, numpy.ndarray
            Number of haplotype groups to assign to each
        **kwargs
            Additional keyword arguments.
        """
        # TODO: implement me
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def from_vcf(fname):
        """
        Does not ensure that data is phased, just reads it as phased.
        """
        # make VCF iterator
        vcf = cyvcf2.VCF(fname)

        # extract taxa names from vcf header
        taxa = numpy.object_(vcf.samples)

        # make empty lists to store extracted values
        mat = []
        vrnt_chrgrp = []
        vrnt_phypos = []
        vrnt_name = []

        # iterate through VCF file and accumulate variants
        for variant in vcf:
            # append chromosome integer
            vrnt_chrgrp.append(int(variant.CHROM))

            # append variant position coordinates
            vrnt_phypos.append(variant.POS)

            # append marker name
            vrnt_name.append(str(variant.ID))

            # extract allele states + whether they are phased or not
            phases = numpy.int8(variant.genotypes)

            # check that they are all phased
            #pybropt.util.check_matrix_all_value(phases[:,2], "is_phased", True)

            # TODO: maybe modify shapes here to avoid transpose and copy below?
            # append genotype states
            mat.append(phases[:,0:2].copy())

        # convert and transpose genotype matrix
        mat = numpy.int8(mat).transpose(2,1,0) # may want to copy()?

        # convert to numpy.ndarray
        vrnt_chrgrp = numpy.int64(vrnt_chrgrp)  # convert to int64 array
        vrnt_phypos = numpy.int64(vrnt_phypos)  # convert to int64 array
        vrnt_name = numpy.object_(vrnt_name)    # convert to object array

        pvm = DensePhasedGenotypeVariantMatrix(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            taxa = taxa
        )

        return pvm



################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhasedGenotypeVariantMatrix(v):
    return isinstance(v, DensePhasedGenotypeVariantMatrix)

def check_is_DensePhasedGenotypeVariantMatrix(v, varname):
    if not isinstance(v, DensePhasedGenotypeVariantMatrix):
        raise TypeError("'%s' must be a DensePhasedGenotypeVariantMatrix." % varname)

def cond_check_is_DensePhasedGenotypeVariantMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DensePhasedGenotypeVariantMatrix(v, varname)
