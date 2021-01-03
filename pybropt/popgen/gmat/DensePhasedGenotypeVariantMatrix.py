# import 3rd party modules we'll need
import cyvcf2
import numpy

# import our libraries
from . import GenotypeVariantMatrix
from . import DensePhasedGenotypeMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.error import check_ndarray_dtype
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_ndim
from pybropt.core.error import cond_check_ndarray_dtype
from pybropt.core.error import cond_check_ndarray_axis_len
from pybropt.core.error import check_ndarray_axis_len


class DensePhasedGenotypeVariantMatrix(DensePhasedGenotypeMatrix,GenotypeVariantMatrix):
    """docstring for PhasedVariantMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, vrnt_chrgrp, vrnt_phypos, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_mask = None, taxa = None, taxa_grp = None, **kwargs):
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
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_mask = vrnt_mask,
            taxa = taxa,
            taxa_grp = taxa_grp
        )

        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_name = vrnt_name
        self.vrnt_genpos = vrnt_genpos
        self.vrnt_xoprob = vrnt_xoprob
        self.vrnt_hapgrp = vrnt_hapgrp
        self.vrnt_mask = vrnt_mask
        self.taxa = taxa
        self.taxa_grp = taxa_grp

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

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
            cond_check_ndarray_dtype(value, "vrnt_name", numpy.string_)
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

    ################# Taxa Data Properites #################
    def taxa():
        doc = "The taxa property."
        def fget(self):
            return self._taxa
        def fset(self, value):
            cond_check_is_ndarray(value, "taxa")
            cond_check_ndarray_dtype(value, "taxa", numpy.string_)
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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def axis_index(self, axis):
        """
        Return an index (unsigned) from a provided axis integer (signed)

        Parameters
        ----------
        axis : int
            Integer representation of the axis. Can be in range (-ndim,ndim).
            If outside this range, will raise an AxisError.

        Returns
        -------
        index : int
            Index representation of the axis. In range [0,ndim).
        """
        # get the number of axis
        naxes = self._mat.ndim

        # handle axis argument
        if (axis >= naxes) or (axis < -naxes):
            raise IndexError("axis {0} is out of bounds for array of dimension {1}".format(axis, naxes))

        # modulo the axis number to get the axis (in the case of negative axis)
        axis %= naxes

        return axis

    def lexsort(self, keys = None, axis = -1):
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
        axis = self.axis_index(axis)                            # transform axis number to an index
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
                    raise RuntimeError("cannot lexsort on axis %s: key %s is incompatible with axis length %s" % (axis, i, l))

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

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

    def sort(self, keys = None, axis = -1):
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
        # get axis
        axis = self.axis_index(axis)

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

    ################### Grouping Methods ###################
    def group(self, axis = -1):
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        # get axis index
        axis = self.axis_index(axis)

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

    def is_grouped(self, axis = -1):
        """
        Determine whether the Matrix has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        # convert axis to index
        axis = self.axis_index(axis)

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

    ################# Interpolation Methods ################
    def interp_genpos(self, gmap):
        """
        Interpolate genetic map postions for variants using a GeneticMap

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        """
        # check if gmap is a GeneticMap
        pybropt.popgen.gmap.check_is_GeneticMap(gmap)

        # interpolate postions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

    def interp_xoprob(self, gmap, gmapfn):
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
        pybropt.popgen.gmap.check_is_GeneticMap(gmap)
        pybropt.popgen.gmap.check_is_GeneticMapFunction(gmapfn)

        # check if self has been sorted and grouped
        if not self.is_grouped():
            raise RuntimeError("must be grouped first before interpolation of crossover probabilities")

        # interpolate genetic positions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

        # interpolate crossover probabilities
        self.vrnt_xoprob = gmapfn.gdist1g(gmap, self._vrnt_chrgrp, self._vrnt_genpos)

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
        taxa = numpy.string_(vcf.samples)

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
            vrnt_name.append(variant.ID)

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
        vrnt_name = numpy.string_(vrnt_name)    # convert to string array

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
