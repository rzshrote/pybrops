class GeneticMap:
    """docstring for GeneticMap."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class GeneticMap.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments.
        """
        super(GeneticMap, self).__init__()

    def __len__(self):
        """Get the number of markers in the genetic map."""
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################### Data Properites ####################
    def vrnt_chrgrp():
        doc = "Variant chromosome group label property."
        def fget(self):
            """Get variant chromosome group lable array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group lable array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group lable array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp = property(**vrnt_chrgrp())

    def vrnt_phypos():
        doc = "Variant physical position property."
        def fget(self):
            """Get variant physical position array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant physical position array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant physical position array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_phypos = property(**vrnt_phypos())

    def vrnt_genpos():
        doc = "Variant genetic position property."
        def fget(self):
            """Get variant genetic position array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant genetic position array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant genetic position array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    ################# Metadata Properites ##################
    def vrnt_chrgrp_name():
        doc = "Variant chromosome group names property."
        def fget(self):
            """Get variant chromosome group name array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group name array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group name array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_name = property(**vrnt_chrgrp_name())

    def vrnt_chrgrp_stix():
        doc = "Variant chromosome group start indices property."
        def fget(self):
            """Get variant chromosome group start indices array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group start indices array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group start indices array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_stix = property(**vrnt_chrgrp_stix())

    def vrnt_chrgrp_spix():
        doc = "Variant chromosome group stop indices property."
        def fget(self):
            """Get variant chromosome group stop indices array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group stop indices array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group stop indices array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_spix = property(**vrnt_chrgrp_spix())

    def vrnt_chrgrp_len():
        doc = "Variant chromosome group length property."
        def fget(self):
            """Get variant chromosome group length array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set variant chromosome group length array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete variant chromosome group length array"""
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_len = property(**vrnt_chrgrp_len())

    ################## Spline Properites ###################
    def spline():
        doc = "Interpolation spline property."
        def fget(self):
            """Get interpolation spline(s)"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set interpolation spline(s)"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete interpolation spline(s)"""
            raise NotImplementedError("method is abstract")
        return locals()
    spline = property(**spline())

    ############# Spline Metadata Properites ###############
    def spline_kind():
        doc = "Spline kind property."
        def fget(self):
            """Get the spline kind"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the spline kind"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the spline kind"""
            raise NotImplementedError("method is abstract")
        return locals()
    spline_kind = property(**spline_kind())

    def spline_fill_value():
        doc = "Spline fill value property."
        def fget(self):
            """Get the spline fill value"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the spline fill value"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the spline fill value"""
            raise NotImplementedError("method is abstract")
        return locals()
    spline_fill_value = property(**spline_fill_value())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys, **kwargs):
        """
        Perform an indirect stable sort using a sequence of keys.

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
        raise NotImplementedError("method is abstract")

    def reorder(self, indices, **kwargs):
        """
        Reorder markers in the GeneticMap using an array of indices.
        Note this modifies the GeneticMap in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def sort(self, keys, **kwargs):
        """
        Sort slements of the GeneticMap using a sequence of keys.
        Note this modifies the GeneticMap in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Grouping Methods ###################
    def group(self, **kwargs):
        """
        Sort the GeneticMap, then populate grouping indices.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def is_grouped(self, **kwargs):
        """
        Determine whether the GeneticMap has been sorted and grouped.

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
        raise NotImplementedError("method is abstract")

    ################ Insert/Delete Methods #################
    def remove(self, indices, **kwargs):
        """
        Remove indices from the GeneticMap. Sort and group internal arrays.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
        **kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def select(self, indices, **kwargs):
        """
        Keep only selected markers, removing all others from the GeneticMap.
        Sort and group internal arrays.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
        **kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def prune(self, nt, M, **kwargs):
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
        **kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    # TODO: def insert(self, ...)

    ################## Integrity Methods ###################
    def congruence(self):
        """
        If not grouped, will group.
        Assess physical and genetic map site congruency.

        Notes:
            This assumes high contiguity between physical and genetic maps
            (i.e. a high quality reference genome). This assumption may cause
            major issues if there are incorrect markers at the beginning of the
            chromosome.
        Returns
        -------
        concordancy : numpy.ndarray
            A boolean matrix of map concordancies where:
            True = the current marker has a map_pos >= the previous position
            False = the current marker has a map_pos < the previous position
        """
        raise NotImplementedError("method is abstract")

    def is_congruent(self):
        """
        Determine if the genetic map is congruent
        Determine if all sites in the genetic map demonstrate congruence with
        their supposed physical and genetic positions.
        """
        raise NotImplementedError("method is abstract")

    def remove_discrepancies(self):
        """
        Remove discrepancies between the physical map and the genetic map.
        In instances of conflict, assume that the physical map is correct.

        Note:
            This assumption may cause major issues if there are incorrect
            markers at the beginning of the chromosome.
        """
        raise NotImplementedError("method is abstract")

    ################# Interpolation Methods ################
    def build_spline(self, kind, fill_value, **kwargs):
        """
        Build a spline for estimating genetic map distances. This is built
        using the marker start indices (self.chr_start)

        Parameters
        ----------
        kind : str
            Specifies the kind of interpolation as a string.
        fill_value : obj
            Fill value for points extrapolated outside the spline.
        **kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def has_spline(self):
        """Return whether or not the GeneticMap has a built spline."""
        raise NotImplementedError("method is abstract")

    def interp_genpos(self, vrnt_chrgrp, vrnt_phypos):
        """
        Interpolate genetic positions given variant physical positions

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
        vrnt_phypos : numpy.ndarray

        Returns
        -------
        out : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")

    def interp_gmap(self, vrnt_chrgrp, vrnt_phypos, **kwargs):
        """
        Interpolate a new genetic map from the current genetic map.
        Associate spline of current GeneticMap with new GeneticMap.
        """
        raise NotImplementedError("method is abstract")

    ############### Genetic Distance Methods ###############
    def gdist1g(self, vrnt_chrgrp, vrnt_genpos, ast, asp):
        """
        Calculate sequential genetic distances using genetic positions.
        Requires vrnt_chrgrp and vrnt_genpos to have been sorted descending.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
        vrnt_genpos : numpy.ndarray
            A 1D array of variant genetic positions.
        ast : int, None
            Optional array start index (inclusive).
        asp : int, None
            Optional array stop index (exclusive).

        Returns
        -------
        gdist : numpy.ndarray
            A 1D array of distances between the marker prior.
        """
        raise NotImplementedError("method is abstract")

    def gdist2g(self, vrnt_chrgrp, vrnt_genpos, rst, rsp, cst, csp):
        """
        Calculate pairwise genetic distances using genetic positions.
        Requires vrnt_chrgrp and vrnt_genpos to have been sorted descending.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
        vrnt_genpos : numpy.ndarray
            A 1D array of variant genetic positions.
        rst : int, None
            Optional row start index (inclusive).
        rsp : int, None
            Optional row stop index (exclusive).
        cst : int, None
            Optional column start index (inclusive).
        csp : int, None
            Optional column stop index (exclusive).

        Returns
        -------
        gdist : numpy.ndarray
            A 2D array of pairwise distances between markers.
        """
        raise NotImplementedError("method is abstract")

    def gdist1p(self, vrnt_chrgrp, vrnt_phypos, ast, asp):
        """
        Calculate sequential genetic distances using physical positions.
        Requires vrnt_chrgrp and vrnt_phypos to have been sorted descending.
        Requires a spline to have been built beforehand.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
        vrnt_phypos : numpy.ndarray
            A 1D array of variant physical positions.
        ast : int, None
            Optional array start index (inclusive).
        asp : int, None
            Optional array stop index (exclusive).

        Returns
        -------
        gdist : numpy.ndarray
            A 1D array of distances between the marker prior.
        """
        raise NotImplementedError("method is abstract")

    def gdist2p(self, vrnt_chrgrp, vrnt_phypos, rst, rsp, cst, csp):
        """
        Calculate pairwise genetic distances using physical positions.
        Requires vrnt_chrgrp and vrnt_phypos to have been sorted descending.
        Requires a spline to have been built beforehand.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
        vrnt_phypos : numpy.ndarray
            A 1D array of variant physical positions.
        rst : int, None
            Optional row start index (inclusive).
        rsp : int, None
            Optional row stop index (exclusive).
        cst : int, None
            Optional column start index (inclusive).
        csp : int, None
            Optional column stop index (exclusive).

        Returns
        -------
        gdist : numpy.ndarray
            A 2D array of pairwise distances between markers.
        """
        raise NotImplementedError("method is abstract")

    #################### Export Methods ####################
    def to_pandas_df(self):
        """
        Convert a GeneticMap object to a pandas DataFrame.

        Returns
        -------
        df : pandas.DataFrame
            A pandas DataFrame containing genetic map data.
        """
        raise NotImplementedError("method is abstract")

    def to_csv(self, fname, sep, header, index, **kwargs):
        """
        Convert a GeneticMap object to a csv file.

        Parameters
        ----------
        fname : str
        sep : str
        header : bool
        index : bool, int
        **kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GeneticMap(v):
    """
    Determine whether an object is a GeneticMap.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GeneticMap object instance.
    """
    return isinstance(v, GeneticMap)

def check_is_GeneticMap(v, vname):
    """
    Check if object is of type GeneticMap. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GeneticMap):
        raise TypeError("variable '{0}' must be a GeneticMap".format(vname))

def cond_check_is_GeneticMap(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type GeneticMap. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        GeneticMap.
    """
    if cond(v):
        check_is_GeneticMap(v, vname)
