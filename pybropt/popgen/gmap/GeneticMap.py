class GeneticMap:
    """docstring for GeneticMap."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GeneticMap, self).__init__()

    def __len__(self):
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################### Data Properites ####################
    def vrnt_chrgrp():
        doc = "The vrnt_chrgrp property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp = property(**vrnt_chrgrp())

    def vrnt_phypos():
        doc = "The vrnt_phypos property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_phypos = property(**vrnt_phypos())

    def vrnt_genpos():
        doc = "The vrnt_genpos property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    ################# Metadata Properites ##################
    def vrnt_chrgrp_name():
        doc = "The vrnt_chrgrp_name property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_name = property(**vrnt_chrgrp_name())

    def vrnt_chrgrp_stix():
        doc = "The vrnt_chrgrp_stix property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_stix = property(**vrnt_chrgrp_stix())

    def vrnt_chrgrp_spix():
        doc = "The vrnt_chrgrp_spix property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_spix = property(**vrnt_chrgrp_spix())

    def vrnt_chrgrp_len():
        doc = "The vrnt_chrgrp_len property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vrnt_chrgrp_len = property(**vrnt_chrgrp_len())

    ################## Spline Properites ###################
    def spline():
        doc = "The spline property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    spline = property(**spline())

    ############# Spline Metadata Properites ###############
    def spline_kind():
        doc = "The spline_kind property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    spline_kind = property(**spline_kind())

    def spline_fill_value():
        doc = "The spline_fill_value property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    spline_fill_value = property(**spline_fill_value())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys):
        raise NotImplementedError("method is abstract")

    def reorder(self, indices):
        raise NotImplementedError("method is abstract")

    def sort(self, keys):
        """
        Set variant chromosome group name, stix, spix, len to None.
        Sort according to keys.
        Preserves spline if it exists.
        """
        raise NotImplementedError("method is abstract")

    def group(self):
        """
        Sort genetic map, then populate grouping indices.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        raise NotImplementedError("method is abstract")

    def is_grouped(self):
        """
        Determine whether the GeneticMap has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        raise NotImplementedError("method is abstract")

    ################ Insert/Delete Methods #################
    def remove(self, indices):
        raise NotImplementedError("method is abstract")

    def select(self, indices):
        raise NotImplementedError("method is abstract")

    def prune(self, nt, M):
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
    def build_spline(self, kind, fill_value):
        raise NotImplementedError("method is abstract")

    def has_spline(self):
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
        raise NotImplementedError("method is abstract")

    def to_csv(self, fname, sep, header, index):
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GeneticMap(v):
    return isinstance(v, GeneticMap)

def check_is_GeneticMap(v, vname):
    if not isinstance(v, GeneticMap):
        raise TypeError("variable '{0}' must be a GeneticMap".format(vname))

def cond_check_is_GeneticMap(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GeneticMap(v, vname)
