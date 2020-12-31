import numpy
import pandas
from scipy.interpolate import interp1d

from pybropt.core.error import check_is_ndarray

from . import GeneticMap

class ExtendedGeneticMap(GeneticMap):
    """docstring for ExtendedGeneticMap."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, vrnt_chrgrp, vrnt_phypos, vrnt_stop, vrnt_genpos, vrnt_name = None, vrnt_fncode = None, **kwargs):
        super(ExtendedGeneticMap, self).__init__(**kwargs)
        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_stop = vrnt_stop
        self.vrnt_genpos = vrnt_genpos
        self.vrnt_name = vrnt_name
        self.vrnt_fncode = vrnt_fncode
        # TODO: check all lengths equivalent

        # sort and group
        self.group()

    def __len__(self):
        return len(self._vrnt_genpos)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################### Data Properites ####################
    def vrnt_chrgrp():
        doc = "The vrnt_chrgrp property."
        def fget(self):
            return self._vrnt_chrgrp
        def fset(self, value):
            check_is_ndarray(value, "vrnt_chrgrp")
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
            check_is_ndarray(value, "vrnt_phypos")
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
            check_is_ndarray(value, "vrnt_genpos")
            pybropt.util.check_matrix_dtype(value, "vrnt_genpos", numpy.float64)
            pybropt.util.check_matrix_ndim(value, "vrnt_genpos", 1)
            self._vrnt_genpos = value
        def fdel(self):
            del self._vrnt_genpos
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    def vrnt_stop():
        doc = "The vrnt_stop property."
        def fget(self):
            return self._vrnt_stop
        def fset(self, value):
            check_is_ndarray(value, "vrnt_stop")
            pybropt.util.check_matrix_dtype(value, "vrnt_stop", numpy.int64)
            pybropt.util.check_matrix_ndim(value, "vrnt_stop", 1)
            self._vrnt_stop = value
        def fdel(self):
            del self._vrnt_stop
        return locals()
    vrnt_stop = property(**vrnt_stop())

    def vrnt_name():
        doc = "The vrnt_name property."
        def fget(self):
            return self._vrnt_name
        def fset(self, value):
            pybropt.util.cond_check_is_matrix(value, "vrnt_name")
            pybropt.util.cond_check_matrix_dtype(value, "vrnt_name", numpy.string_)
            pybropt.util.cond_check_matrix_ndim(value, "vrnt_name", 1)
            self._vrnt_name = value
        def fdel(self):
            del self._vrnt_name
        return locals()
    vrnt_name = property(**vrnt_name())

    def vrnt_fncode():
        doc = "The vrnt_fncode property."
        def fget(self):
            return self._vrnt_fncode
        def fset(self, value):
            pybropt.util.cond_check_is_matrix(value, "vrnt_fncode")
            pybropt.util.cond_check_matrix_dtype(value, "vrnt_fncode", numpy.string_)
            pybropt.util.cond_check_matrix_ndim(value, "vrnt_fncode", 1)
            self._vrnt_fncode = value
        def fdel(self):
            del self._vrnt_fncode
        return locals()
    vrnt_fncode = property(**vrnt_fncode())

    ################# Metadata Properites ##################
    def vrnt_chrgrp_name():
        doc = "The vrnt_chrgrp_name property."
        def fget(self):
            return self._vrnt_chrgrp_name
        def fset(self, value):
            pybropt.util.cond_check_is_matrix(value, "vrnt_chrgrp_name")
            pybropt.util.cond_check_matrix_dtype(value, "vrnt_chrgrp_name", numpy.int64)
            pybropt.util.cond_check_matrix_ndim(value, "vrnt_chrgrp_name", 1)
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
            pybropt.util.cond_check_is_matrix(value, "vrnt_chrgrp_stix")
            pybropt.util.cond_check_matrix_dtype(value, "vrnt_chrgrp_stix", numpy.int64)
            pybropt.util.cond_check_matrix_ndim(value, "vrnt_chrgrp_stix", 1)
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
            pybropt.util.cond_check_is_matrix(value, "vrnt_chrgrp_spix")
            pybropt.util.cond_check_matrix_dtype(value, "vrnt_chrgrp_spix", numpy.int64)
            pybropt.util.cond_check_matrix_ndim(value, "vrnt_chrgrp_spix", 1)
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
            pybropt.util.cond_check_is_matrix(value, "vrnt_chrgrp_len")
            pybropt.util.cond_check_matrix_dtype(value, "vrnt_chrgrp_len", numpy.int64)
            pybropt.util.cond_check_matrix_ndim(value, "vrnt_chrgrp_len", 1)
            self._vrnt_chrgrp_len = value
        def fdel(self):
            del self._vrnt_chrgrp_len
        return locals()
    vrnt_chrgrp_len = property(**vrnt_chrgrp_len())

    ################## Spline Properites ###################
    def spline():
        doc = "The spline property."
        def fget(self):
            return self._spline
        def fset(self, value):
            pybropt.util.cond_check_is_dict(value, "spline")
            self._spline = value
        def fdel(self):
            del self._spline
        return locals()
    spline = property(**spline())

    ############# Spline Metadata Properites ###############
    def spline_kind():
        doc = "The spline_kind property."
        def fget(self):
            return self._spline_kind
        def fset(self, value):
            pybropt.util.cond_check_is_string(value, "spline_kind")
            self._spline_kind = value
        def fdel(self):
            del self._spline_kind
        return locals()
    spline_kind = property(**spline_kind())

    def spline_fill_value():
        doc = "The spline_fill_value property."
        def fget(self):
            return self._spline_fill_value
        def fset(self, value):
            self._spline_fill_value = value
        def fdel(self):
            del self._spline_fill_value
        return locals()
    spline_fill_value = property(**spline_fill_value())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys = None):
        """
        Generate indices for lexsort.
        """
        # if no keys were provided, set a default
        if keys is None:
            keys = (
                self._vrnt_genpos,  # 3rd priority
                self._vrnt_phypos,  # 2nd priority
                self._vrnt_chrgrp   # 1st priority
            )
        else:
            # check that matrix lengths are the same
            for i,k in enumerate(keys):
                pybropt.util.check_matrix_size(k, "key"+i, self.__len__())

        # build tuple
        keys = tuple(k for k in keys if k is not None)

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def reorder(self, indices):
        """
        Reorder the genetic map.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        """
        # sort internal self
        self._vrnt_chrgrp = self._vrnt_chrgrp[indices]
        self._vrnt_phypos = self._vrnt_phypos[indices]
        self._vrnt_stop = self._vrnt_stop[indices]
        self._vrnt_genpos = self._vrnt_genpos[indices]
        if self._vrnt_name is not None:
            self._vrnt_name = self._vrnt_name[indices]
        if self._vrnt_fncode is not None:
            self._vrnt_fncode = self._vrnt_fncode[indices]

    def sort(self, keys = None):
        """
        Set variant chromosome group name, stix, spix, len to None.
        Sort according to keys.
        Preserves spline if it exists.
        """
        # reset variant chromosome group name, stix, spix, len
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len = None

        # get indices for sort
        indices = self.lexsort(keys)

        # reorder internals
        self.reorder(indices)

    def group(self):
        """
        Sort genetic map, then populate grouping indices.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        # sort the genetic map using default keys
        self.sort()

        # get unique group names, start indices, group lengths
        uniq = numpy.unique(self._vrnt_chrgrp, return_index = True, return_counts = True)

        # make assignments
        self._vrnt_chrgrp_name, self._vrnt_chrgrp_stix, self._vrnt_chrgrp_len = uniq

        # calculate chr_grp_spix (stop indices)
        self._vrnt_chrgrp_spix = self._vrnt_chrgrp_stix + self._vrnt_chrgrp_len

    def is_grouped(self):
        """
        Determine whether the GeneticMap has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        return (
            (self._vrnt_chrgrp_name is not None) and
            (self._vrnt_chrgrp_stix is not None) and
            (self._vrnt_chrgrp_spix is not None) and
            (self._vrnt_chrgrp_len is not None)
        )

    ################ Insert/Delete Methods #################
    def remove(self, indices):
        """
        Remove indices from the MarkerSet.
        Sort and group internal arrays.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
        """
        # delete indices from self
        self._vrnt_chrgrp = numpy.delete(self._vrnt_chrgrp, indices)
        self._vrnt_phypos = numpy.delete(self._vrnt_phypos, indices)
        self._vrnt_stop = numpy.delete(self._vrnt_stop, indices)
        self._vrnt_genpos = numpy.delete(self._vrnt_genpos, indices)
        if self._vrnt_name is not None:
            self._vrnt_name = numpy.delete(self._vrnt_name, indices)
        if self._vrnt_fncode is not None:
            self._vrnt_fncode = numpy.delete(self._vrnt_fncode, indices)

        # sort and group
        self.group()

    def select(self, indices):
        """
        Keep only selected markers, removing all others from the MarkerSet.
        Sort and group internal arrays.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
        """
        # keep only selected markers.
        self._vrnt_chrgrp = self._vrnt_chrgrp[indices]
        self._vrnt_phypos = self._vrnt_phypos[indices]
        self._vrnt_stop = self._vrnt_stop[indices]
        self._vrnt_genpos = self._vrnt_genpos[indices]
        if self._vrnt_name is not None:
            self._vrnt_name = self._vrnt_name[indices]
        if self._vrnt_fncode is not None:
            self._vrnt_fncode = self._vrnt_fncode[indices]

        # sort and group
        self.group()

    def prune(self, nt = None, M = None):
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
        # check if we have acceptible inputs
        if (nt is None) and (M is None):
            raise ValueError("'nt' and 'M' cannot both be None")

        # if not grouped, sort and group
        if not is_grouped():
            self.group()

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

        # select the desired marker indices
        self.select(indices)

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
            Assumes the first marker on the chromosome is placed correctly.
        Returns
        -------
        out : numpy.ndarray
            A boolean matrix of map concordancies where:
            True = the current marker has a map_pos >= the previous position
            False = the current marker has a map_pos < the previous position
        """
        # group as a prerequisite
        if not self.is_grouped():
            self.group()

        # make empty bool vector
        out = numpy.zeros(len(self._vrnt_phypos), dtype = 'bool')

        # for each linkage group
        for st,sp in zip(self._vrnt_chrgrp_stix, self._vrnt_chrgrp_spix):
            # we assume the first element is correctly placed
            out[st] = True
            # test if the current element is <= the previous element
            out[st+1:sp] = self._vrnt_genpos[st:sp-1] <= self._vrnt_genpos[st+1:sp]

        # return results
        return out

    def is_congruent(self):
        """
        Determine if the genetic map is congruent
        Determine if all sites in the genetic map demonstrate congruence with
        their supposed physical and genetic positions.
        """
        # determine if all sites are congruent
        out = numpy.all(self.congruence())
        return out

    def remove_discrepancies(self):
        """
        Remove discrepancies between the physical map and the genetic map.
        In instances of conflict, assume that the physical map is correct.

        Note:
            This assumption may cause major issues if there are incorrect
            markers at the beginning of the chromosome.
        """
        # get boolean mask for concordant map positions
        mask = self.congruence()

        # if all are concordant, exit this function
        if( numpy.all(concordancy) ):
            return

        # select (and group) only sites that are congruent
        self.select(mask)

    ################ Interpolation Methods #################
    # TODO: do not require sorted internals to build spline
    def build_spline(self, kind = 'linear', fill_value = 'extrapolate'):
        """
        Build a spline for estimating genetic map distances. This is built
        using the marker start indices (self.chr_start)

        Parameters
        ----------
        kind : str, default = 'linear'
            Specifies the kind of interpolation as a string ('linear',
            'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous',
            'next', where 'zero', 'slinear', 'quadratic' and 'cubic' refer to a
            spline interpolation of zeroth, first, second or third order;
            'previous' and 'next' simply return the previous or next value of
            the point) or as an integer specifying the order of the spline
            interpolator to use.
        fill_value : array-like, {'extrapolate'}, default = 'extrapolate'
            If 'extrapolate', then points outside the data range will be
            extrapolated.
            If a ndarray (or float), this value will be used to fill in for
            requested points outside of the data range. If not provided, then
            the default is NaN. The array-like must broadcast properly to the
            dimensions of the non-interpolation axes.
            If a two-element tuple, then the first element is used as a fill
            value for x_new < x[0] and the second element is used for
            x_new > x[-1]. Anything that is not a 2-element tuple (e.g., list
            or ndarray, regardless of shape) is taken to be a single array-like
            argument meant to be used for both bounds as below,
            above = fill_value, fill_value.
        """
        # get unique chromosome group labels
        uniq = numpy.unique(self._vrnt_chrgrp)

        # make an empty dictionary for map splines
        self._spline = {}

        # set the spline type
        self._spline_kind = kind
        self._spline_fill_value = fill_value

        # iterate through unique groups and build spline
        for grp in uniq:
            # get mask of elements belonging to grp
            mask = (self._vrnt_chrgrp == grp)

            # build spline and assign to dictionary
            self._spline[grp] = interp1d(
                x = self._vrnt_phypos[mask],    # chromosome positions
                y = self._vrnt_genpos[mask],    # map positions
                kind = kind,                    # type of interpolation
                fill_value = fill_value,        # default fill value
                assume_sorted = False           # arrays are not sorted
            )

    def has_spline(self):
        return (
            (self._spline is not None) and
            (self._spline_kind is not None) and
            (self._spline_fill_value is not None)
        )

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
        # raise error if no spline is found
        if not self.has_spline():
            raise RuntimeError("interpolation spline not built")

        # allocate empty memory
        out = numpy.empty(vrnt_phypos.shape, dtype='float64')

        # for each chromosome-position pair
        for i,(chrgrp,phypos) in enumerate(zip(vrnt_chrgrp, vrnt_phypos)):
            try:                                # try to index dict
                model = self._spline[chr_grp]   # try to get model
                out[i] = model(phypos)          # interpolate genetic position
            except KeyError:                    # if dict key not found
                out[i] = numpy.nan              # set to NaN

        return out

    def interp_gmap(self, vrnt_chrgrp, vrnt_phypos, vrnt_stop, vrnt_name = None, vrnt_fncode = None, **kwargs):
        """
        Interpolate a new genetic map from the current genetic map.
        Associate spline of current GeneticMap with new GeneticMap.
        """
        # interpolate genetic positions
        vrnt_genpos = self.interp_genpos(vrnt_chrgrp, vrnt_phypos)

        # create a new genetic map using interpolations from self
        gmap = self.__class__(
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_stop = vrnt_stop,
            vrnt_genpos = vrnt_genpos,
            vrnt_name = vrnt_name,
            vrnt_fncode = vrnt_fncode,
            **kwargs
        )

        # copy pointers to spline and spline metadata
        gmap._spline = self._spline
        gmap._spline_kind = self._spline_kind
        gmap._spline_fill_value = self._spline_fill_value

        # return new GeneticMap instance
        return gmap

    ############### Genetic Distance Methods ###############
    def gdist1g(self, vrnt_chrgrp, vrnt_genpos, ast = None, asp = None):
        """
        Calculate sequential genetic distances using genetic positions.
        Requires vrnt_chrgrp and vrnt_genpos to have been sorted descending.

        Sequential distance arrays will start every chromosome with numpy.inf!

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
        # get views of only sections we'll be looking at
        view_chrgrp = vrnt_chrgrp[ast:asp]
        view_genpos = vrnt_genpos[ast:asp]

        # get chromosome groups, start index, length of region
        uniq, start, counts = numpy.unique(view_chrgrp, return_index = True, return_counts = True)

        # get stop index
        stop = start + counts

        # allocate empty array
        gdist = numpy.empty(view_genpos.shape, dtype='float64')

        # for each chromosome group
        for st,sp in zip(start, stop):
            # infinite distance between chromosomes
            gdist[st] = numpy.inf
            # subtract distances
            gdist[st+1:sp] = view_genpos[st+1:sp] - view_genpos[st:sp-1]

        return gdist

    def gdist2g(self, vrnt_chrgrp, vrnt_genpos, rst = None, rsp = None, cst = None, csp = None):
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
            A 2D array of distances between marker pairs.
        """
        # TODO: # OPTIMIZE: fill with inf, then add distances? - assumes sort
        # create sparse meshgrid indicating where chromosomes differ
        mi, mj = numpy.meshgrid(vrnt_chrgrp[rst:rsp], vrnt_chrgrp[cst:csp], indexing='ij', sparse=True)

        # create sparse meshgrid indicating where genetic positions are
        gi, gj = numpy.meshgrid(vrnt_genpos[rst:rsp], vrnt_genpos[cst:csp], indexing='ij', sparse=True)

        # calculate absolute distances
        gdist = numpy.abs(gi - gj)

        # where chromosomes are different, set to numpy.inf
        gdist[mi != mj] = numpy.inf

        return gdist

    def gdist1p(self, vrnt_chrgrp, vrnt_phypos, ast = None, asp = None):
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

        Returns
        -------
        gdist : numpy.ndarray
            A 1D array of distances between the marker prior.
        """
        # TODO: # OPTIMIZE: don't interpolate all markers
        # interpolate genetic positions
        vrnt_genpos = self.interp_genpos(vrnt_chrgrp, vrnt_phypos)

        # calculate genetic distances
        gdist = self.gdist1g(vrnt_chrgrp, vrnt_genpos, ast, asp)

        return gdist

    def gdist2p(self, vrnt_chrgrp, vrnt_phypos, rst = None, rsp = None, cst = None, csp = None):
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
            A 2D array of distances between marker pairs.
        """
        # TODO: # OPTIMIZE: don't interpolate all markers
        # interpolate genetic positions
        vrnt_genpos = self.interp_genpos(vrnt_chrgrp, vrnt_phypos)

        # calculate genetic distances
        gdist = self.gdist2g(vrnt_chrgrp, vrnt_genpos, rst, rsp, cst, csp)

        return gdist

    #################### Export Methods ####################
    def to_pandas_df(self):
        """
        Convert a GeneticMap object to a pandas DataFrame.

        Returns
        -------
        df : pandas.DataFrame
            A pandas DataFrame containing genetic map data.
        """
        # create a dictionary of our values
        df_dict = {
            "chr_grp" : numpy.char.decode(self._chr_grp, 'utf-8'),
            "chr_start" : self._chr_start,
            "chr_stop" : self._chr_stop,
            "map_pos" : self._map_pos,
            "mkr_name" : numpy.char.decode(self._mkr_name, 'utf-8'),
            "map_fncode" : numpy.char.decode(self._map_fncode, 'utf-8')
        }

        # make the DataFrame
        df = pandas.DataFrame(df_dict)

        # return it
        return df

    def to_csv(self, fname, sep = ',', header = True, index = False):
        # convert GeneticMap to DataFrame
        df = self.to_pandas_df()

        # write to csv
        df.to_csv(
            fname,
            sep = sep,
            header = header,
            index = index
        )

    def to_egmap(self, fname):
        # write as special instance of csv
        self.to_csv(
            fname = fname,
            sep = '\t',
            header = True,
            index = False
        )


    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def from_pandas_df(pandas_df, vrnt_chrgrp_ix = 0, vrnt_phypos_ix = 1, vrnt_stop_ix = 2, vrnt_genpos_ix = 3, vrnt_name_ix = None, vrnt_fncode_ix = None):
        """
        Read genetic map data from a CSV-like file.

        ------------------------------------------------------------=
        Genetic map file (*.gmap) format (similar to BED file format)
        ------------------------------------------------------------=
        Genetic map assumptions:
            This file format assumes that we have a high quality genome assembly
            with near complete chromosome pseudomolecules and low number of errors.
            This file format also assumes that we have a high quality genetic map
            with minimal inversions and mis-assemblies. The number of linkage groups
            in the genetic map should equal the number of whole chromosomes in our
            genome assembly.
            ***For discrepancies in the physical map vs. the genetic map, we assume
            that the physical map is correct.***
        File delimiters:
            *.gmap files are tab delimited
        File headers:
            *.gmap files are headerless.
        File fields:
            1) chrom (REQUIRED)
                Name of the chromosome; equivalent to the linkage group.
                This is of type 'str'.
            2) vrnt_phypos (REQUIRED)
                The start position of the feature on the chromosome or scaffold.
                This is 0-indexed (e.g. the first base of a chromosome is 0) and
                inclusive (e.g. vrnt_phypos <= sequence < vrnt_stop).
                This is an integer type.
            3) vrnt_stop (REQUIRED)
                The stop position of the feature on the chromosome or scaffold. This
                is 0-indexed (e.g. the first base of a chromosome is 0) and
                exclusive (e.g. vrnt_phypos <= sequence < vrnt_stop).
                This is an integer type.
            4) vrnt_genpos (REQUIRED)
                The genetic map position in Morgans. (NOT centiMorgans!)
                This is an floating type.
            5) vrnt_name (optional)
                The name of the marker on the genetic map.
                This is of type 'str'.
            6) vrnt_fncode (optional)
                The mapping function code used to create this gentic map.
                This is of type 'str'.
                Mapping function codes:
                    Haldane: 'H'
                    Kosambi: 'K'
                    Unknown: 'U'
                    Custom: <str of any length>

        Parameters
        ----------
        vrnt_chrgrp_ix : int, default=0
            Column index to specify 'chrom' (chromosome) field.
        vrnt_phypos_ix : int, default=1
            Column index to specify 'vrnt_phypos' (chromosome start) field.
        vrnt_stop_ix : int, default=2
            Column index to specify 'vrnt_stop' (chromosome end) field.
        vrnt_genpos_ix : int, default=3
            Column index to specify 'vrnt_genpos' (map position in Morgans) field.
        vrnt_name_ix : int, default=None
            Column index to specify 'vrnt_name' (marker name) field.
        vrnt_fncode_ix : int, default=None
            Column index to specify 'vrnt_fncode' (mapping function code) field.

        Returns
        -------
        genetic_map : ExtendedGeneticMap
            An ExtendedGeneticMap object containing all data required for
            genetic inferences.
        """
        # check to make sure several indices aren't None.
        pybropt.util.check_is_not_none(vrnt_chrgrp_ix, "vrnt_chrgrp_ix")
        pybropt.util.check_is_not_none(vrnt_phypos_ix, "vrnt_phypos_ix")
        pybropt.util.check_is_not_none(vrnt_stop_ix, "vrnt_stop_ix")
        pybropt.util.check_is_not_none(vrnt_genpos_ix, "vrnt_genpos_ix")

        # get data
        vrnt_chrgrp = pandas_df.iloc[:,vrnt_chrgrp_ix].values
        vrnt_phypos = pandas_df.iloc[:,vrnt_phypos_ix].values
        vrnt_stop   = pandas_df.iloc[:,vrnt_stop_ix].values
        vrnt_genpos = pandas_df.iloc[:,vrnt_genpos_ix].values
        vrnt_name   = pandas_df.iloc[:,vrnt_name_ix].values if vrnt_name_ix is not None else None
        vrnt_fncode = pandas_df.iloc[:,vrnt_fncode_ix].values if vrnt_fncode_ix is not None else None

        # force data conversion
        vrnt_chrgrp = numpy.int64([str(e) for e in vrnt_chrgrp])
        vrnt_phypos = numpy.int64(vrnt_phypos)
        vrnt_stop = numpy.int64(vrnt_stop)
        vrnt_genpos = numpy.float64(vrnt_genpos)
        vrnt_name = numpy.string_([str(e) for e in vrnt_name]) if vrnt_name is not None else None
        vrnt_fncode = numpy.string_([str(e) for e in vrnt_fncode]) if vrnt_fncode is not None else None

        # construct the gmap object
        genetic_map = ExtendedGeneticMap(
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_stop = vrnt_stop,
            vrnt_genpos = vrnt_genpos,
            vrnt_name = vrnt_name,
            vrnt_fncode = vrnt_fncode
        )

        return genetic_map

    @staticmethod
    def from_csv(fpath, sep = ',', header=0, vrnt_chrgrp_ix = 0, vrnt_phypos_ix = 1, vrnt_stop_ix = 2, vrnt_genpos_ix = 3, vrnt_name_ix = None, vrnt_fncode_ix = None):
        """
        Create an ExtendedGeneticMap object from a csv or delimited file.

        Parameters
        ----------
        fpath : str, path object, or file-like object
            Any valid string path, including URLs. Valid URL schemes include http,
            ftp, s3, and file. For file URLs, a host is expected. (see pandas docs)
        sep : str, default ','
            Delimiter to use.
        header : int, list of int, default=0
            Row number(s) to use as the column names, and the start of the data.
        vrnt_chrgrp_ix : int, default=0
            Column index to specify 'chrom' (chromosome) field.
        vrnt_phypos_ix : int, default=1
            Column index to specify 'vrnt_phypos' (chromosome start) field.
        vrnt_stop_ix : int, default=2
            Column index to specify 'vrnt_stop' (chromosome end) field.
        vrnt_genpos_ix : int, default=3
            Column index to specify 'vrnt_genpos' (map position in Morgans) field.
        vrnt_name_ix : int, default=None
            Column index to specify 'vrnt_name' (marker name) field.
        vrnt_fncode_ix : int, default=None
            Column index to specify 'vrnt_fncode' (mapping function code) field.

        Returns
        -------
        genetic_map : ExtendedGeneticMap
            An ExtendedGeneticMap object containing all data required for
            genetic inferences.

        """
        # read file using pandas.read_csv
        df = pandas.read_csv(
            fpath,
            sep=sep,
            header=header
        )

        genetic_map = ExtendedGeneticMap.from_pandas_df(
            pandas_df = df,
            vrnt_chrgrp_ix = vrnt_chrgrp_ix,
            vrnt_phypos_ix = vrnt_phypos_ix,
            vrnt_stop_ix = vrnt_stop_ix,
            vrnt_genpos_ix = vrnt_genpos_ix,
            vrnt_name_ix = vrnt_name_ix,
            vrnt_fncode_ix = vrnt_fncode_ix
        )

        return genetic_map

    @staticmethod
    def from_egmap(fpath):
        """
        Read extended genetic map file (*.egmap).

        Assumes that the file is correctly formatted.

        ------------------------------------------------------------=
        Genetic map file (*.gmap) format (similar to BED file format)
        ------------------------------------------------------------=
        Genetic map assumptions:
            This file format assumes that we have a high quality genome assembly
            with near complete chromosome pseudomolecules and low number of
            errors. This file format also assumes that we have a high quality
            genetic map with minimal inversions and mis-assemblies. The number
            of linkage groups in the genetic map should equal the number of
            whole chromosomes in our genome assembly.
            ***For discrepancies in the physical map vs. the genetic map, we
            assume that the physical map is correct.***
        File delimiters:
            *.egmap files are tab delimited
        File headers:
            *.gmap files are headerless.
        File fields:
            1) chrom (REQUIRED)
                Name of the chromosome; equivalent to the linkage group.
                This is of type 'str'.
            2) chr_start (REQUIRED)
                The start position of the feature on the chromosome or scaffold.
                This is 1-indexed (e.g. the first base of a chromosome is 1) and
                inclusive (e.g. chr_start <= sequence <= chr_stop).
                This is an integer type.
            3) chr_stop (REQUIRED)
                The stop position of the feature on the chromosome or scaffold. This
                is 1-indexed (e.g. the first base of a chromosome is 1) and
                inclusive (e.g. chr_start <= sequence <= chr_stop).
                This is an integer type.
            4) map_pos (REQUIRED)
                The genetic map position in Morgans. (NOT centiMorgans!)
                This is an floating type.
            5) mkr_name (optional)
                The name of the marker on the genetic map.
                This is of type 'str'.
            6) map_fncode (optional)
                The mapping function code used to create this gentic map.
                This is of type 'str'.
                Mapping function codes:
                    Haldane: 'H'
                    Kosambi: 'K'
                    Unknown: 'U'
                    Custom: <str of any length>

        Parameters
        ----------
        fpath : str, path object, or file-like object
            Any valid string path, including URLs. Valid URL schemes include http,
            ftp, s3, and file. For file URLs, a host is expected. (see pandas docs)

        Returns
        -------
        genetic_map : gmap
            A gmap object containing all data required for genetic inferences.
        """
        # read the file using pandas
        df = pandas.read_csv(
            fpath,          # file path
            sep = '\t',     # tab delimited
            header = 0      # there is a header
        )

        genetic_map = ExtendedGeneticMap.from_pandas_df(
            pandas_df = df,
            vrnt_chrgrp_ix = 0,
            vrnt_phypos_ix = 1,
            vrnt_stop_ix = 2,
            vrnt_genpos_ix = 3,
            vrnt_name_ix = 4 if "mkr_name" in df.columns else None,
            vrnt_fncode_ix = 5 if "map_fncode" in df.columns else None
        )

        return genetic_map



################################################################################
################################## Utilities ###################################
################################################################################
def is_ExtendedGeneticMap(v):
    return isinstance(v, ExtendedGeneticMap)

def check_is_ExtendedGeneticMap(v, varname):
    if not isinstance(v, ExtendedGeneticMap):
        raise TypeError("'%s' must be an ExtendedGeneticMap." % varname)

def cond_check_is_ExtendedGeneticMap(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_ExtendedGeneticMap(v, varname)
