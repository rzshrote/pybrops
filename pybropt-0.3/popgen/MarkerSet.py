# 3rd party libraries
import numpy
import math

# import out libraries
import pybropt.util

class MarkerSet:
    """docstring for MarkerSet."""
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, chr_grp, chr_start, chr_stop, mkr_name = None,
            auto_sort = False, auto_mkr_rename = False):
        # check for matrices
        pybropt.util.check_is_matrix(chr_grp, "chr_grp")
        pybropt.util.check_is_matrix(chr_start, "chr_start")
        pybropt.util.check_is_matrix(chr_stop, "chr_stop")
        pybropt.util.cond_check_is_matrix(mkr_name, "mkr_name")

        # check number of dimensions == 1
        pybropt.util.check_matrix_ndim(chr_grp, "chr_grp", 1)
        pybropt.util.check_matrix_ndim(chr_start, "chr_start", 1)
        pybropt.util.check_matrix_ndim(chr_stop, "chr_stop", 1)
        pybropt.util.cond_check_matrix_ndim(mkr_name, "mkr_name", 1)

        # check length of matrix == chr_grp.size
        pybropt.util.check_matrix_size(chr_start, "chr_start", chr_grp.size)
        pybropt.util.check_matrix_size(chr_stop, "chr_stop", chr_grp.size)
        pybropt.util.cond_check_matrix_size(mkr_name, "mkr_name", chr_grp.size)

        # check matrix dtypes
        pybropt.util.check_matrix_dtype_is_string_(chr_grp, "chr_grp")
        pybropt.util.check_matrix_dtype_is_integer(chr_start, "chr_start")
        pybropt.util.check_matrix_dtype_is_integer(chr_stop, "chr_stop")
        pybropt.util.cond_check_matrix_dtype_is_string_(mkr_name, "mkr_name")

        # set private variables
        self._chr_grp = chr_grp
        self._chr_start = chr_start
        self._chr_stop = chr_stop
        self._mkr_name = mkr_name
        self._chr_grp_name = None
        self._chr_grp_len = None
        self._chr_grp_stix = None
        self._chr_grp_spix = None

        # give default marker names if needed
        if auto_mkr_rename:
            self.mkr_rename()

        # sort if needed
        if auto_sort:
            self.sort()

    def __len__(self):
        """
        Get length of marker set.
        """
        return len(self._chr_grp)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def chr_grp():
        doc = "Chromosome group array access. See chr_grp.{fget,fset,fdel} function docs."
        def fget(self):
            """
            Retrieve a pointer to the chromosome group array.

            Returns
            -------
            chr_grp : numpy.array
                Pointer to the chromosome group array.
            """
            return self._chr_grp
        def fset(self, chr_grp):
            """
            Set the pointer to the chromosome group array. Perform error checks.

            Parameters
            ----------
            chr_grp : numpy.array
                A numpy.string_ matrix 
            """
            self._chr_grp = chr_grp
        def fdel(self):
            del self._chr_grp
        return locals()
    chr_grp = property(**chr_grp())

    def chr_start():
        doc = "The chr_start property."
        def fget(self):
            return self._chr_start
        def fset(self, value):
            self._chr_start = value
        def fdel(self):
            del self._chr_start
        return locals()
    chr_start = property(**chr_start())

    def chr_stop():
        doc = "The chr_stop property."
        def fget(self):
            return self._chr_stop
        def fset(self, value):
            self._chr_stop = value
        def fdel(self):
            del self._chr_stop
        return locals()
    chr_stop = property(**chr_stop())

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

    def chr_grp_stix():
        doc = "The chr_grp_stix property."
        def fget(self):
            return self._chr_grp_stix
        def fset(self, value):
            self._chr_grp_stix = value
        def fdel(self):
            del self._chr_grp_stix
        return locals()
    chr_grp_stix = property(**chr_grp_stix())

    def chr_grp_spix():
        doc = "The chr_grp_spix property."
        def fget(self):
            return self._chr_grp_spix
        def fset(self, value):
            self._chr_grp_spix = value
        def fdel(self):
            del self._chr_grp_spix
        return locals()
    chr_grp_spix = property(**chr_grp_spix())

    def chr_grp_name():
        doc = "The chr_grp_name property."
        def fget(self):
            return self._chr_grp_name
        def fset(self, value):
            self._chr_grp_name = value
        def fdel(self):
            del self._chr_grp_name
        return locals()
    chr_grp_name = property(**chr_grp_name())

    def chr_grp_len():
        doc = "The chr_grp_len property."
        def fget(self):
            return self._chr_grp_len
        def fset(self, value):
            self._chr_grp_len = value
        def fdel(self):
            del self._chr_grp_len
        return locals()
    chr_grp_len = property(**chr_grp_len())

    def sorted():
        doc = "The sorted property."
        def fget(self):
            return self._sorted
        def fset(self, value):
            self._sorted = value
        def fdel(self):
            del self._sorted
        return locals()
    sorted = property(**sorted())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys = None):
        """
        Perform an indirect stable sort using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using chr_grp as primary key, chr_start as
            secondary key, chr_stop as tertiary key, mkr_name as quaternary key.

        Returns
        -------
        indices : numpy.ndarray
            Array of indices that sort the keys.
            Remark: does not check length of 'indices' with length of 'self'.
        """
        # if no keys were provided, set a default
        if keys is None:
            # chr_grp 1st, chr_start 2nd, chr_stop 3rd, mkr_name 4th
            keys = [self._mkr_name, self._chr_stop, self._chr_start, self._chr_grp]

        # filter out None (if mkr_name is None); build list
        keys = list(k for k in keys if k is not None)

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def reorder(self, indices):
        """
        Reorder the marker set.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements. Length must equal of
        """
        # make sure len(indices) == len(self)
        pybropt.util.check_equal_len(indices, "indices", self, "self")

        # reorder internal arrays
        self._chr_grp = self._chr_grp[indices]
        self._chr_start = self._chr_start[indices]
        self._chr_stop = self._chr_stop[indices]
        if self._mkr_name is not None:
            self._mkr_name = self._mkr_name[indices]

    def group(self):
        """
        Calculate chromosome grouping indices (group by chr_grp).
        If not sorted, raise RuntimeError.
        """
        # raise RuntimeError if not sorted.
        if not self._sorted:
            raise RuntimeError("cannot group unsorted marker set")

        # get unique group names, start indices, group lengths
        uniq = numpy.unique(self._chr_grp, return_index = True, return_counts = True)

        # make assignments
        self._chr_grp_name, self._chr_grp_stix, self._chr_grp_len = uniq

        # calculate chr_grp_spix (stop indices)
        self._chr_grp_spix = self._chr_grp_stix + self._chr_grp_len

    def sort(self, keys = None):
        """
        Sort the MarkerSet object using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using chr_grp as primary key, chr_start as
            secondary key, chr_stop as tertiary key, mkr_name as quaternary key.
        """
        # get indices for sort
        indices = self.lexsort(keys)

        # reorder internals
        self.reorder(indices)

        # indicate that we've sorted
        self._sorted = True

        # calculate grouping indices
        self.group()

    ################ Insert/Delete Methods #################
    def remove(self, indices, auto_sort = False):
        """
        Remove indices from the MarkerSet.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
        auto_sort : bool, default = False
            Whether to sort the array or not after removing indices.
        """
        # delete indices from self
        self._chr_grp = numpy.delete(self._chr_grp, indices)
        self._chr_start = numpy.delete(self._chr_start, indices)
        self._chr_stop = numpy.delete(self._chr_stop, indices)
        if self._mkr_name is not None:
            self._mkr_name = numpy.delete(self._mkr_name, indices)

        # sort if needed
        if auto_sort:
            self.sort()

    def select(self, indices, auto_sort = False):
        """
        Keep only selected markers, removing all others from the MarkerSet.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
        auto_sort : bool, default = False
            Whether to sort the array or not after selecting indices.
        """
        # keep only selected markers.
        self._chr_grp = self._chr_grp[indices]
        self._chr_start = self._chr_start[indices]
        self._chr_stop = self._chr_stop[indices]
        if self._mkr_name is not None:
            self._mkr_name = self._mkr_name[indices]

        # sort if needed
        if auto_sort:
            self.sort()

    def prune(self, nt, auto_sort = False):
        """
        Prune markers evenly across all chromosomes.

        Parameters
        ----------
        nt : int
            Target distance between each selected marker in nucleotides.
        auto_sort : bool, default = False
            Whether to sort the array or not after selecting indices.
        """
        # if not sorted, sort
        if not self._sorted:
            self.sort()

        # make empty index list to store selected marker indices
        indices = []

        # for each chromosome
        for st,sp in zip(self._chr_grp_stix, self._chr_grp_spix):
            # calculate chromosome length given start, end marker
            chr_len = self._chr_start[sp-1] - self._chr_start[st]

            # calculate the target distance between each marker (float)
            step = chr_len / int(math.ceil(chr_len / nt))

            # force addition of the first marker on the chromosome
            indices.append(st)

            # target site; we want markers as close to this value (float)
            target =  self._chr_start[st] + step

            # for each locus index in the chromosome
            for i in range(st+1, sp):
                # if position exceeds target, determine which marker to add
                if self._chr_start[i] >= target:
                    # get distance between target and previous marker
                    downstream = target - self._chr_start[i-1]

                    # get distance between target and current marker
                    upstream = self._chr_start[i] - target

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

        # convert indices into an array
        indices = numpy.array(indices)

        # select the desired marker indices
        self.select(indices, auto_sort)

    # TODO: def insert(self, ...)

    ################### Distance Methods ###################
    # FIXME: will not work with distances that span different chromosomes.
    def physical_dist(self, rst, rsp, cst, csp):
        raise NotImplementedError

    #################### Marker Methods ####################
    # TODO: accept a pattern
    def mkr_rename(self, new_mkr_name = None):
        # if new names have not been provided, auto generate them
        if new_mkr_name is None:
            new_mkr_name = numpy.core.defchararray.add(
                self._chr_grp,
                numpy.string_(["_%s" % e for e in self._chr_start])
            )

        # check for correct type and dimensions
        pybropt.util.check_matrix_ndim(new_mkr_name, "new_mkr_name", 1)
        pybropt.util.check_matrix_size(new_mkr_name, "new_mkr_name", self.__len__())
        pybropt.util.check_matrix_dtype_is_string_(new_mkr_name, "new_mkr_name")

        # set new marker names
        self._mkr_name = new_name

    def mkr_mask(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None, invert = False):
        """
        Make a mask of entries in the MarkerSet that have the designated
        markers. Markers fitting a *single* search criteria are logical and'ed
        together into a single mask.

        Parameters
        ----------
        chr_grp : str, list of str
            String or list of strings indicating which chromosome group to
            select.
        chr_start : int, list of int
            Integer or list of integers indicating start position(s) to select.
            Length of chr_start and chr_stop must be equal.
            Positions are compared to marker starting positions:
                self._chr_start >= chr_start
        chr_stop : int, list of int
            Integer or list of integers indicating stop position(s) to select.
            Length of chr_start and chr_stop must be equal.
            Positions are compared to marker starting positions:
                self._chr_start <= chr_stop
        mkr_name : str, list of str
            String or list of strings indicating which marker names to select.
        invert : bool
            Whether to invert mask or not after mask construction.

        Returns
        -------
        mask : numpy.ndarray
            A mask array of search results.
        """
        # check input data types; raise error if they're bad.
        pybropt.util.cond_check_is_string_or_iterable(chr_grp, "chr_grp")
        pybropt.util.cond_check_is_integer_or_iterable(chr_start, "chr_start")
        pybropt.util.cond_check_is_integer_or_iterable(chr_stop, "chr_stop")
        pybropt.util.cond_check_is_string_or_iterable(mkr_name, "mkr_name")
        pybropt.util.check_is_bool(invert, "invert")

        # convert any integers to list
        if (chr_start is not None) and numpy.issubdtype(type(chr_start), numpy.integer):
            chr_start = [chr_start]
        if (chr_stop is not None) and numpy.issubdtype(type(chr_stop), numpy.integer):
            chr_stop = [chr_stop]

        # make an empty list of masks
        masks = []

        # find matches for chr_grp
        if chr_grp is not None:
            masks.append(numpy.in1d(self._chr_grp, chr_grp))

        # find matches for chr_start, chr_stop ranges
        if (chr_start is not None) and (chr_stop is not None):
            # make sure length of provided lists or iterables are equivalent
            pybropt.util.check_equal_len(chr_start, "chr_start", chr_stop, "chr_stop")

            # for each pair of start and stop positions
            for st,sp in zip(chr_start, chr_stop):
                # make mask for current start, stop pair
                m = numpy.logical_and(self._chr_start >= st, self._chr_start <= sp)

                # append to list of masks
                masks.append(m)

        # find matches for mkr_name
        if mkr_name is not None:
            masks.append(numpy.in1d(self._mkr_name, mkr_name))

        # default value of None (for cases where len(masks) == 0)
        mask = None

        # if the len(masks) > 0, logical_and merge them together
        if len(masks) > 0:
            mask = numpy.logical_and.reduce(masks)

            # invert mask if we need to
            if invert:
                mask = ~mask

        # return mask
        return mask

    # TODO: make this so that indices are matched with the arguments?
    def mkr_index(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None, invert = False):
        # find markers
        mask = self.mkr_mask(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            invert = invert
        )

        # get indices
        ix = numpy.flatnonzero(mask) if mask is not None else None

        return ix

    def mkr_remove(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None, invert = False, auto_sort = False):
        """
        Remove markers matching a unique signature from the marker set.
        """
        # get markers
        mask = self.mkr_mask(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            invert = False
        )

        self.remove(mask, auto_sort)

    # TODO: accept chr_grp, chr_start, chr_stop; support for invert
    def mkr_where(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None):
        """
        Find indices for where a marker is. Indices match the order of the
        input.
        """
        # initialize tuple of empty variables
        # these represent chr_grp_ix, chr_start_ix, chr_stop_ix, mkr_name_ix
        finds = [None, None, None, None]

        # if chromosome group given, build locations arrays
        if chr_grp is not None:
            finds[0] = [numpy.flatnonzero(self._chr_grp == e) for e in chr_grp]

        # if chromosome start given, build start location arrays
        if chr_start is not None:
            finds[1] = [numpy.flatnonzero(self._chr_start == e) for e in chr_start]

        # if chromosome stop given, build stop location arrays
        if chr_stop is not None:
            finds[2] = [numpy.flatnonzero(self._chr_stop == e) for e in chr_stop]

        # if marker name given, build marker name arrays
        if mkr_name is not None:
            finds[3] = [numpy.flatnonzero(self._mkr_name == e) for e in mkr_name]

        # remove None from list of finds
        finds = tuple(f for f in finds if f is not None)

        # overlap the finds
        overlap = None
        if len(finds) > 1:
            overlap = finds[0]
            for i in range(1, len(finds)):
                overlap = [numpy.intersect1d(overlap[j], e) for j,e in enumerate(finds[i])]
            # if we are searchign for overlapping keys, intersect
            # overlap = [functools.reduce(numpy.intersect1d, e) for e in zip(*finds)]
        elif len(finds) == 1:
            # if we are just searching for one key, grab only element
            overlap = finds[0]
        else:
            # otherwise, return None
            return overlap

        # check for multi-mapping + build where array
        where = numpy.empty(len(overlap), dtype = 'int64')
        for i,o in enumerate(overlap):
            # get length of overlap numpy.ndarray
            l = len(o)

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
            where[i] = l if l > 0 else -1

        # return where array
        return where

    # TODO: work on this since chr_start could be on multiple chromosomes
    def mkr_has(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None):
        """
        Determine if a genetic map has a specific marker in it.
        If multiple arguments are specified, a logical and is applied to all
        fields.

        Parameters
        ----------
        chr_grp : numpy.ndarray, None
            An array of chromosome groups to determine if in self. Optional.
        chr_start : numpy.ndarray
            An array of chromosome start positions to determine if in self.
            Optional.
        chr_stop : numpy.ndarray
            An array of chromosome start positions to determine if in self.
            Optional.
        map_pos : numpy.ndarray
            An array of genetic map positions to determine if in self. Optional.
        mkr_name : numpy.ndarray
            An array of marker names to determine if in self. Optional.
        map_fncode : numpy.ndarray
            An array of genetic map function codes to determine if in self.
            Optional.

        Returns
        -------
        mask : numpy.ndarray
            A boolean mask of whether the provided information corresponds to
            a marker within the genetic map object.
        """
        # make a tuple of None values (one for each argument)
        masks = (None, None, None, None)

        # test whether chr_grp is in self._chr_grp_name
        if chr_grp is not None:
            masks[0] = numpy.in1d(chr_grp, self._chr_grp_name)
        # test whether chr_start is in self._chr_start
        if chr_start is not None:
            masks[1] = numpy.in1d(chr_start, self._chr_start)
        # test whether chr_stop is in self._chr_stop
        if chr_stop is not None:
            masks[2] = numpy.in1d(chr_stop, self._chr_stop)
        # test whether mkr_name is in self._mkr_name
        if mkr_name is not None:
            masks[3] = numpy.in1d(mkr_name, self._mkr_name)

        # filter out None
        masks = tuple(m for m in masks if m is not None)

        # default value of None (for cases where len(masks) == 0)
        mask = None

        # if the len(masks) > 0, logical_and merge them together
        if len(masks) > 0:
            mask = numpy.logical_and.reduce(masks)

        # return mask
        return mask

    ################### Copying Methods ####################
    def copy(self, deepcopy = True):
        return self.from_self(deepcopy = True)

    def from_self(self, deepcopy = True,
        chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None,
        select_ix = None, remove_ix = None,
        auto_sort = False, auto_mkr_rename = False):
        """
        Create a copy of self and optionally alter select variables.

        Parameters
        ----------
        deepcopy : bool, default = True
            Perform a deep copy of internal arrays if True. If False, only copy
            pointers (other objects sharing the same array pointer could alter
            array values in an undesireable manner).
            If chr_grp, chr_start, chr_stop, mkr_name are provided, elements
            from these arrays are deep copied.
        chr_grp : numpy.ndarray
            Optional chr_grp array to replace in the new MarkerSet.
        chr_start : numpy.ndarray
            Optional chr_start array to replace in the new MarkerSet.
        chr_stop : numpy.ndarray
            Optional chr_stop array to replace in the new MarkerSet.
        mkr_name : numpy.ndarray
            Optional mkr_name array to replace in the new MarkerSet.
        select_ix : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
            If chr_grp, chr_start, chr_stop, mkr_name are provided, indices from
            these arrays are selected at the corresponding location(s).
        remove_ix : numpy.ndarray, slice, int
            Array of shape (a,), slice or int of item(s) to remove.
            Where:
                'a' is the number of indices to remove.
            If chr_grp, chr_start, chr_stop, mkr_name are provided, indices from
            these arrays are removed at the corresponding location(s).

        Returns
        -------
        marker_set : MarkerSet
            A MarkerSet object.
        """
        if chr_grp is None:
            chr_grp = self._chr_grp
        if chr_start is None:
            chr_start = self._chr_start
        if chr_stop is None:
            chr_stop = self._chr_stop
        if mkr_name is None:
            mkr_name = self._mkr_name

        # test if we are selecting, removing, and finally deep copying.
        if select_ix is not None:
            # if remove_ix has also been provided, raise error
            if remove_ix is not None:
                raise ValueError("'select_ix' and 'remove_ix' cannot be provided at the same time")

            # select indices; this copies data leaving the original intact
            chr_grp = chr_grp[select_ix]
            chr_start = chr_start[select_ix]
            chr_stop = chr_stop[select_ix]
            mkr_name = mkr_name[select_ix]

        elif remove_ix is not None:
            # remove indices; this copies data leaving the original intact
            chr_grp = numpy.delete(chr_stop, remove_ix)
            chr_start = numpy.delete(chr_start, remove_ix)
            chr_stop = numpy.delete(chr_stop, remove_ix)
            mkr_name = numpy.delete(mkr_name, remove_ix)

        elif deepcopy:
            # copy arrays if we are not selecting or removing indices
            chr_grp = chr_grp.copy()
            chr_start = chr_start.copy()
            chr_stop = chr_stop.copy()
            mkr_name = mkr_name.copy()

        # make a new object instance
        marker_set = self.__class__(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            auto_sort = auto_sort,
            auto_mkr_rename = auto_mkr_rename
        )

        return marker_set