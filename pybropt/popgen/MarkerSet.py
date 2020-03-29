# 3rd party libraries
import numpy

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
        doc = "The chr_grp property."
        def fget(self):
            return self._chr_grp
        def fset(self, value):
            self._chr_grp = value
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
    def lexsort(self, keys = None):
        # if no keys were provided, set a default
        if keys is None:
            keys = (self._mkr_name, self._chr_stop, self._chr_start, self._chr_grp)
        else:
            for i,k in enumerate(keys):
                pybropt.util.check_matrix_size(k, "key"+i, self.__len__())

        # filter out None; build tuple
        keys = tuple(k for k in keys if k is not None)

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
            Indices of where to place elements.
        """
        # sort internal self
        self._chr_grp = self._chr_grp[indices]
        self._chr_start = self._chr_start[indices]
        self._chr_stop = self._chr_stop[indices]
        if self._mkr_name is not None:
            self._mkr_name = self._mkr_name[indices]

    def group(self):
        """
        Calculate grouping indices.
        """
        if not self._sorted:
            raise RuntimeError("cannot group unsorted marker set")

        # get unique group names, start indices, group lengths
        uniq = numpy.unique(
            self._chr_grp,
            return_index = True,
            return_counts = True
        )

        # make assignments
        self._chr_grp_name, self._chr_grp_stix, self._chr_grp_len = uniq

        # calculate chr_grp_spix (stop indices)
        self._chr_grp_spix = self._chr_grp_stix + self._chr_grp_len

    def sort(self, keys = None):
        """
        Sort marker set.
        """
        # get indices for sort
        indices = self.lexsort(keys)

        # reorder internals
        self.reorder(indices)

        # indicate that we've sorted
        self._sorted = True

        # calculate grouping indices
        self.group()

    # TODO: accept a pattern
    def mkr_rename(self):
        new_name = numpy.core.defchararray.add(
            self._chr_grp,
            numpy.string_(["_%s" % e for e in self._chr_start])
        )

        self._mkr_name = new_name

    def has(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None):
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

    def physical_dist(self, rst, rsp, cst, csp):
        raise NotImplementedError("Physical distance not implemented yet.")
