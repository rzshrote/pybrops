# import 3rd party
import math
import numpy
import pandas
from scipy.interpolate import interp1d

# import our libraries
from . import MarkerSet
# import pybropt.popgen.MarkerSet # causes problems
import pybropt.util

class GeneticMap(MarkerSet):
    """
    A class to represent genetic maps and perform associated operations.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, chr_grp, chr_start, chr_stop, map_pos,
            mkr_name = None, map_fncode = None, mapfn = 'haldane',
            auto_sort = False, auto_mkr_rename = False, auto_fncode = False,
            auto_spline = False, kind = 'linear', fill_value = 'extrapolate'
            ):
        # calls super constructor
        super(GeneticMap, self).__init__(
            chr_grp,
            chr_start,
            chr_stop,
            mkr_name,
            auto_sort = False,
            auto_mkr_rename = False
        )

        # check for matricies
        pybropt.util.check_is_matrix(map_pos, "map_pos")
        pybropt.util.cond_check_is_matrix(map_fncode, "map_fncode")

        # check number of dimensions == 1
        pybropt.util.check_matrix_ndim(map_pos, "map_pos", 1)
        pybropt.util.cond_check_matrix_ndim(map_fncode, "map_fncode", 1)

        # check length of matrix == chr_grp.size
        pybropt.util.check_matrix_size(map_pos, "map_pos", chr_grp.size)
        pybropt.util.cond_check_matrix_size(map_fncode, "map_fncode", chr_grp.size)

        # check matrix dtypes
        pybropt.util.check_matrix_dtype_is_numeric(map_pos, "map_pos")
        pybropt.util.cond_check_matrix_dtype_is_string_(map_fncode, "map_fncode")

        # if mapfn is a string, make it a lowercase
        mapfn = pybropt.util.cond_str_lower(mapfn)

        # set private variables or fill with None
        self._map_pos = map_pos
        self._map_fncode = map_fncode
        self._chr_grp_name = None
        self._chr_grp_len = None
        self._chr_grp_stix = None
        self._chr_grp_spix = None
        self._map_spline = None
        self._map_spline_kind = None
        self._map_spline_fill_value = None
        self._xo_prob = None
        self._xo_prob_len = None
        self._xo_prob_stix = None
        self._xo_prob_spix = None
        self._sorted = False

        # set mapfn: convert to callable, otherwise set to None
        self._mapfn = GeneticMap.key_to_mapfn(mapfn)

        # rename markers if desired
        if auto_mkr_rename:
            self.mkr_rename()

        # if True, auto extract code from mapfn using KEY_TO_MAP_FNCODE dict
        if auto_fncode is True:
            auto_fncode = GeneticMap.KEY_TO_MAP_FNCODE.get(mapfn, 'u') # default == 'u'

        # if auto_fncode is a string, fill map_fncode with its value
        if isinstance(auto_fncode, str):
            self._map_fncode = numpy.repeat(
                numpy.string_(auto_fncode),
                self.__len__()
            )

        # sort if needed
        if auto_sort:
            self.sort() # sort
            self.remove_discrepancies() # remove any discrepancies in map
            self.calc_xo_prob() # calculate crossover probabilities

        # build spline if needed
        if auto_spline:
            self.build_spline(kind, fill_value)

    def __len__(self):
        """
        Get the length of the genetic map
        """
        return len(self._map_pos)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def map_pos():
        doc = "The map_pos property."
        def fget(self):
            return self._map_pos
        def fset(self, value):
            self._map_pos = value
        def fdel(self):
            del self._map_pos
        return locals()
    map_pos = property(**map_pos())

    def map_fncode():
        doc = "The map_fncode property."
        def fget(self):
            return self._map_fncode
        def fset(self, value):
            self._map_fncode = value
        def fdel(self):
            del self._map_fncode
        return locals()
    map_fncode = property(**map_fncode())

    def map_spline():
        doc = "The map_spline property."
        def fget(self):
            return self._map_spline
        def fset(self, value):
            self._map_spline = value
        def fdel(self):
            del self._map_spline
        return locals()
    map_spline = property(**map_spline())

    def map_spline_kind():
        doc = "The map_spline_kind property."
        def fget(self):
            return self._map_spline_kind
        def fset(self, value):
            self._map_spline_kind = value
        def fdel(self):
            del self._map_spline_kind
        return locals()
    map_spline_kind = property(**map_spline_kind())

    def map_spline_fill_value():
        doc = "The map_spline_fill_value property."
        def fget(self):
            return self._map_spline_fill_value
        def fset(self, value):
            self._map_spline_fill_value = value
        def fdel(self):
            del self._map_spline_fill_value
        return locals()
    map_spline_fill_value = property(**map_spline_fill_value())

    def mapfn():
        doc = "The mapfn property."
        def fget(self):
            return self._mapfn
        def fset(self, value):
            self._mapfn = value
        def fdel(self):
            del self._mapfn
        return locals()
    mapfn = property(**mapfn())

    def xo_prob():
        doc = "The xo_prob property."
        def fget(self):
            return self._xo_prob
        def fset(self, value):
            self._xo_prob = value
        def fdel(self):
            del self._xo_prob
        return locals()
    xo_prob = property(**xo_prob())

    def xo_prob_len():
        doc = "The xo_prob_len property."
        def fget(self):
            return self._xo_prob_len
        def fset(self, value):
            self._xo_prob_len = value
        def fdel(self):
            del self._xo_prob_len
        return locals()
    xo_prob_len = property(**xo_prob_len())

    def xo_prob_stix():
        doc = "The xo_prob_stix property."
        def fget(self):
            return self._xo_prob_stix
        def fset(self, value):
            self._xo_prob_stix = value
        def fdel(self):
            del self._xo_prob_stix
        return locals()
    xo_prob_stix = property(**xo_prob_stix())

    def xo_prob_spix():
        doc = "The xo_prob_spix property."
        def fget(self):
            return self._xo_prob_spix
        def fset(self, value):
            self._xo_prob_spix = value
        def fdel(self):
            del self._xo_prob_spix
        return locals()
    xo_prob_spix = property(**xo_prob_spix())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys = None):
        # if no keys were provided, set a default
        if keys is None:
            keys = (
                self._map_fncode,
                self._mkr_name,
                self._chr_stop,
                self._map_pos,      # 3rd priority
                self._chr_start,    # 2nd priority
                self._chr_grp       # 1st priority
            )
        else:
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
        self._chr_grp = self._chr_grp[indices]
        self._chr_start = self._chr_start[indices]
        self._chr_stop = self._chr_stop[indices]
        self._map_pos = self._map_pos[indices]
        if self._mkr_name is not None:
            self._mkr_name = self._mkr_name[indices]
        if self._map_fncode is not None:
            self._map_fncode = self._map_fncode[indices]

    # group(self) not in need of being overridden

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

        # remove any discrepancies in map
        self.remove_discrepancies()

        # calculate crossover probabilities
        self.calc_xo_prob()

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
        self._map_pos = numpy.delete(self._map_pos, indices)
        if self._mkr_name is not None:
            self._mkr_name = numpy.delete(self._mkr_name, indices)
        if self._map_fncode is not None:
            self._map_fncode = numpy.delete(self._map_fncode, indices)

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
        self._map_pos = self._map_pos[indices]
        if self._mkr_name is not None:
            self._mkr_name = self._mkr_name[indices]
        if self._map_fncode is not None:
            self._map_fncode = self._map_fncode[indices]

        # sort if needed
        if auto_sort:
            self.sort()

    def prune(self, nt = None, M = None, auto_sort = False):
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
        auto_sort : bool, default = False
            Whether to sort the array or not after selecting indices.
        """
        # check if we have acceptible inputs
        if (nt is None) and (M is None):
            raise ValueError("'nt' and 'M' cannot both be None")

        # if not sorted, sort
        if not self._sorted:
            self.sort()

        # make empty index list to store selected marker indices
        indices = []

        # make a generic array pointer to a position; this position can be a
        # physical position (self._chr_start) or a genetic position
        # (self._map_pos); this is used in initial marker selection below.
        # genetic position takes dominance
        position = self._map_pos if M is not None else self._chr_start

        # generic spacing variable
        spacing = M if M is not None else nt

        # for each chromosome
        for st,sp in zip(self._chr_grp_stix, self._chr_grp_spix):
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
                if self._chr_grp[up] == self._chr_grp[down]:
                    # calculate physical distance between two selected markers
                    dist = self._chr_start[down] - self._chr_start[up]

                    # if we exceed 'nt' distance
                    if dist > nt:
                        # calculate the target distance between each marker (float)
                        step = dist / int(math.ceil(dist / nt))

                        # target site; we want markers as close to this value (float)
                        target = self._chr_start[up] + step

                        # for each locus between upstream and downstream markers
                        for i in range(up+1, down):
                            # if position exceeds target, determine which marker to add
                            if self._chr_start[i] >= target:
                                # get distance between target and previous marker
                                downstream = target - self._chr_start[i-1]

                                # get distance between target and current marker
                                upstream = self._chr_start[i] - target

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
        self.select(indices, auto_sort)

    # TODO: def insert(self, ...)

    ################## Integrity Methods ###################
    def map_concordancy(self):
        """
        Assess physical and genetic map concordancies.

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
        # make empty bool vector
        concordancy = numpy.zeros(len(self._chr_start), dtype = 'bool')

        # for each linkage group
        for st,sp in zip(self._chr_grp_stix, self._chr_grp_spix):
            # we assume the first element is correctly placed
            concordancy[st] = True
            # test if the current element is <= the previous element
            concordancy[st+1:sp] = self._map_pos[st:sp-1] <= self._map_pos[st+1:sp]

        # return results
        return concordancy

    def all_concordant(self):
        concordant = numpy.all(self.map_concordancy())
        return concordant

    def remove_discrepancies(self):
        """
        Remove discrepancies between the physical map and the genetic map.
        In instances of conflict, assume that the physical map is correct.

        Note:
            This assumption may cause major issues if there are incorrect
            markers at the beginning of the chromosome.
        """
        # get boolean mask for concordant map positions
        concordancy = self.map_concordancy()

        # if all are concordant, exit this function
        if( numpy.all(concordancy) ):
            return

        # for each chr_grp, set the number of concordant positions (new len)
        for i,(st,sp) in enumerate(zip(self._chr_grp_stix, self._chr_grp_spix)):
            self._chr_grp_len[i] = concordancy[st:sp].sum()

        # next we physically remove discrepancies
        self._chr_grp_stix  = self._chr_grp_len.cumsum() - self._chr_grp_len
        self._chr_grp_spix  = self._chr_grp_len.cumsum()
        self._chr_grp       = self._chr_grp[concordancy]
        self._chr_start     = self._chr_start[concordancy]
        self._chr_stop      = self._chr_stop[concordancy]
        self._map_pos       = self._map_pos[concordancy]
        self._mkr_name      = self._mkr_name[concordancy]
        self._map_fncode    = self._map_fncode[concordancy]

    ############# Distance/Probability Methods #############
    # FIXME: will not work with probabilities that span different chromosomes.
    def genetic_dist(self, rst, rsp, cst, csp):
        """
        Get a Morgan distance matrix chunk given coordinates.

        Parameters
        ----------
        rst : int
            Row start index (inclusive).
            Must fall within the length of the GeneticMap.
        rsp : int
            Row stop position (exclusive).
            Must fall within the length of the GeneticMap.
        cst : int
            Column start position (inclusive).
            Must fall within the length of the GeneticMap.
        csp : int
            Column stop position (exclusive).
            Must fall within the length of the GeneticMap.

        Returns
        -------
        d : numpy.ndarray
            A ((rst-rsp), (cst-csp)) distance matrix of every pair of markers
            in the designated GeneticMap regions.

        Example
        -------
            [in]  > genetic_map.map_pos
            [out] > [0,1,2,3,4,5,6]
            [in]  > genetic_map.dist(0,2,3,5)
                  # [[abs(0-3), abs(1-3)],
                  #  [abs(0-4), abs(1-4)]]
            [out] > [[3, 2],
                    [4, 3]]
        """
        # make an mesh grid for row and column gmap's
        rowmesh, colmesh = numpy.meshgrid(
            self._map_pos[rst:rsp],
            self._map_pos[cst:csp],
            indexing='ij'
        )

        # take the difference of the two components in the grid and return abs()
        d = numpy.absolute(rowmesh - colmesh)

        # return distances
        return d

    # FIXME: will not work with probabilities that span different chromosomes.
    def recomb_prob(self, rst, rsp, cst, csp, mapfn = None):
        """
        Get a Morgan distance matrix chunk given coordinates.

        Parameters
        ----------
        rst : int
            Row start index (inclusive).
            Must fall within the length of the GeneticMap.
        rsp : int
            Row stop position (exclusive).
            Must fall within the length of the GeneticMap.
        cst : int
            Column start position (inclusive).
            Must fall within the length of the GeneticMap.
        csp : int
            Column stop position (exclusive).
            Must fall within the length of the GeneticMap.
        mapfn : str, callable
            Name of mapping function or custom mapping function.

        Returns
        -------
        r : numpy.ndarray
            A ((rst-rsp), (cst-csp)) recombination probability matrix of every
            pair of markers in the designated GeneticMap regions.
        """
        # convert mapfn to a callable function
        if mapfn is None:
            mapfn = self._mapfn
        elif callable(mapfn):
            pass
        else:
            mapfn = GeneticMap.KEY_TO_MAPFN.get(mapfn, None)

        # check if we have a callable function after the chain above
        if not callable(mapfn):
            raise ValueError("Unrecognized type for 'mapfn'")

        # calculate genetic distance
        d = self.genetic_dist(rst, rsp, cst, csp)

        # convert to recombination probability
        r = mapfn(d)

        return r

    def calc_xo_prob(self, mapfn = None):
        """
        Calculate vector of crossover probabilities for between neighbors.
        This sets the internal xo_prob, xo_prob_len, xo_prob_stix, and
        xo_prob_spix variables.

        Parameters
        ----------
        mapfn : str, callable
            String name for the mapping function to use, or a callable mapping
            function.

        Returns
        -------
        xo_prob : numpy.ndarray
            A vector of crossover probabilities between each neighbor. This is
            the same vector as self.xo_prob.
        """
        # prerequisite of sorted internals
        if not self._sorted:
            self.sort()

        # convert mapfn to a callable function
        if mapfn is None:
            mapfn = self._mapfn
        else:
            mapfn = GeneticMap.key_to_mapfn(mapfn)

        # since probabilities are between elements, we subtract 1 from chr len
        self._xo_prob_len = self._chr_grp_len - 1

        # get cumsum for stop indices
        self._xo_prob_spix = self._xo_prob_len.cumsum()

        # subtract lengths to get the start indices
        self._xo_prob_stix = self._xo_prob_spix - self._xo_prob_len

        # allocate empty vector for our crossover probability vector
        self._xo_prob = numpy.empty(
            self._xo_prob_len.sum(),
            dtype = 'float64'
        )

        # calculate probabilities for recombination between neighbors
        for gst, gsp, pst, psp in zip(
                self._chr_grp_stix, self._chr_grp_spix,
                self._xo_prob_stix, self._xo_prob_spix
        ):
            # Step 1) take n+1 map positions and subtract n map positions
            # Step 2) then take the probability of the genetic distance.
            self._xo_prob[pst:psp] = mapfn(
                self._map_pos[gst+1:gsp] - self._map_pos[gst:gsp-1]
            )

        return self._xo_prob

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
        # get unique chr_grp and their positions
        uniq_chr_grp, uniq_inv = numpy.unique(
            self._chr_grp,
            return_inverse = True
        )

        # make an empty dictionary for map splines
        self._map_spline = {}

        # set the spline type
        self._map_spline_kind = kind
        self._map_spline_fill_value = fill_value

        # iterate through unique groups and build spline
        for i,grp in enumerate(uniq_chr_grp):
            # get mask of elements in the chr_grp
            mask = uniq_inv == i

            # build spline and assign to dictionary
            self._map_spline[grp] = interp1d(
                x = self._chr_start[mask],  # chromosome positions
                y = self._map_pos[mask],    # map positions
                kind = kind,                # type of interpolation
                fill_value = fill_value,    # default fill value
                assume_sorted = False       # arrays are not sorted
            )

    def interpolate(self, chr_grp, chr_start, chr_stop, mkr_name = None,
            map_fncode = None, auto_sort = False, auto_mkr_rename = False,
            kind = None, fill_value = None):
        """
        Interpolate a genetic map from the current one.
        """
        # if kind and fill_value are none, set them to what is internal.
        if kind is None:
            kind = self._map_spline_kind
        if fill_value is None:
            fill_value = self._map_spline_fill_value

        # test if either kind or fill_value is None

        # build spline if we need it
        if (kind != self._map_spline_kind) or (fill_value != self._map_spline_fill_value):
            self.build_spline(kind, fill_value)
        elif (kind is None) or (fill_value is None):
            raise RuntimeError("'kind' and 'fill_value' cannot be None")

        # allocate memory for map_pos vector
        map_pos = numpy.empty(len(chr_grp), dtype = 'float64')

        # for each position, interpolate
        for i,(grp,start) in enumerate(zip(chr_grp, chr_start)):
            # get spline function, else return lambda that just return NaN
            splinefn = self._map_spline.get(grp, lambda x: numpy.nan)

            # make assignment
            map_pos[i] = splinefn(start)

        # if there is no map_fncode attached, we assign it to
        if map_fncode is None:
            map_fncode = numpy.repeat(
                numpy.string_([GeneticMap.KEY_TO_MAP_FNCODE["interpolate"]]),
                len(chr_grp)
            )

        # make the genetic map
        genetic_map = GeneticMap(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            map_pos = map_pos,
            mkr_name = mkr_name,
            map_fncode = map_fncode,
            mapfn = self._mapfn,
            auto_sort = auto_sort,
            auto_mkr_rename = auto_mkr_rename,
            auto_fncode = False,
            auto_spline = False,
            kind = kind,
            fill_value = fill_value
        )

        return genetic_map

    # TODO: rewrite me
    # def interpolate_pandas_df(self, pandas_df, chr_grp_ix = 0, chr_start_ix = 1,
    #         chr_stop_ix = 2, mkr_name_ix = None, map_fncode_ix = None,
    #         kind = 'linear', fill_value = 'extrapolate'):

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

    # TODO: remove me? this has been replaced by mkr_remove
    def has(self, chr_grp = None, chr_start = None, chr_stop = None,
            map_pos = None, mkr_name = None, map_fncode = None
            ):
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
        # get mask from super function
        supermask = super(GeneticMap, self).has(
            chr_grp,
            chr_start,
            chr_stop,
            mkr_name
        )

        # make a tuple of None values (one for each argument)
        masks = (supermask, None, None)

        # test whether map_pos is in self._map_pos
        if map_pos is not None:
            masks[1] = numpy.in1d(map_pos, self._map_pos)
        # test whether map_fncode is in self._map_fncode
        if map_fncode is not None:
            masks[2] = numpy.in1d(map_fncode, self._map_fncode)

        # filter out Nones
        masks = tuple(filter(None, masks))

        # default value of None (for cases where len(masks) == 0)
        mask = None

        # if the len(masks) > 0, logical_and merge them together
        if len(masks) > 0:
            mask = numpy.logical_and.reduce(tuple(filter(None, masks)))

        # return mask
        return mask

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def from_egmap(fpath, mapfn = 'haldane'):
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

        genetic_map = GeneticMap.from_pandas_df(
            pandas_df = df,
            chr_grp_ix = 0,
            chr_start_ix = 1,
            chr_stop_ix = 2,
            map_pos_ix = 3,
            mkr_name_ix = 4 if "mkr_name" in df.columns else None,
            map_fncode_ix = 5 if "map_fncode" in df.columns else None,
            mapfn = mapfn
        )

        return genetic_map

    @staticmethod
    def from_pandas_df(pandas_df, chr_grp_ix = 0, chr_start_ix = 1,
                       chr_stop_ix = 2, map_pos_ix = 3, mkr_name_ix = None,
                       map_fncode_ix = None, mapfn = 'haldane'):
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
            2) chr_start (REQUIRED)
                The start position of the feature on the chromosome or scaffold.
                This is 0-indexed (e.g. the first base of a chromosome is 0) and
                inclusive (e.g. chr_start <= sequence < chr_stop).
                This is an integer type.
            3) chr_stop (REQUIRED)
                The stop position of the feature on the chromosome or scaffold. This
                is 0-indexed (e.g. the first base of a chromosome is 0) and
                exclusive (e.g. chr_start <= sequence < chr_stop).
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
        chr_grp_ix : int, default=0
            Column index to specify 'chrom' (chromosome) field.
        chr_start_ix : int, default=1
            Column index to specify 'chr_start' (chromosome start) field.
        chr_stop_ix : int, default=2
            Column index to specify 'chr_stop' (chromosome end) field.
        map_pos_ix : int, default=3
            Column index to specify 'map_pos' (map position in Morgans) field.
        mkr_name_ix : int, default=None
            Column index to specify 'mkr_name' (marker name) field.
        map_fncode_ix : int, default=None
            Column index to specify 'map_fncode' (mapping function code) field.

        Returns
        -------
        genetic_map : gmap
            A gmap object containing all data required for genetic inferences.
        """
        # check to make sure several indices aren't None.
        pybropt.util.check_is_not_none(chr_grp_ix, "chr_grp_ix")
        pybropt.util.check_is_not_none(chr_start_ix, "chr_start_ix")
        pybropt.util.check_is_not_none(chr_stop_ix, "chr_stop_ix")
        pybropt.util.check_is_not_none(map_pos_ix, "map_pos_ix")

        # get data
        chr_grp     = pandas_df.iloc[:,chr_grp_ix].values
        chr_start   = pandas_df.iloc[:,chr_start_ix].values
        chr_stop    = pandas_df.iloc[:,chr_stop_ix].values
        map_pos     = pandas_df.iloc[:,map_pos_ix].values
        mkr_name    = pandas_df.iloc[:,mkr_name_ix].values if mkr_name_ix is not None else None
        map_fncode  = pandas_df.iloc[:,map_fncode_ix].values if map_fncode_ix is not None else None

        # force data conversion
        chr_grp = numpy.string_([str(e) for e in chr_grp])
        chr_start = numpy.int64(chr_start)
        chr_stop = numpy.int64(chr_stop)
        map_pos = numpy.float64(map_pos)
        mkr_name = numpy.string_([str(e) for e in mkr_name]) if mkr_name is not None else None
        map_fncode = numpy.string_([str(e) for e in map_fncode]) if map_fncode is not None else None

        # construct the gmap object
        genetic_map = GeneticMap(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            map_pos = map_pos,
            mkr_name = mkr_name,
            map_fncode = map_fncode,
            mapfn = mapfn
        )

        return genetic_map

    @staticmethod
    def from_csv(fpath, sep = ',', header=0, chr_grp_ix = 0, chr_start_ix = 1,
            chr_stop_ix = 2, map_pos_ix = 3, mkr_name_ix = None,
            map_fncode_ix = None, mapfn = 'haldane'):
        """
        Parameters
        ----------
        fpath : str, path object, or file-like object
            Any valid string path, including URLs. Valid URL schemes include http,
            ftp, s3, and file. For file URLs, a host is expected. (see pandas docs)
        sep : str, default ','
            Delimiter to use.
        header : int, list of int, default=0
            Row number(s) to use as the column names, and the start of the data.

        """
        # read file using pandas.read_csv
        df = pandas.read_csv(
            fpath,
            sep=sep,
            header=header
        )

        genetic_map = GeneticMap.from_pandas_df(
            pandas_df = df,
            chr_grp_ix = chr_grp_ix,
            chr_start_ix = chr_start_ix,
            chr_stop_ix = chr_stop_ix,
            map_pos_ix = map_pos_ix,
            mkr_name_ix = mkr_name_ix,
            map_fncode_ix = map_fncode_ix,
            mapfn = mapfn
        )

        return genetic_map

    @staticmethod
    def cM2d(cM):
        """
        Convert centiMorgan units to genetic map distance units (equivalent to
        Morgan units).

        Parameters
        ----------
        cM : numpy.ndarray
            An array of centiMorgan map distance values.

        Returns
        -------
        d : numpy.ndarray
            An array of genetic map distance units (Morgans).
        """
        d = 0.01 * cM   # convert cM to d
        return d

    @staticmethod
    def haldane_to_r(d):
        """
        Convert genetic distances in Morgans to recombination probabilities using
        the Haldane mapping function (Haldane, 1919).

        This is a bare-bones function and does no parameter checking for errors.

        Parameters
        ----------
        d : numpy.ndarray
            An array of genetic distances in Morgans. This can be an array of any
            shape.

        Returns
        -------
        r : numpy.ndarray
            An array of recombination probabilities. The shape of the array is
            the same shape as that of 'd'.
        """
        # convert d to r
        r = 0.5 * (1.0 - numpy.exp(-2.0 * d))
        return r

    @staticmethod
    def r_to_haldane(r):
        """
        Convert recombination probabilities between two loci to genetic map
        distances using the Haldane mapping function (Haldane, 1919).

        This is a bare-bones function and does no parameter checking for errors.

        Parameters
        ----------
        r : numpy.ndarray
            An array of recombination probabilities between two loci.

        Returns
        -------
        d : numpy.ndarray
            An array of genetic map distances as defined by the Haldane mapping
            function.
        """
        # convert r to d
        d = -0.5 * numpy.log(1.0 - (2.0 * r))
        return d

    @staticmethod
    def kosambi_to_r(d):
        """
        Convert genetic map distances to recombination probabilities using the
        Kosambi mapping function (Kosambi, 1944).

        Parameters
        ----------
        d : numpy.ndarray
            An array of genetic map distances to convert to recombination
            probabilities.

        Returns
        -------
        r : numpy.ndarray
            An array of recombination probabilities between loci as defined by 'd'.
        """
        # convert d to r
        r = 0.5 * numpy.tanh(2.0 * d)
        return r

    @staticmethod
    def r_to_kosambi(r):
        """
        Convert recombination probabilities between two loci to genetic map
        distances using the Kosambi mapping function (Kosambi, 1944).

        Parameters
        ----------
        r : numpy.ndarray
            An array of recombination probabilities between two loci.

        Returns
        -------
        d : numpy.ndarray
            An array of genetic map distances as defined by the Kosambi mapping
            function.
        """
        # convert r to d
        d = numpy.log(1.0 + (2.0 * r)) / (4.0 - (8.0 * r))
        return d

    @staticmethod
    def key_to_mapfn(key):
        return key if callable(key) else GeneticMap.KEY_TO_MAPFN[key]

################################################################################
############################# Bind Class Constants #############################
################################################################################

# table of map function codes
GeneticMap.KEY_TO_MAP_FNCODE = {
    # haldane functions
    "haldane": "h",
    GeneticMap.haldane_to_r: 'h',
    GeneticMap.haldane_to_r.__name__: 'h',

    # kosambi functions
    "kosambi": "k",
    GeneticMap.kosambi_to_r: 'k',
    GeneticMap.kosambi_to_r.__name__: 'k',

    # other keys
    "unknown": "u",
    "None": 'u',
    None: 'u',
    "interpolate": "i"
}

GeneticMap.MAP_FNCODE_TO_KEY = {
    "h": "haldane",
    "i": "interpolate",
    "k": "kosambi",
    "u": "unknown"
}

GeneticMap.KEY_TO_MAPFN = {
    # haldane keys
    "h": GeneticMap.haldane_to_r,
    "haldane" : GeneticMap.haldane_to_r,

    # kosambi keys
    "k": GeneticMap.kosambi_to_r,
    "kosambi": GeneticMap.kosambi_to_r
}