import numpy
import pandas
from scipy.interpolate import interp1d
from .const import GMAP_MAP2R_DICT

class gmap:
    """
    A class to represent genetic maps and perform associated operations.
    """
    #################
    # Instance Data #
    #################
    # new Instance data
    col_present = None  # tuple of booleans to indicate if column is present
    chr_grpls = None    # chromosome group name list (align to chr_grplen)
    chr_grplen = None   # chromosome group sizes (for indexing/iterating)
    chr_grpstix = None  # chromosome group start indices
    chr_grpspix = None  # chromosome group stop indices
    chr_grp = None      # chromosome group array
    chr_start = None    # marker start position (inclusive)
    chr_stop = None     # marker stop position (exclusive)
    map_pos = None      # map position for each marker
    map_name = None     # marker name
    map_fncode = None   # mapping function code for each marker
    map_spline = None   # genetic map spline for estimating genetic distances
    map_spline_t = None # genetic map spline type

    @classmethod
    def __init__(self, chr_grpls, chr_grplen, chr_grp, chr_start, chr_stop,
                 map_pos, map_name = None, map_fncode = None):
        """
        Constructor for a gmap object. The gmap class represents genetic maps.

        Parameters
        ==========
        chr_grpls : numpy.ndarray
            A numpy.ndarray of string objects (dtype="object") representing the
            names of the linkage groups.
            Example:
                chr_grpls = numpy.array(["chr1", "chr2", "chr3"])
        chr_grplen : numpy.ndarray
            A numpy.ndarray of int64 holding the number of markers
        """
        # check for None, raise ValueError if needed.
        if chr_grpls is None:
            raise ValueError("'chr_grpls' cannot be None.")
        if chr_grplen is None:
            raise ValueError("'chr_grplen' cannot be None.")
        if chr_grp is None:
            raise ValueError("'chr_grp' cannot be None.")
        if chr_start is None:
            raise ValueError("'chr_start' cannot be None.")
        if chr_stop is None:
            raise ValueError("'chr_stop' cannot be None.")
        if map_pos is None:
            raise ValueError("'map_pos' cannot be None.")

        # set values
        self.col_present = tuple(arg is not None for arg in (
                chr_grp, chr_start, chr_stop, map_pos, map_name, map_fncode
            )
        )
        self.chr_grpls = chr_grpls
        self.chr_grplen = chr_grplen
        self.chr_grpstix = numpy.cumsum(chr_grplen) - chr_grplen[0]
        self.chr_grpspix = numpy.cumsum(chr_grplen)
        self.chr_grp = chr_grp
        self.chr_start = chr_start
        self.chr_stop = chr_stop
        self.map_pos = map_pos
        self.map_name = map_name
        self.map_fncode = map_fncode

    @classmethod
    def __call__(self, chr = None, start = None, stop = None, mstart = None,
                 mstop = None, nstart = None, nstop = None, elements = False):
        """
        Return an iterator of indices for accessing specific ranges of the
        genetic map.

        Parameters
        ==========
        chr : None, str
            Chromosome name string. If not provided, iterate over the entire
            genetic map.
        start : None, int
            Chromosome start position (inclusive) in nucleotides.
            ***Must specify 'chr' in addition to this argument.***
        stop : None, int
            Chromosome stop position (exclusive) in nucleotides.
            ***Must specify 'chr' in addition to this argument.***
        mstart : None, int
            Genetic distance start position (inclusive) in Morgans.
            ***Must specify 'chr' in addition to this argument.***
        mstop : None, int
            Genetic distance stop position (exclusive) in Morgans.
            ***Must specify 'chr' in addition to this argument.***
        nstart : None, int
            Begin iterating (inclusive) at a starting marker name.
            ***Must have provided names with the genetic map.***
        nstop : None, int
            Stop iterating (inclusive) at a starting marker name.
            ***Must have provided names with the genetic map.***
        elements : boolean, default=False
            If False, return indices for accessing elements in gmap.
            If True, return the 'chrom', 'chr_start', 'chr_stop', 'map_pos',
            'map_name', and 'map_fncode' for the specified interval.

        Returns
        =======
        mapiter : iter
            An iterator object to iterate over the gmap.
        """
        # declare iter start and stop positions; assigned in logic tree below
        iter_start = 0
        iter_stop = len(self.map_pos)

        # if either nstart or nstop have been provided
        if (nstart is not None) or (nstop is not None):
            # if any additional arguments have been passed, raise an error.
            if any(arg is not None for arg in (start, stop, mstart, mstop)):
                raise TypeError(
                    "'nstart' and 'nstop' options are incompatible with "\
                    "'start', 'stop', 'mstart', and 'mstop' options."
                )

            # find the positions of matching strings
            nstart_pos = numpy.flatnonzero(self.map_name == nstart)
            nstop_pos = numpy.flatnonzero(self.map_name == nstop)

            # if there are multiple matching strings for either, raise error
            if (len(nstart_pos) > 1) or (len(nstop_pos) > 1):
                raise KeyError(
                    "Multiple keys found for 'nstart' and 'nstop'.\n"\
                    "    'nstart' keys found: %d\n    'nstop' keys found: %d\n"
                    % (len(nstart_pos), len(nstop_pos))
                )

            # if both the start and stop positions were not found, raise error
            if (len(nstart_pos) == 0) and (len(nstop_pos) == 0):
                raise KeyError(
                    "Unable to find keys for 'nstart' and 'nstop'.\n"\
                    "    'nstart' keys found: %d\n    'nstop' keys found: %d\n"
                    % (len(nstart_pos), len(nstop_pos))
                )

            # calculate the start and stop indices
            iter_start = nstart_pos[0] if len(nstart_pos) == 1 else 0
            iter_stop = nstop_pos[0] if len(nstop_pos) == 1 else len(self.map_pos)

        # else if the chrom is specified (must be specified for other args)
        elif chrom is not None:
            # find the 'chrom' in self.chrom
            chrom_pos = numpy.flatnonzero(self.chrom == chrom)

            # if we do not have a match, raise an error
            if len(chrom_pos) != 1:
                raise KeyError("Unable to find key: %s" % chrom)

            # gather chrom view indices
            csp = self.chr_grplen[:chrom_pos[0]+1].sum() # chrom start index
            cst = csp - self.chr_grplen[chrom_pos[0]]    # chrom stop index

            # find the start and stop positions
            start_pos = self.chr_start[cst:csp].searchsorted(start) if start is not None else cst
            stop_pos = self.chr_start[cst:csp].searchsorted(stop) if stop is not None else csp
            mstart_pos = self.map_pos[cst:csp].searchsorted(mstart) if mstart is not None else cst
            mstop_pos = self.map_pos[cst:csp].searchsorted(mstop) if mstop is not None else csp

            # assign the union of the positions
            iter_start = start_pos if start_pos < mstart_pos else mstart_pos
            iter_stop = stop_pos if stop_pos > mstop_pos else mstop_pos

        # else if any args are not None raise an error b
        elif any(arg is not None for arg in (start, stop, mstart, mstop)):
                raise ValueError("'chrom' cannot be None.")

        # build the iterator; if elements == True, return elements
        it = None
        if elements:
            it = ((chrom, chr_start[i], chr_stop[i], map_pos[i], map_name[i],
            map_fncode[i]) for i in range(iter_start, iter_stop))
        else:
            it = iter(range(iter_start, iter_stop))

        return it

    @classmethod
    def __len__(self):
        return len(self.map_pos)

    @classmethod
    def dist(self, rst, rsp, cst, csp):
        """
        Get distance matrix.
        """
        # make an mesh grid for row and column gmap's
        mesh = numpy.meshgrid(
            self.map_pos[rst:rsp],
            self.map_pos[cst:csp],
            indexing='ij'
        )

        # take the difference of the two components in the grid and return abs()
        d = numpy.absolute(mesh[0] - mesh[1])

        # return distances
        return d


    @classmethod
    def build_spline(self, kind = 'linear', fill_value = 'extrapolate'):
        """
        Build a spline for estimating genetic map distances. This is built
        using the marker start indices (self.chr_start)

        Parameters
        ==========
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
        # allocate empty object array for splines
        self.map_spline = numpy.empty(
            len(self.chr_grplen),
            dtype = 'object'
        )

        # set the spline type
        self.map_spline_t = kind

        # for each chromosome, build its spline
        for i,(st,sp) in enumerate(zip(self.chr_grpstix, self.chr_grpspix)):
            # stash a spline object into self.map_spline
            self.map_spline[i] = interp1d(
                x = self.chr_start[st:sp],  # chromsome positions
                y = self.map_pos[st:sp],    # map positions
                kind = kind,                # type of interpolation
                fill_value = fill_value,    # default fill value
                assume_sorted = True        # assume that the array is sorted
            )

    @classmethod
    def to_pandas_df(self, header = False):
        """
        Convert genetic map to a pandas.DataFrame.

        Parameters
        ==========
        header : boolean, default = True
            Include non-numeric indices.

        Returns
        =======
        df : pandas.DataFrame
        """
        # declare header dictionary variable
        df_dict = None

        # if we want a named header
        if header:
            # fill with required fields
            df_dict = {
                "chr_grp" : self.chr_grp if self.col_present[0] else None,
                "chr_start" : self.chr_start if self.col_present[1] else None,
                "chr_stop" : self.chr_stop if self.col_present[2] else None,
                "map_pos" : self.map_pos if self.col_present[3] else None,
                "map_name" : self.map_name if self.col_present[4] else None,
                "map_fncode" : self.map_fncode if self.col_present[5] else None
            }

        # else we do not want a named header (only integer)
        else:
            # fill with required fields
            df_dict = {
                0 : self.chr_grp if self.col_present[0] else None,
                1 : self.chr_start if self.col_present[1] else None,
                2 : self.chr_stop if self.col_present[2] else None,
                3 : self.map_pos if self.col_present[3] else None,
                4 : self.map_name if self.col_present[4] else None,
                5 : self.map_fncode if self.col_present[5] else None
            }

        # create DataFrame
        df = pandas.DataFrame(df_dict)

        # return DataFrame
        return df

    @classmethod
    def to_gmap(self, fname):
        """
        Convert gmap to *.gmap file.
        """
        df = self.to_pandas_df()

        df.to_csv(fname, sep='\t', header=False, index=None)


    @classmethod
    def interpolate(self, pandas_df, chr_grp_ix = 0, chr_start_ix = 1,
                    chr_stop_ix = 2, map_name_ix = None, map_fncode_ix = None,
                    kind = 'linear', fill_value = 'extrapolate'):
        """
        Interpolate genetic map positions given vectors of coordinates.

        Parameters
        ==========
        pandas_df : pandas.DataFrame
            A pandas DataFrame containing columns for 'chrom', 'chr_start',
            'chr_stop', and optionally 'map_name' and 'map_fncode'
        chr_grp_ix : int, default = 0
            An index for the pandas_df column containing 'chrom' positions.
        chr_start_ix : int, default = 1
            An index for the pandas_df column containing 'chr_start' positions.
        chr_stop_ix : int, default = 2
            An index for the pandas_df column containing 'chr_stop' positions.
        map_name_ix : int, default = None
            An index for the pandas_df column containing names.
        map_fncode_ix : int, default = None
            An index for the pandas_df column containing 'map_fncode' codes.
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

        Returns
        =======
        genetic_map : gmap
            A genetic map containing interpolated positions
        """
        # check to make sure several indices aren't None.
        if chr_grp_ix is None:
            raise ValueError("chr_grp_ix cannot be None.")
        if chr_start_ix is None:
            raise ValueError("chr_start_ix cannot be None.")
        if chr_stop_ix is None:
            raise ValueError("chr_stop_ix cannot be None.")

        # get input indices as tuple
        col_ix = (
            chr_grp_ix,
            chr_start_ix,
            chr_stop_ix,
            True,       # set true because we will be calculating these
            map_name_ix,
            map_fncode_ix
        )

        # determine presence-absence of columns
        col_present = tuple(ix is not None for ix in col_ix)

        # construct DataFrame subset (put Nones where needed)
        df_dict = {
            0: pandas_df.iloc[:,chr_grp_ix].values if col_present[0] else None,
            1: pandas_df.iloc[:,chr_start_ix].values if col_present[1] else None,
            2: pandas_df.iloc[:,chr_stop_ix].values if col_present[2] else None,
            3: None,
            4: pandas_df.iloc[:,map_name_ix].values if col_present[4] else None,
            5: pandas_df.iloc[:,map_fncode_ix].values if col_present[5] else None
        }

        # create DataFrame
        df = pandas.DataFrame(df_dict)

        # if our splines are not present/correct, (re)build them
        if self.map_spline_t != kind:
            self.build_spline(kind = kind, fill_value = fill_value)

        # group the input by first column (chromosome)
        df_groups = df.groupby(0)

        # declare several variables
        new_chr_grpls = []  # list to store string variables
        new_chr_grplen = [] # list to store linkage group sizes
        new_chr_grp = []    # list to store chr_grp's
        new_chr_start = []  # list to store chr_start's
        new_chr_stop = []   # list to store chr_stop's
        new_map_pos = []    # list to store map_pos's
        new_map_name = []   # list to store map_name's
        new_map_fncode = [] # list to store map_fncode codes

        # for each chromosome in self,
        for i, chr_name in enumerate(self.chr_grpls):
            # if the input DataFrame has our chromosome group, interpolate
            if chr_name in df_groups.groups:
                # get the group and sort by chr_start and chr_stop
                df_group = df_groups.get_group(chr_name).sort_values([1,2])

                #####################
                # add data to lists #
                #####################
                chr_grpls.append(str(chr_name))                         # append chrom_name to chrom list
                chr_grplen.append(len(df_group))                        # add linkage group length
                chr_grp.append(numpy.object_(df_group[0].values))       # append chr_grp
                chr_start.append(numpy.int64(df_group[1].values))       # append chr_start
                chr_stop.append(numpy.int64(df_group[2].values))        # append chr_stop
                map_pos.append(numpy.float64(self.map_spline[i](        # append spline approximations
                    df_group[1].values
                )))
                map_name.append(numpy.object_(df_group[4].values))      # append map_names (even if None)
                map_fncode.append(numpy.object_(df_group[5].values))    # append map_fncode (even if None)

        ###################################
        # convert lists to numpy.ndarrays #
        ###################################
        chr_grpls = numpy.object_(chr_grpls)
        chr_grplen = numpy.int64(chr_grplen)
        chr_grp = numpy.concatenate(chr_grp)
        chr_start = numpy.concatenate(chr_start)
        chr_stop = numpy.concatenate(chr_stop)
        map_pos = numpy.concatenate(map_pos)
        map_name = numpy.concatenate(map_name)
        map_fncode = numpy.concatenate(map_fncode)

        # construct the gmap object
        genetic_map = gmap(
            chr_grpls   = chr_grpls,
            chr_grplen  = chr_grplen,
            chr_grp     = chr_grp if col_present[0] else None,
            chr_start   = chr_start if col_present[1] else None,
            chr_stop    = chr_stop if col_present[2] else None,
            map_pos     = map_pos if col_present[3] else None,
            map_name    = map_name if col_present[4] else None,
            map_fncode  = map_fncode if col_present[5] else None
        )

        return genetic_map

    @classmethod
    def remove_discrepancies(self):
        """
        Remove discrepancies between the physical map and the genetic map.
        In instances of conflict, assume that the physical map is correct.

        Note:
            This assumption may cause major issues if there are incorrect
            markers at the beginning of the chromosome.
        """
        # make empty bool vector
        correct_pos = numpy.zeros(len(self.chr_start), dtype = 'bool')

        new_chr_grplen = numpy.ones(len(self.chr_grplen), dtype = 'int64')

        # counter for correct_pos indexing
        k = 0

        # for each linkage group
        for i,(stix,spix) in enumerate(zip(self.chr_grpstix, self.chr_grpspix)):
            prev_pos = self.map_pos[stix]   # get the first map position
            correct_pos[k] = True           # we assume it is true
            for j in range(stix+1, spix):   # for each marker after
                k += 1
                if prev_pos <= self.map_pos[j]: # if the current is <= prev
                    correct_pos[k] = True       # the position is consistent
                    prev_pos = self.map_pos[j]  # move the previous position up
                    new_chr_grplen[i] += 1      # increase correct length by one
                else:
                    correct_pos[k] = False      # else, it is inconsistent

        # next we physically remove discrepancies
        self.chr_grplen = new_chr_grplen
        self.chr_grpstix = numpy.cumsum(new_chr_grplen) - new_chr_grplen[0]
        self.chr_grpspix = numpy.cumsum(new_chr_grplen)
        self.chr_grp = self.chr_grp[correct_pos]
        self.chr_start = self.chr_start[correct_pos]
        self.chr_stop = self.chr_stop[correct_pos]
        self.map_pos = self.map_pos[correct_pos]
        self.map_fncode = self.map_fncode[correct_pos]


    @staticmethod
    def read_gmap(fpath):
        """
        Read genetic map file (*.gmap).

        Assumes that the file is correctly formatted.

        =============================================================
        Genetic map file (*.gmap) format (similar to BED file format)
        =============================================================
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
            5) map_name (optional)
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
        ==========
        fpath : str, path object, or file-like object
            Any valid string path, including URLs. Valid URL schemes include http,
            ftp, s3, and file. For file URLs, a host is expected. (see pandas docs)

        Returns
        =======
        genetic_map : gmap
            A gmap object containing all data required for genetic inferences.
        """
        # read the file using pandas
        df = pandas.read_csv(fpath, sep='\t', header=None)

        # get lists of chromosome names and chromosome group lengths
        chr_grpls, chr_grplen = zip(*[(str(n),len(g)) for n,g in df.groupby(0)])

        # get values as numpy.ndarray's
        chr_grpls = numpy.object_(chr_grpls)    # convert chrom to object
        chr_grplen = numpy.int64(chr_grplen)    # convert chr_grplen to int64
        chr_grp = numpy.object_(df[0].values)   # get chr_grp values
        chr_start = numpy.int64(df[1].values)   # get chr_start values
        chr_stop = numpy.int64(df[2].values)    # get chr_stop values
        map_pos = numpy.float64(df[3].values)   # get map_pos values

        # get name values if they exist
        map_name = numpy.object_(df[4].values) if 4 in df else None
        # get map_fncode values if they exist
        map_fncode = numpy.object_(df[5].values) if 5 in df else None

        # construct the gmap object
        genetic_map = gmap(
            chr_grpls = chr_grpls,
            chr_grplen = chr_grplen,
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            map_pos = map_pos,
            map_name = map_name,
            map_fncode = map_fncode
        )

        return genetic_map

    @staticmethod
    def from_pandas_df(pandas_df, chr_grp_ix = 0, chr_start_ix = 1,
                       chr_stop_ix = 2, map_pos_ix = 3, map_name_ix = None,
                       map_fncode_ix = None):
        """
        Read genetic map data from a CSV-like file.

        =============================================================
        Genetic map file (*.gmap) format (similar to BED file format)
        =============================================================
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
            5) map_name (optional)
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
        ==========
        chr_grp_ix : int, default=0
            Column index to specify 'chrom' (chromosome) field.
        chr_start_ix : int, default=1
            Column index to specify 'chr_start' (chromosome start) field.
        chr_stop_ix : int, default=2
            Column index to specify 'chr_stop' (chromosome end) field.
        map_pos_ix : int, default=3
            Column index to specify 'map_pos' (map position in Morgans) field.
        map_name_ix : int, default=None
            Column index to specify 'map_name' (marker name) field.
        map_fncode_ix : int, default=None
            Column index to specify 'map_fncode' (mapping function code) field.

        Returns
        =======
        genetic_map : gmap
            A gmap object containing all data required for genetic inferences.
        """
        # check to make sure several indices aren't None.
        if chr_grp_ix is None:
            raise ValueError("chr_grp_ix cannot be None.")
        if chr_start_ix is None:
            raise ValueError("chr_start_ix cannot be None.")
        if chr_stop_ix is None:
            raise ValueError("chr_stop_ix cannot be None.")
        if map_pos_ix is None:
            raise ValueError("map_pos_ix cannot be None.")

        # get input indices as tuple
        col_ix = (
            chr_grp_ix,
            chr_start_ix,
            chr_stop_ix,
            map_pos_ix,
            map_name_ix,
            map_fncode_ix
        )

        # determine presence-absence of columns
        col_present = tuple(ix is not None for ix in col_ix)

        # construct DataFrame subset (put Nones where needed)
        df_dict = {
            0: pandas_df.iloc[:,chr_grp_ix].values if col_present[0] else None,
            1: pandas_df.iloc[:,chr_start_ix].values if col_present[1] else None,
            2: pandas_df.iloc[:,chr_stop_ix].values if col_present[2] else None,
            3: pandas_df.iloc[:,map_pos_ix].values if col_present[3] else None,
            4: pandas_df.iloc[:,map_name_ix].values if col_present[4] else None,
            5: pandas_df.iloc[:,map_fncode_ix].values if col_present[5] else None
        }

        # create DataFrame
        df = pandas.DataFrame(df_dict)

        # declare several variables
        chr_grpls = []  # list to store string variables
        chr_grplen = [] # list to store linkage group sizes
        chr_grp = []    # list to store chr_grp's
        chr_start = []  # list to store chr_start's
        chr_stop = []   # list to store chr_stop's
        map_pos = []    # list to store map_pos's
        map_name = []   # list to store map_name's
        map_fncode = [] # list to store map_fncode codes

        # iterate though each group of markers by chromosome (col index = 0)
        for chr_name, chr_df in df.groupby(0):
            # sort the data inplace
            chr_df = chr_df.sort_values(
                by=[1,2],                   # by chr_start, then by chr_stop
                ascending=[True, True]      # ascending for both columns
            )

            #####################
            # add data to lists #
            #####################
            chr_grpls.append(str(chr_name)) # append chrom_name to chrom list
            chr_grplen.append(len(chr_df))  # add linkage group length

            # extract data as numpy.ndarray's and add to lists
            chr_grp.append(numpy.object_(chr_df[0].values))     # append chr_grp
            chr_start.append(numpy.int64(chr_df[1].values))     # append chr_start
            chr_stop.append(numpy.int64(chr_df[2].values))      # append chr_stop
            map_pos.append(numpy.float64(chr_df[3].values))     # append map_pos
            map_name.append(numpy.object_(chr_df[4].values))    # append map_names (even if None)
            map_fncode.append(numpy.object_(chr_df[5].values))  # append map_fncode (even if None)

        ###################################
        # convert lists to numpy.ndarrays #
        ###################################
        chr_grpls = numpy.object_(chr_grpls)
        chr_grplen = numpy.int64(chr_grplen)
        chr_grp = numpy.concatenate(chr_grp)
        chr_start = numpy.concatenate(chr_start)
        chr_stop = numpy.concatenate(chr_stop)
        map_pos = numpy.concatenate(map_pos)
        map_name = numpy.concatenate(map_name)
        map_fncode = numpy.concatenate(map_fncode)

        # construct the gmap object
        genetic_map = gmap(
            chr_grpls   = chr_grpls,
            chr_grplen  = chr_grplen,
            chr_grp     = chr_grp if col_present[0] else None,
            chr_start   = chr_start if col_present[1] else None,
            chr_stop    = chr_stop if col_present[2] else None,
            map_pos     = map_pos if col_present[3] else None,
            map_name    = map_name if col_present[4] else None,
            map_fncode  = map_fncode if col_present[5] else None
        )

        return genetic_map

    @staticmethod
    def from_csv(fpath, sep = ',', header=0,
                 chr_grp_ix = 0, chr_start_ix = 1, chr_stop_ix = 2,
                 map_pos_ix = 3, map_name_ix = None, map_fncode_ix = None):
        """
        Parameters
        ==========
        fpath : str, path object, or file-like object
            Any valid string path, including URLs. Valid URL schemes include http,
            ftp, s3, and file. For file URLs, a host is expected. (see pandas docs)
        sep : str, default ','
            Delimiter to use.
        header : int, list of int, default=0
            Row number(s) to use as the column names, and the start of the data.

        """
        # read file using pandas.read_csv
        df = pandas.read_csv(fpath, sep=sep, header=header)

        genetic_map = gmap.from_pandas_df(
            pandas_df = df,
            chr_grp_ix = chr_grp_ix,
            chr_start_ix = chr_start_ix,
            chr_stop_ix = chr_stop_ix,
            map_pos_ix = map_pos_ix,
            map_name_ix = map_name_ix,
            map_fncode_ix = map_fncode_ix
        )

        return genetic_map

    @classmethod
    def to_r(self, mapfn = 'haldane', mtype = 'sym',
             dtype = numpy.dtype('float64')):
        """
        Convert genetic map positions to a matrix of recombination probabilities
        using a provided mapping function.

        Parameters
        ----------
        mapfn : callable, {'haldane', 'kosambi'}, default = 'haldane'
            Mapping function to use.
            Options:
                callable    : A custom callable Python function.
                              Must take the following arguments:
                                mapfn(d)
                                Parameters:
                                    d : numpy.ndarray
                                        An array of genetic distances in Morgans.
                                        This can be an array of any shape.
                                Returns:
                                    r : numpy.ndarray
                                        An array of recombination probabilities.
                                        The shape of the array is the same shape
                                        as that of 'd'.
                'haldane'   : Haldane (1919) mapping function
                'kosambi'   : Kosambi (1944) mapping function
        mtype : {'sym', 'triu', 'tril'}, default = 'sym'
            Matrix type to return.
            Options:
                'sym'   : Symmetric matrix
                'triu'  : Upper triangular matrix
                'tril'  : Lower triangular matrix
        dtype : numpy.dtype, default = numpy.dtype('float64')
            The type of the output array. If dtype is None, infer the data
            type from the 'gmap' argument.

        Returns
        -------
        r : numpy.ndarray
            An array of recombination probabilities between loci as defined by 'd'.
        """
        if not callable(mapfn):
            mapfn = GMAP_MAP2R_DICT[mapfn]

        # if the mtype is not in the list of options, raise an error
        if mtype not in ['sym', 'triu', 'tril']:
            raise ValueError(
                "Invalid 'mtype'. Options are:\n"\
                "    'sym'   : Symmetric matrix\n"\
                "    'triu'  : Upper triangular matrix\n"\
                "    'tril'  : Lower triangular matrix\n"
            )

        ########################### Compute whole matrix ###########################
        # allocate an empty matrix for recombination probabilities
        r = numpy.empty(            # make empty matrix
            (len(self.map_pos), len(self.map_pos)), # make a square matrix of width len(gmap)
            dtype=dtype             # set the data type as dtype
        )

        # for each start, stop index
        for st,sp in zip(self.chr_grpstix, self.chr_grpspix):
            # make mesh for that array region
            mesh = numpy.meshgrid(self.map_pos[st:sp], self.map_pos[st:sp])

            # calculate recombination for that region
            r[st:sp,st:sp] = mapfn(numpy.abs(mesh[0] - mesh[1]))

            # fill everything else assuming independent assortment
            r[:st,st:sp].fill(0.5)  # fill matrix above with 0.5
            r[sp:,st:sp].fill(0.5)  # fill matrix below with 0.5

        if mtype == 'triu':             # if 'triu', multiply by a triu mask
            r *= numpy.logical_not(     # invert triangle mask for triu (fast)
                numpy.tri(              # make inverse of triu mask
                    *r.shape[-2:],      # get first two dimensions
                    k=-1,               # offset to get lower triangle
                    dtype=numpy.bool    # use bool (1 byte) data type (fast)
                )
            )
        elif mtype == 'tril':           # if 'tril', multiply by a tril mask
            r *= numpy.tri(             # make lower triangle mask
                *r.shape[-2:],          # get first two dimensions
                k=0,                    # offset of zero to get lower triangle mask
                dtype=numpy.bool        # use bool (1 byte) data type (fast)
            )

        # return probability matrix
        return r
