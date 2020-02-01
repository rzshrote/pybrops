import numpy
import pandas


class gmap:
    """
    A class to represent genetic maps and perform associated operations.
    """
    #################
    # Instance Data #
    #################
    chrom = None        # chromosome names (align to chromSize)
    chromSize = None    # chromosome sizes (for indexing/iterating)
    chromStart = None   # marker start positions (inclusive)
    chromEnd = None     # marker stop positions (exclusive)
    mapPos = None       # map positions for each marker
    name = None         # name for each marker
    mapFn = None        # mapping function code for each marker

    def __init__(self, chrom, chromSize, chromStart, chromEnd, mapPos,
                 name = None, mapFn = None):
        """
        Constructor for a gmap object. The gmap class represents genetic maps.
        """
        # check for None, raise ValueError if needed.
        if chrom is None:
            raise ValueError("'chrom' cannot be None.")
        if chromSize is None:
            raise ValueError("'chromSize' cannot be None.")
        if chromStart is None:
            raise ValueError("'chromStart' cannot be None.")
        if chromEnd is None:
            raise ValueError("'chromEnd' cannot be None.")
        if mapPos is None:
            raise ValueError("'mapPos' cannot be None.")

        # set values
        self.chrom = chrom
        self.chromSize = chromSize
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.mapPos = mapPos
        self.name = name
        self.mapFn = mapFn

    def chrom_iter(chrom_name):
        

def read_gmap(fpath):
    """
    Read genetic map file (*.gmap).

    Assumes that the file is correctly formatted.

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
        2) chromStart (REQUIRED)
            The start position of the feature on the chromosome or scaffold.
            This is 0-indexed (e.g. the first base of a chromosome is 0) and
            inclusive (e.g. chromStart <= sequence < chromEnd).
            This is an integer type.
        3) chromEnd (REQUIRED)
            The stop position of the feature on the chromosome or scaffold. This
            is 0-indexed (e.g. the first base of a chromosome is 0) and
            exclusive (e.g. chromStart <= sequence < chromEnd).
            This is an integer type.
        4) mapPos (REQUIRED)
            The genetic map position in Morgans. (NOT centiMorgans!)
            This is an floating type.
        5) name (optional)
            The name of the marker on the genetic map.
            This is of type 'str'.
        6) mapFn (optional)
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
    chrom, chromSize = zip(*[(str(n),len(g)) for n,g in df.groupby(0)])

    # get values as numpy
    chrom_arr = numpy.object_(chrom)        # convert chrom to object ndarray
    chromSize = numpy.int64(chromSize)      # convert chromSize to int64 ndarray
    chromStart_arr = numpy.int64(df[1].values)  # get chromStart values
    chromEnd_arr = numpy.int64(df[2].values)    # get chromEnd values
    mapPos_arr = numpy.float64(df[3].values)    # get mapPos values

    # get name values if they exist
    name_arr = numpy.object_(df[4].values) if 4 in df else None
    # get mapFn values if they exist
    mapFn_arr = numpy.object_(df[5].values) if 5 in df else None

    # construct the gmap object
    genetic_map = gmap(
        chrom = chrom_arr,
        chromSize = chromSize_arr,
        chromStart = chromStart_arr,
        chromEnd = chromEnd_arr,
        mapPos = mapPos_arr,
        name = name_arr,
        mapFn = mapFn_arr
    )

    return genetic_map

def read_csv(fpath, sep = ',', header=0,
             chrom_ix = 0, chromStart_ix = 1, chromEnd_ix = 2, mapPos_ix = 3,
             name_ix = None, mapFn_ix = None):
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
        2) chromStart (REQUIRED)
            The start position of the feature on the chromosome or scaffold.
            This is 0-indexed (e.g. the first base of a chromosome is 0) and
            inclusive (e.g. chromStart <= sequence < chromEnd).
            This is an integer type.
        3) chromEnd (REQUIRED)
            The stop position of the feature on the chromosome or scaffold. This
            is 0-indexed (e.g. the first base of a chromosome is 0) and
            exclusive (e.g. chromStart <= sequence < chromEnd).
            This is an integer type.
        4) mapPos (REQUIRED)
            The genetic map position in Morgans. (NOT centiMorgans!)
            This is an floating type.
        5) name (optional)
            The name of the marker on the genetic map.
            This is of type 'str'.
        6) mapFn (optional)
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
    sep : str, default ','
        Delimiter to use.
    header : int, list of int, default=0
        Row number(s) to use as the column names, and the start of the data.
    chrom_ix : int, default=0
        Column index to specify 'chrom' (chromosome) field.
    chromStart_ix : int, default=1
        Column index to specify 'chromStart' (chromosome start) field.
    chromEnd_ix : int, default=2
        Column index to specify 'chromEnd' (chromosome end) field.
    mapPos_ix : int, default=3
        Column index to specify 'mapPos' (map position in Morgans) field.
    name_ix : int, default=None
        Column index to specify 'name' (marker name) field.
    mapFn_ix : int, default=None
        Column index to specify 'mapFn' (mapping function code) field.

    Returns
    =======
    genetic_map : gmap
        A gmap object containing all data required for genetic inferences.
    """
    # check to make sure several indices aren't None.
    if chrom_ix is None:
        raise ValueError("chrom_ix cannot be None.")
    if chromStart_ix is None:
        raise ValueError("chromStart_ix cannot be None.")
    if chromEnd_ix is None:
        raise ValueError("chromEnd_ix cannot be None.")
    if mapPos_ix is None:
        raise ValueError("mapPos_ix cannot be None.")

    # read file using pandas.read_csv
    df = pandas.read_csv(fpath, sep=sep, header=header)

    # build up a new DataFrame index and its corresponding integer indices
    new_Index = ["chrom", "chromStart", "chromEnd", "mapPos"]
    df_index = [chrom_ix, chromStart_ix, chromEnd_ix, mapPos_ix]
    if name_ix is not None:         # if name_ix was provided
        new_Index += ["name"]       # add string to indices
        df_index += [name_ix]       # add int to indices
    if mapFn_ix is not None:        # if mapFn_ix was provided
        new_Index += ["mapFn"]      # add string to indices
        df_index += [mapFn_ix]      # add int to indices

    # subset indices in the correct order
    df = df.iloc[:,df_index]

    # assign string indices to the columns for grouping (below)
    df.columns = pandas.Index(new_Index)

    # declare several variables
    chrom = []          # list to store string variables
    chromSize = []      # list to store linkage group sizes
    chromStart = []     # list to store chromStart's
    chromEnd = []       # list to store chromEnd's
    mapPos = []         # list to store mapPos's
    name = []           # list to store name's
    mapFn = []          # list to store mapFn codes

    # iterate though each group of markers by chromosome
    for chrom_name, chrom_df in df.groupby("chrom"):
        # sort the data inplace
        chrom_df.sort_values(
            by=["chromStart","chromEnd"],   # by chromStart, then by chromEnd
            ascending=[True, True],         # ascending for both columns
            inplace=True                    # modify dataframe inplace
        )

        #####################
        # add data to lists #
        #####################
        chrom.append(chrom_name)        # append chrom_name to chrom list
        chromSize.append(len(chrom_df)) # add linkage group length

        # extract data as numpy.ndarray's and add to lists
        chromStart.append(chrom_df["chromStart"].values)    # append chromStart
        chromEnd.append(chrom_df["chromEnd"].values)        # append chromEnd
        mapPos.append(chrom_df["mapPos"].values)            # append mapPos

        # conditional data fields
        if name_ix is not None:                     # if name given
            name.append(chrom_df["name"].values)    # append name
        if mapFn_ix is not None:                    # if mapFn given
            name.append(chrom_df["mapFn"].values)   # append mapFn

    ###################################
    # convert lists to numpy.ndarrays #
    ###################################
    chrom_arr = numpy.array(        # convert chrom list to numpy.ndarray
        [str(e) for e in chrom],    # force conversion to string
        dtype='object'              # store as object, not string type
    )

    chromSize_arr = numpy.int64(    # convert chromSize list to numpy.ndarray
        chromSize                   # chromSize is a list
    )

    chromStart_arr = numpy.int64(   # convert chromStart list to numpy.ndarray
        numpy.concatenate(          # concatenate list of numpy.ndarrays
            chromStart              # chromStart is a list of numpy.ndarrays
        )
    )

    chromEnd_arr = numpy.int64(     # convert chromEnd list to numpy.ndarray
        numpy.concatenate(          # concatenate list of numpy.ndarrays
            chromEnd                # chromEnd is a list of numpy.ndarrays
        )
    )

    mapPos_arr = numpy.float64(     # convert mapPos list to numpy.ndarray
        numpy.concatenate(          # concatenate list of numpy.ndarrays
            mapPos                  # mapPos is a list of numpy.ndarrays
        )
    )

    # conditional conversions
    name_arr = None                     # declare name array
    if name_ix is not None:             # if name is provided
        name_arr = numpy.concatenate(   # concatenate list of numpy.ndarrays
            name                        # name is a list of object (str) arrays
        )

    mapFn_arr = None                    # declare mapFn array
    if mapFn_ix is not None:            # if name is provided
        mapFn_arr = numpy.concatenate(  # concatenate list of numpy.ndarrays
            mapFn                       # mapFn is a list of object (str) arrays
        )

    # construct the gmap object
    genetic_map = gmap(
        chrom = chrom_arr,
        chromSize = chromSize_arr,
        chromStart = chromStart_arr,
        chromEnd = chromEnd_arr,
        mapPos = mapPos_arr,
        name = name_arr,
        mapFn = mapFn_arr
    )

    return genetic_map
