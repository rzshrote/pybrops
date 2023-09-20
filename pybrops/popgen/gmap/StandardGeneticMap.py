"""
Module implementing a standard genetic map format and associated error checking
routines.
"""

__all__ = [
    "StandardGeneticMap",
    "check_is_StandardGeneticMap",
]

import copy
import math
from numbers import Integral
from typing import Optional, Sequence, Union
import warnings
import numpy
from scipy.interpolate import interp1d
from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_floating, check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_ndim
from pybrops.core.error.error_type_python import check_is_dict, check_is_str
from pybrops.popgen.gmap.GeneticMap import GeneticMap


class StandardGeneticMap(GeneticMap):
    """
    A concrete class for representing a standard genetic map format.

    The purpose of this concrete class is to implement functionality for:
        1) Genetic map representation.
        2) Genetic map metadata.
        3) Genetic map routines.
        4) Genetic map interpolation spline construction.
        5) Genetic map spline interpolation.
        6) Import and export of genetic maps.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            vrnt_genpos_units: str = "M",
            auto_group: bool = True,
            auto_build_spline = True,
            **kwargs: dict
        ) -> None:
        """
        Constructor for creating a standard genetic map object.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            Chromosome or linkage group assignment array of shape ``(n,)`` where 
            ``n`` is the number of markers.
        
        vrnt_phypos : numpy.ndarray
            Physical positions array of shape ``(n,)`` where ``n`` is the number 
            of markers. This array contains the physical positions on the 
            chromosome or linkage group for each marker.
        
        vrnt_genpos : numpy.ndarray
            Genetic positions array of shape ``(n,)`` where ``n`` is the number 
            of markers. This array contains the genetic positions on the 
            chromosome or linkage group for each marker.

        vrnt_genpos_units : str, default = "M"
            Units in which genetic positions in the ``vrnt_genpos`` array are 
            stored. Options are:
            
            - ``"M"`` - genetic position units are in Morgans
            - ``"Morgans"`` - genetic position units are in Morgans
            - ``"cM"`` - genetic position units are in centiMorgans
            - ``"centiMorgans"`` - genetic position units are in centiMorgans
            
            Internally, all genetic positions are stored in Morgans. Providing 
            the units of the input  

        auto_group : bool
            Whether to automatically sort and group variants into chromosome groups.
        
        auto_build_spline : bool
            Whether to automatically construct a 

        kwargs : dict
            Additional keyword arguments.
        """
        # convert genetic positions to Morgans
        units = vrnt_genpos_units.lower()
        if units in ("m","morgans"):
            pass
        elif units in ("cm","centimorgans"):
            vrnt_genpos = 0.01 * vrnt_genpos
        else:
            raise ValueError(
                "argument '{0}' must be of value '{1}' or '{2}' but received value '{3}'".format(
                    "vrnt_genpos_units",
                    "Morgans",
                    "centiMorgans",
                    vrnt_genpos_units
                )
            )

        # order dependent assignments!
        self.vrnt_chrgrp = vrnt_chrgrp # must be assigned first!
        self.vrnt_phypos = vrnt_phypos # depends on vrnt_chrgrp assignment
        self.vrnt_genpos = vrnt_genpos # depends on vrnt_chrgrp assignment


    def __len__(self):
        """Get the number of markers in the genetic map."""
        return len(self._vrnt_genpos)

    ################## GeneticMap copying ##################
    def __copy__(
            self
        ) -> 'StandardGeneticMap':
        """
        Make a shallow copy of the StandardGeneticMap.

        Returns
        -------
        out : StandardGeneticMap
            A shallow copy of the original StandardGeneticMap.
        """
        # get the class of self
        cls = self.__class__

        # create a new object instance
        out = cls.__new__(cls)

        # initialize the object instance
        out.__init__(
            vrnt_chrgrp = copy.copy(self.vrnt_chrgrp), 
            vrnt_phypos = copy.copy(self.vrnt_phypos), 
            vrnt_genpos = copy.copy(self.vrnt_genpos),
        )

        # copy variant metadata to the new object

        # copy the spline metadata to the new object

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'StandardGeneticMap':
        """
        Make a deep copy of the StandardGeneticMap.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : StandardGeneticMap
            A deep copy of the original StandardGeneticMap.
        """
        # get the class of self
        cls = self.__class__

        # create a new object instance
        out = cls.__new__(cls)

        # initialize the object instance
        out.__init__(
            vrnt_chrgrp = copy.deepcopy(self.vrnt_chrgrp), 
            vrnt_phypos = copy.deepcopy(self.vrnt_phypos), 
            vrnt_genpos = copy.deepcopy(self.vrnt_genpos),
        )

        # copy variant metadata to the new object

        # copy the spline metadata to the new object

        return out

    ############################ Object Properties #############################
    @property
    def nvrnt(self) -> Integral:
        """Number of variants in the GeneticMap."""
        return len(self._vrnt_genpos)

    @property
    def vrnt_chrgrp(self) -> numpy.ndarray:
        """Variant chromosome group label."""
        return self._vrnt_chrgrp
    @vrnt_chrgrp.setter
    def vrnt_chrgrp(self, value: numpy.ndarray) -> None:
        """Set variant chromosome group label array"""
        check_is_ndarray(value, "vrnt_chrgrp")
        check_ndarray_dtype_is_integer(value, "vrnt_chrgrp")
        check_ndarray_ndim(value, "vrnt_chrgrp", 1)
        self._vrnt_chrgrp = value

    @property
    def vrnt_phypos(self) -> numpy.ndarray:
        """Variant physical position."""
        return self._vrnt_phypos
    @vrnt_phypos.setter
    def vrnt_phypos(self, value: numpy.ndarray) -> None:
        """Set variant physical position array"""
        check_is_ndarray(value, "vrnt_phypos")
        check_ndarray_dtype_is_integer(value, "vrnt_phypos")
        check_ndarray_ndim(value, "vrnt_phypos", 1)
        check_ndarray_len_eq(value, "vrnt_phypos", len(self.vrnt_chrgrp))
        self._vrnt_phypos = value

    @property
    def vrnt_genpos(self) -> numpy.ndarray:
        """Variant genetic position in Morgans."""
        return self._vrnt_genpos
    @vrnt_genpos.setter
    def vrnt_genpos(self, value: numpy.ndarray) -> None:
        """Set variant genetic position array"""
        check_is_ndarray(value, "vrnt_genpos")
        check_ndarray_dtype_is_floating(value, "vrnt_genpos")
        check_ndarray_ndim(value, "vrnt_genpos", 1)
        check_ndarray_len_eq(value, "vrnt_genpos", len(self.vrnt_chrgrp))
        self._vrnt_genpos = value

    ################# Metadata Properites ##################
    @property
    def vrnt_chrgrp_name(self) -> Union[numpy.ndarray,None]:
        """Variant chromosome group names."""
        return self._vrnt_chrgrp_name
    @vrnt_chrgrp_name.setter
    def vrnt_chrgrp_name(self, value: Union[numpy.ndarray,None]) -> None:
        """Set variant chromosome group name array"""
        if value is not None:
            check_is_ndarray(value, "vrnt_chrgrp_name")
            check_ndarray_dtype_is_integer(value, "vrnt_chrgrp_name")
            check_ndarray_ndim(value, "vrnt_chrgrp_name", 1)
        self._vrnt_chrgrp_name = value

    @property
    def vrnt_chrgrp_stix(self) -> Union[numpy.ndarray,None]:
        """Variant chromosome group start indices."""
        return self._vrnt_chrgrp_stix
    @vrnt_chrgrp_stix.setter
    def vrnt_chrgrp_stix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set variant chromosome group start indices array"""
        if value is not None:
            check_is_ndarray(value, "vrnt_chrgrp_stix")
            check_ndarray_dtype_is_integer(value, "vrnt_chrgrp_stix")
            check_ndarray_ndim(value, "vrnt_chrgrp_stix", 1)
        self._vrnt_chrgrp_stix = value

    @property
    def vrnt_chrgrp_spix(self) -> Union[numpy.ndarray,None]:
        """Variant chromosome group stop indices."""
        return self._vrnt_chrgrp_spix
    @vrnt_chrgrp_spix.setter
    def vrnt_chrgrp_spix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set variant chromosome group stop indices array"""
        if value is not None:
            check_is_ndarray(value, "vrnt_chrgrp_spix")
            check_ndarray_dtype_is_integer(value, "vrnt_chrgrp_spix")
            check_ndarray_ndim(value, "vrnt_chrgrp_spix", 1)
        self._vrnt_chrgrp_spix = value

    @property
    def vrnt_chrgrp_len(self) -> Union[numpy.ndarray,None]:
        """Variant chromosome group length."""
        return self._vrnt_chrgrp_len
    @vrnt_chrgrp_len.setter
    def vrnt_chrgrp_len(self, value: Union[numpy.ndarray,None]) -> None:
        """Set variant chromosome group length array"""
        if value is not None:
            check_is_ndarray(value, "vrnt_chrgrp_len")
            check_ndarray_dtype_is_integer(value, "vrnt_chrgrp_len")
            check_ndarray_ndim(value, "vrnt_chrgrp_len", 1)
        self._vrnt_chrgrp_len = value

    ################## Spline Properites ###################
    @property
    def spline(self) -> Union[dict,None]:
        """Interpolation spline(s)."""
        return self._spline
    @spline.setter
    def spline(self, value: Union[dict,None]) -> None:
        """Set interpolation spline(s)"""
        if value is not None:
            check_is_dict(value, "spline")
        self._spline = value

    ############# Spline Metadata Properites ###############
    @property
    def spline_kind(self) -> Union[str,None]:
        """Spline kind."""
        return self._spline_kind
    @spline_kind.setter
    def spline_kind(self, value: Union[str,None]) -> None:
        """Set the spline kind"""
        if value is not None:
            check_is_str(value, "spline_kind")
        self._spline_kind = value

    @property
    def spline_fill_value(self) -> object:
        """Default spline fill value."""
        return self._spline_fill_value
    @spline_fill_value.setter
    def spline_fill_value(self, value: object) -> None:
        """Set the default spline fill value"""
        self._spline_fill_value = value

    ############################## Object Methods ##############################

    ################## GeneticMap copying ##################
    def copy(
            self
        ) -> 'StandardGeneticMap':
        """
        Make a shallow copy of the StandardGeneticMap.

        Returns
        -------
        out : StandardGeneticMap
            A shallow copy of the original StandardGeneticMap.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: dict
        ) -> 'StandardGeneticMap':
        """
        Make a deep copy of the StandardGeneticMap.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : StandardGeneticMap
            A deep copy of the original StandardGeneticMap.
        """
        return copy.deepcopy(self, memo)

    ################### Sorting Methods ####################
    def lexsort(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a sequence of keys.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        # if no keys were provided, set a default
        if keys is None:
            # key priority: (1) vrnt_chrgrp (2) vrnt_phypos (3) vrnt_genpos
            keys = (self.vrnt_genpos, self.vrnt_phypos, self.vrnt_chrgrp)
        
        # if a tuple of keys was provided, check that all keys are ndarrays of 
        # acceptable length, then remove any None
        elif isinstance(keys, tuple):
            for i,key in enumerate(keys):
                check_is_ndarray(key, "keys[{0}]".format(i))
                check_ndarray_len_eq(key, "keys[{0}]".format(i), self.nvrnt)
        
        # if a ndarray was provided, check that its length is acceptable, 
        # then convert to tuple
        elif isinstance(keys, numpy.ndarray):
            check_ndarray_len_eq(keys, "keys", self.nvrnt)
            keys = (keys,)
        
        # otherwise raise type error
        else:
            raise TypeError(
                "argument '{0}' must be of type '{1}', '{2}', or '{3}' but received type '{4}'".format(
                    "keys",
                    tuple.__name__,
                    numpy.ndarray.__name__,
                    "None",
                    type(keys).__name__
                )
            )

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def reorder(
            self, 
            indices: Union[Sequence,numpy.ndarray]
        ) -> None:
        """
        Reorder markers in the GeneticMap using an array of indices.
        Note this modifies the GeneticMap in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # sort internal self
        self.vrnt_chrgrp = self.vrnt_chrgrp[indices]
        self.vrnt_phypos = self.vrnt_phypos[indices]
        self.vrnt_genpos = self.vrnt_genpos[indices]

    def sort(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None
        ) -> None:
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

    def group(
            self
        ) -> None:
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

    def is_grouped(
            self
        ) -> bool:
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
    def remove(
            self, 
            indices: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove indices from the GeneticMap. Sort and group internal arrays.

        Parameters
        ----------
        indices : int, slice, Sequence
            Array of shape ``(a,)``, ``slice`` or ``int`` of item(s) to remove.

            Where:

            - ``a`` is the number of indices to remove.
        kwargs : dict
            Additional keyword arguments.
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

    def select(
            self, 
            indices: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Keep only selected markers, removing all others from the GeneticMap.
        Sort and group internal arrays.

        Parameters
        ----------
        indices : int, slice, Sequence
            Array of shape ``(a,)``, ``slice`` or ``int`` of item(s) to remove.

            Where:

            - ``a`` is the number of indices to remove.
        kwargs : dict
            Additional keyword arguments.
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

    ################## Integrity Methods ###################
    def congruence(
            self
        ) -> numpy.ndarray:
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

    def is_congruent(
            self
        ) -> bool:
        """
        Determine if the genetic map is congruent
        Determine if all sites in the genetic map demonstrate congruence with
        their supposed physical and genetic positions.
        """
        # determine if all sites are congruent
        out = numpy.all(self.congruence())
        return out

    def remove_discrepancies(
            self
        ) -> None:
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
        if( numpy.all(mask) ):
            return

        # select (and group) only sites that are congruent
        self.select(mask)

    ################ Interpolation Methods #################
    # TODO: do not require sorted internals to build spline
    def build_spline(
            self, 
            kind: str = 'linear', 
            fill_value: Union[numpy.ndarray,str] = 'extrapolate', 
            **kwargs: dict
        ) -> None:
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
        kwargs : dict
            Additional keyword arguments.
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

    def has_spline(
            self
        ) -> bool:
        """Return whether or not the GeneticMap has a built spline."""
        return (
            (self._spline is not None) and
            (self._spline_kind is not None) and
            (self._spline_fill_value is not None)
        )

    def interp_genpos(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray
        ) -> numpy.ndarray:
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

        # raise warning if map is not congruent
        if not self.is_congruent():
            warnings.warn("genetic map is not congruent: markers are out of order", RuntimeWarning)

        # allocate empty memory
        out = numpy.empty(vrnt_phypos.shape, dtype='float64')

        # for each chromosome-position pair
        for i,(chrgrp,phypos) in enumerate(zip(vrnt_chrgrp, vrnt_phypos)):
            try:                                # try to index dict
                model = self._spline[chrgrp]    # try to get model
                out[i] = model(phypos)          # interpolate genetic position
            except KeyError:                    # if dict key not found
                out[i] = numpy.nan              # set to NaN

        return out

    def interp_gmap(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_stop: numpy.ndarray, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_fncode: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'StandardGeneticMap':
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
    def gdist1g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            ast: Optional[int] = None, 
            asp: Optional[int] = None
        ) -> numpy.ndarray:
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

    def gdist2g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            rst: Optional[int] = None, 
            rsp: Optional[int] = None, 
            cst: Optional[int] = None, 
            csp: Optional[int] = None
        ) -> numpy.ndarray:
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

    def gdist1p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            ast: Optional[int] = None, 
            asp: Optional[int] = None
        ) -> numpy.ndarray:
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

    def gdist2p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            rst: Optional[int] = None, 
            rsp: Optional[int] = None, 
            cst: Optional[int] = None, 
            csp: Optional[int] = None
        ) -> numpy.ndarray:
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

    ############################## Class Methods ###############################

    ############################# Static Methods ###############################



################################## Utilities ###################################
def check_is_StandardGeneticMap(v: object, vname: str) -> None:
    """
    Check if object is of type StandardGeneticMap. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, StandardGeneticMap):
        raise TypeError("'{0}' must be of type StandardGeneticMap.".format(vname))


class A:
    def __init__(self, x):
        self.x = x
    def __copy__(self):
        cls = self.__class__
        out = cls.__new__(cls)
        out.__init__(x = copy.copy(self.x))
        return out
        # return self.__class__(x = copy.copy(self.x))
    def copy(self):
        return copy.copy(self)

class B(A):
    pass

a = A(4)
b = B(5)
c = copy.copy(a)
d = copy.copy(b)
e = a.copy()
f = b.copy()
type(a)
type(b)
type(c)
type(d)
