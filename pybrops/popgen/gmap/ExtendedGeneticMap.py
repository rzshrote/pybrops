"""
Module implementing a custom, extended genetic map format and associated error
checking routines.
"""

__all__ = [
    "ExtendedGeneticMap",
    "check_is_ExtendedGeneticMap",
]

import copy
from numbers import Integral
from typing import Optional, Sequence, Union
import numpy
import math
import pandas
import warnings
from scipy.interpolate import interp1d

from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_floating, check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_python import check_is_not_None
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_numpy import check_ndarray_size
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_object
from pybrops.popgen.gmap.GeneticMap import GeneticMap

class ExtendedGeneticMap(GeneticMap):
    """
    A concrete class for representing an extended genetic map format.

    The purpose of this concrete class is to implement functionality for:
        1) Extended genetic map representation.
        2) Extended genetic map metadata.
        3) Extended genetic map routines.
        4) Extended genetic map interpolation spline construction.
        5) Extended genetic map spline interpolation.
        6) Import and export of extended genetic maps.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_stop: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_fncode: Optional[numpy.ndarray] = None, 
            auto_group: bool = True,
            auto_build_spline: bool = True,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class ExtendedGeneticMap.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
        vrnt_phypos : numpy.ndarray
        vrnt_stop : numpy.ndarray
        vrnt_genpos : numpy.ndarray
        vrnt_name : numpy.ndarray, None
        vrnt_fncode : numpy.ndarray, None
        kwargs : dict
            Additional keyword arguments.
        """
        super(ExtendedGeneticMap, self).__init__(**kwargs)
        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_stop = vrnt_stop
        self.vrnt_genpos = vrnt_genpos
        self.vrnt_name = vrnt_name
        self.vrnt_fncode = vrnt_fncode
        # TODO: check all lengths equivalent

        # set variant metadata to None
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len  = None

        # set spline metadata to None
        self.spline = None
        self.spline_kind = None
        self.spline_fill_value = None

        # optionally automatically group variants
        if auto_group:
            self.group()
        
        # optionally automatically build interpolation splines
        if auto_build_spline:
            self.build_spline(**kwargs)

    def __len__(self):
        """Get the number of markers in the genetic map."""
        return len(self._vrnt_genpos)

    ################## GeneticMap copying ##################
    def __copy__(
            self
        ) -> 'ExtendedGeneticMap':
        """
        Make a shallow copy of the ExtendedGeneticMap.

        Returns
        -------
        out : ExtendedGeneticMap
            A shallow copy of the original ExtendedGeneticMap.
        """
        # construct a copy of the genetic map
        out = self.__class__(
            vrnt_chrgrp       = copy.copy(self.vrnt_chrgrp), 
            vrnt_phypos       = copy.copy(self.vrnt_phypos), 
            vrnt_stop         = copy.copy(self.vrnt_stop), 
            vrnt_genpos       = copy.copy(self.vrnt_genpos), 
            vrnt_name         = copy.copy(self.vrnt_name), 
            vrnt_fncode       = copy.copy(self.vrnt_fncode), 
            auto_group        = False,
            auto_build_spline = False,
        )

        # copy the variant metadata to the new object
        out.vrnt_chrgrp_name = copy.copy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.copy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.copy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len  = copy.copy(self.vrnt_chrgrp_len)

        # copy the spline metadata to the new object
        out.spline            = copy.copy(self.spline)
        out.spline_kind       = copy.copy(self.spline_kind)
        out.spline_fill_value = copy.copy(self.spline_fill_value)

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'ExtendedGeneticMap':
        """
        Make a deep copy of the ExtendedGeneticMap.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : ExtendedGeneticMap
            A deep copy of the original ExtendedGeneticMap.
        """
        # construct a deep copy of the genetic map
        out = self.__class__(
            vrnt_chrgrp       = copy.deepcopy(self.vrnt_chrgrp, memo), 
            vrnt_phypos       = copy.deepcopy(self.vrnt_phypos, memo), 
            vrnt_stop         = copy.deepcopy(self.vrnt_stop  , memo), 
            vrnt_genpos       = copy.deepcopy(self.vrnt_genpos, memo), 
            vrnt_name         = copy.deepcopy(self.vrnt_name  , memo), 
            vrnt_fncode       = copy.deepcopy(self.vrnt_fncode, memo), 
            auto_group        = False,
            auto_build_spline = False,
        )

        # deep copy the variant metadata to the new object
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name, memo)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix, memo)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix, memo)
        out.vrnt_chrgrp_len  = copy.deepcopy(self.vrnt_chrgrp_len , memo)

        # deep copy the spline metadata to the new object
        out.spline            = copy.copy(self.spline           , memo)
        out.spline_kind       = copy.copy(self.spline_kind      , memo)
        out.spline_fill_value = copy.copy(self.spline_fill_value, memo)

        return out

    ############################ Object Properties #############################

    ################### Data Properites ####################
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
        self._vrnt_phypos = value

    @property
    def vrnt_genpos(self) -> numpy.ndarray:
        """Variant genetic position."""
        return self._vrnt_genpos
    @vrnt_genpos.setter
    def vrnt_genpos(self, value: numpy.ndarray) -> None:
        """Set variant genetic position array"""
        check_is_ndarray(value, "vrnt_genpos")
        check_ndarray_dtype_is_floating(value, "vrnt_genpos")
        check_ndarray_ndim(value, "vrnt_genpos", 1)
        self._vrnt_genpos = value

    @property
    def vrnt_stop(self) -> numpy.ndarray:
        """Variant physical position stop position."""
        return self._vrnt_stop
    @vrnt_stop.setter
    def vrnt_stop(self, value: numpy.ndarray) -> None:
        """Set variant physical position stop position array"""
        check_is_ndarray(value, "vrnt_stop")
        check_ndarray_dtype_is_integer(value, "vrnt_stop")
        check_ndarray_ndim(value, "vrnt_stop", 1)
        self._vrnt_stop = value

    @property
    def vrnt_name(self) -> Union[numpy.ndarray,None]:
        """Variant names."""
        return self._vrnt_name
    @vrnt_name.setter
    def vrnt_name(self, value: Union[numpy.ndarray,None]) -> None:
        """Set variant names."""
        if value is not None:
            check_is_ndarray(value, "vrnt_name")
            check_ndarray_dtype_is_object(value, "vrnt_name")
            check_ndarray_ndim(value, "vrnt_name", 1)
        self._vrnt_name = value

    @property
    def vrnt_fncode(self) -> Union[numpy.ndarray,None]:
        """Variant function codes."""
        return self._vrnt_fncode
    @vrnt_fncode.setter
    def vrnt_fncode(self, value: Union[numpy.ndarray,None]) -> None:
        """Set variant function codes."""
        if value is not None:
            check_is_ndarray(value, "vrnt_fncode")
            check_ndarray_dtype_is_object(value, "vrnt_fncode")
            check_ndarray_ndim(value, "vrnt_fncode", 1)
        self._vrnt_fncode = value

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
        ) -> 'ExtendedGeneticMap':
        """
        Make a shallow copy of the ExtendedGeneticMap.

        Returns
        -------
        out : ExtendedGeneticMap
            A shallow copy of the original ExtendedGeneticMap.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: dict
        ) -> 'ExtendedGeneticMap':
        """
        Make a deep copy of the ExtendedGeneticMap.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : ExtendedGeneticMap
            A deep copy of the original ExtendedGeneticMap.
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
            keys = (
                self._vrnt_genpos,  # 3rd priority
                self._vrnt_phypos,  # 2nd priority
                self._vrnt_chrgrp   # 1st priority
            )
        else:
            # check that matrix lengths are the same
            for i,k in enumerate(keys):
                check_ndarray_size(k, "key"+i, len(self))

        # build tuple
        keys = tuple(k for k in keys if k is not None)

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
        self._vrnt_chrgrp = self._vrnt_chrgrp[indices]
        self._vrnt_phypos = self._vrnt_phypos[indices]
        self._vrnt_stop = self._vrnt_stop[indices]
        self._vrnt_genpos = self._vrnt_genpos[indices]
        if self._vrnt_name is not None:
            self._vrnt_name = self._vrnt_name[indices]
        if self._vrnt_fncode is not None:
            self._vrnt_fncode = self._vrnt_fncode[indices]

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

    def prune(
            self, 
            nt: int = None, 
            M: float = None
        ) -> None:
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
        kwargs : dict
            Additional keyword arguments.
        """
        # check if we have acceptible inputs
        if (nt is None) and (M is None):
            raise ValueError("'nt' and 'M' cannot both be None")

        # if not grouped, sort and group
        if not self.is_grouped():
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
        ) -> 'ExtendedGeneticMap':
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

    #################### Export Methods ####################
    def to_pandas_df(
            self
        ) -> pandas.DataFrame:
        """
        Convert a GeneticMap object to a pandas DataFrame.

        Returns
        -------
        df : pandas.DataFrame
            A pandas DataFrame containing genetic map data.
        """
        # create a dictionary of our values
        df_dict = {
            "chr_grp" : self.vrnt_chrgrp,
            "chr_start" : self.vrnt_phypos,
            "chr_stop" : self.vrnt_stop,
            "map_pos" : self.vrnt_genpos,
            "mkr_name" : self.vrnt_name,
            "map_fncode" : self.vrnt_fncode
        }

        # make the DataFrame
        df = pandas.DataFrame(df_dict)

        # return it
        return df

    def to_csv(
            self, 
            fname: str, 
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Convert a GeneticMap object to a csv file.
        """
        # convert GeneticMap to DataFrame
        df = self.to_pandas_df()

        # write to csv
        df.to_csv(
            fname,
            sep = sep,
            header = header,
            index = index
        )

    def to_egmap(
            self, 
            fname: str
        ) -> None:
        # write as special instance of csv
        self.to_csv(
            fname = fname,
            sep = '\t',
            header = True,
            index = False
        )

    ############################## Class Methods ###############################
    @classmethod
    def from_pandas_df(
            cls, 
            pandas_df: pandas.DataFrame, 
            vrnt_chrgrp_ix: int = 0, 
            vrnt_phypos_ix: int = 1, 
            vrnt_stop_ix: int = 2, 
            vrnt_genpos_ix: int = 3, 
            vrnt_name_ix: Optional[int] = None, 
            vrnt_fncode_ix: Optional[int] = None,
            auto_group: bool = True,
            auto_build_spline: bool = True
        ) -> 'ExtendedGeneticMap':
        """
        Read genetic map data from a Pandas DataFrame.

        Parameters
        ----------
        pandas_df : pandas.DataFrame
            A Pandas DataFrame.
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
        auto_group : bool
            Whether to automatically group and sort markers into linkage groups after 
            loading the ``ExtendedGeneticMap``.
        auto_build_spline : bool
            Whether to automatically construct an interpolation spline after 
            loading the ``ExtendedGeneticMap``. Interpolation spline type is the default
            value for the ``build_spline`` method.

        Returns
        -------
        genetic_map : ExtendedGeneticMap
            An ExtendedGeneticMap object containing all data required for
            genetic inferences.

        Notes
        -----
        Genetic map assumptions:

        - This assumes that we have a high quality genome assembly
          with near complete chromosome pseudomolecules and low number of
          errors.
        - This also assumes that we have a high quality genetic map
          with minimal inversions and mis-assemblies.
        - The number of linkage groups in the genetic map should equal the
          number of whole chromosomes in our genome assembly.
        - For discrepancies in the physical map vs. the genetic map, we assume
          that the physical map is correct.

        Pandas DataFrame field specifications:

        1) vrnt_chrgrp (REQUIRED)
            Name of the chromosome; equivalent to the linkage group.
            This is of type 'str'.
        2) vrnt_phypos (REQUIRED)
            The start position of the feature on the chromosome or scaffold.
            This is 1-indexed (e.g. the first base of a chromosome is 1) and
            inclusive (e.g. chr_start <= sequence <= chr_stop).
            This is an integer type.
        3) vrnt_stop (REQUIRED)
            The stop position of the feature on the chromosome or scaffold. This
            is 1-indexed (e.g. the first base of a chromosome is 1) and
            inclusive (e.g. chr_start <= sequence <= chr_stop).
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

            - Haldane: 'H'
            - Kosambi: 'K'
            - Unknown: 'U'
            - Custom: <str of any length>
        """
        # check to make sure several indices aren't None.
        check_is_not_None(vrnt_chrgrp_ix, "vrnt_chrgrp_ix")
        check_is_not_None(vrnt_phypos_ix, "vrnt_phypos_ix")
        check_is_not_None(vrnt_stop_ix, "vrnt_stop_ix")
        check_is_not_None(vrnt_genpos_ix, "vrnt_genpos_ix")

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
        vrnt_name = numpy.object_([str(e) for e in vrnt_name]) if vrnt_name is not None else None
        vrnt_fncode = numpy.object_([str(e) for e in vrnt_fncode]) if vrnt_fncode is not None else None

        # construct the gmap object
        genetic_map = cls(
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_stop = vrnt_stop,
            vrnt_genpos = vrnt_genpos,
            vrnt_name = vrnt_name,
            vrnt_fncode = vrnt_fncode,
            auto_group = auto_group,
            auto_build_spline = auto_build_spline
        )

        return genetic_map

    @classmethod
    def from_csv(
            cls, 
            fpath: str, 
            sep: str = ',', 
            header: int = 0, 
            vrnt_chrgrp_ix: int = 0, 
            vrnt_phypos_ix: int = 1, 
            vrnt_stop_ix: int = 2, 
            vrnt_genpos_ix: int = 3, 
            vrnt_name_ix: Optional[int] = None, 
            vrnt_fncode_ix: Optional[int] = None,
            auto_group: bool = True,
            auto_build_spline: bool = True
        ) -> 'ExtendedGeneticMap':
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
        auto_group : bool
            Whether to automatically group and sort markers into linkage groups on 
            loading the ``ExtendedGeneticMap``.
        auto_build_spline : bool
            Whether to automatically construct an interpolation spline after 
            loading the ``ExtendedGeneticMap``. Interpolation spline type is the default
            value for the ``build_spline`` method.

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

        genetic_map = cls.from_pandas_df(
            pandas_df = df,
            vrnt_chrgrp_ix = vrnt_chrgrp_ix,
            vrnt_phypos_ix = vrnt_phypos_ix,
            vrnt_stop_ix = vrnt_stop_ix,
            vrnt_genpos_ix = vrnt_genpos_ix,
            vrnt_name_ix = vrnt_name_ix,
            vrnt_fncode_ix = vrnt_fncode_ix,
            auto_group = auto_group,
            auto_build_spline = auto_build_spline
        )

        return genetic_map

    @classmethod
    def from_egmap(
            cls, 
            fpath: str,
            auto_group: bool = True,
            auto_build_spline: bool = True
        ) -> 'ExtendedGeneticMap':
        """
        Read an extended genetic map file (.egmap).

        Parameters
        ----------
        fpath : str, path object, or file-like object
            Any valid string path, including URLs. Valid URL schemes include http,
            ftp, s3, and file. For file URLs, a host is expected. (see pandas docs)
        auto_group : bool
            Whether to automatically group and sort markers into linkage groups on 
            loading the ``ExtendedGeneticMap``.
        auto_build_spline : bool
            Whether to automatically construct an interpolation spline after 
            loading the ``ExtendedGeneticMap``. Interpolation spline type is the default
            value for the ``build_spline`` method.

        Returns
        -------
        genetic_map : ExtendedGeneticMap
            A gmap object containing all data required for genetic inferences.

        Notes
        -----
        Extended genetic map file (.egmap) format (similar to BED file format).

        Genetic map assumptions:

        - This file format assumes that we have a high quality genome assembly
          with near complete chromosome pseudomolecules and low number of
          errors.
        - This file format also assumes that we have a high quality genetic map
          with minimal inversions and mis-assemblies.
        - The number of linkage groups in the genetic map should equal the
          number of whole chromosomes in our genome assembly.
        - For discrepancies in the physical map vs. the genetic map, we assume
          that the physical map is correct.

        General Extended Genetic Map File (.egmap) specifications:

        - This file format is headerless.
        - This file format is tab delimited.

        Extended Genetic Map File (.egmap) field specifications:

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

            - Haldane: 'H'
            - Kosambi: 'K'
            - Unknown: 'U'
            - Custom: <str of any length>
        """
        # read the file using pandas
        df = pandas.read_csv(
            fpath,          # file path
            sep = '\t',     # tab delimited
            header = 0      # there is a header
        )

        genetic_map = cls.from_pandas_df(
            pandas_df = df,
            vrnt_chrgrp_ix = 0,
            vrnt_phypos_ix = 1,
            vrnt_stop_ix = 2,
            vrnt_genpos_ix = 3,
            vrnt_name_ix = 4 if "mkr_name" in df.columns else None,
            vrnt_fncode_ix = 5 if "map_fncode" in df.columns else None,
            auto_group = auto_group,
            auto_build_spline = auto_build_spline
        )

        return genetic_map



################################## Utilities ###################################
def check_is_ExtendedGeneticMap(v: object, vname: str) -> None:
    if not isinstance(v, ExtendedGeneticMap):
        raise TypeError("'%s' must be an ExtendedGeneticMap." % vname)
