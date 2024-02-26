"""
Module implementing a standard genetic map format and associated error checking
routines.
"""

__all__ = [
    "StandardGeneticMap",
    "check_is_StandardGeneticMap",
]

import copy
from numbers import Integral
from typing import Optional, Sequence, Tuple, Union
import warnings
import numpy
import pandas
from scipy.interpolate import interp1d
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_is_str_or_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_floating
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_integer
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_python import check_is_str_or_Integral
from pybrops.core.error.error_value_python import check_tuple_len_eq
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
            spline: Optional[dict] = None,
            spline_kind: str = "linear",
            spline_fill_value: Union[str,numpy.ndarray] = "extrapolate",
            vrnt_genpos_units: str = "M",
            auto_group: bool = True,
            auto_build_spline: bool = True,
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

        spline : dict, None, default = None
            Pre-built interpolation spline to associate with the genetic map.
        
        spline_kind : str, default = "linear"
            In automatic building of splines, the spline kind to be built.
            Specifies the kind of interpolation as a string ('linear',
            'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous',
            'next', where 'zero', 'slinear', 'quadratic' and 'cubic' refer to a
            spline interpolation of zeroth, first, second or third order;
            'previous' and 'next' simply return the previous or next value of
            the point) or as an integer specifying the order of the spline
            interpolator to use.
            
        spline_fill_value : str, numpy.ndarray, default = "extrapolate"
            In automatic building of splines, the spline fill value to use.
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

        vrnt_genpos_units : str, default = "M"
            Units in which genetic positions in the ``vrnt_genpos`` array are 
            stored. Options are listed below and are case-sensitive:
            
            - ``"M"`` - genetic position units are in Morgans
            - ``"Morgans"`` - genetic position units are in Morgans
            - ``"cM"`` - genetic position units are in centiMorgans
            - ``"centiMorgans"`` - genetic position units are in centiMorgans
            
            Internally, all genetic positions are stored in Morgans. Providing 
            the units of the input  

        auto_group : bool
            Whether to automatically sort and group variants into chromosome groups.
        
        auto_build_spline : bool
            Whether to automatically construct a spline on object construction.
            If ``spline`` is provided, then this spline is overwritten.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments!
        self.vrnt_chrgrp       = vrnt_chrgrp # must be assigned first!
        self.vrnt_phypos       = vrnt_phypos
        self.vrnt_genpos       = (vrnt_genpos,vrnt_genpos_units)
        self.spline            = spline
        self.spline_kind       = spline_kind
        self.spline_fill_value = spline_fill_value

        # set variant metadata to None
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len  = None

        # automatically group if requested
        if auto_group:
            self.group()
        
        # automatically build spline if requested
        if auto_build_spline:
            self.build_spline(self.spline_kind, self.spline_fill_value)

    def __len__(self) -> int:
        """
        Get the number of markers in the genetic map.
        
        Returns
        -------
        out : int
            The number of markers in the genetic map.
        """
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
            vrnt_chrgrp       = copy.copy(self.vrnt_chrgrp), 
            vrnt_phypos       = copy.copy(self.vrnt_phypos), 
            vrnt_genpos       = copy.copy(self.vrnt_genpos),
            spline            = copy.copy(self.spline),
            spline_kind       = copy.copy(self.spline_kind),
            spline_fill_value = copy.copy(self.spline_fill_value),
            vrnt_genpos_units = "M",
            auto_group        = False,
            auto_build_spline = False,
        )

        # copy variant metadata to the new object
        out.vrnt_chrgrp_name = copy.copy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.copy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.copy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len  = copy.copy(self.vrnt_chrgrp_len)

        return out

    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
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
            vrnt_chrgrp       = copy.deepcopy(self.vrnt_chrgrp,       memo), 
            vrnt_phypos       = copy.deepcopy(self.vrnt_phypos,       memo), 
            vrnt_genpos       = copy.deepcopy(self.vrnt_genpos,       memo),
            spline            = copy.deepcopy(self.spline,            memo),
            spline_kind       = copy.deepcopy(self.spline_kind,       memo),
            spline_fill_value = copy.deepcopy(self.spline_fill_value, memo),
            vrnt_genpos_units = "M",
            auto_group        = False,
            auto_build_spline = False,
        )

        # copy variant metadata to the new object
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name, memo)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix, memo)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix, memo)
        out.vrnt_chrgrp_len  = copy.deepcopy(self.vrnt_chrgrp_len,  memo)

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
    def vrnt_genpos(self, value: Union[numpy.ndarray,Tuple[numpy.ndarray,str]]) -> None:
        """Set variant genetic position array"""
        array = value
        units = "M"
        if isinstance(value, tuple):
            check_tuple_len_eq(value, "vrnt_genpos", 2)
            check_is_ndarray(value[0], "vrnt_genpos[0]")
            check_is_str(value[1], "vrnt_genpos[1]")
            array = value[0]
            units = value[1]
        if units in ("M","Morgans"):
            pass
        elif units in ("cM","centiMorgans"):
            array = 0.01 * array
        else:
            raise ValueError(
                "vrnt_genpos units must be 'M', 'Morgans', 'cM', or 'centiMorgans' but received value '{0}'".format(
                    units
                )
            )
        check_is_ndarray(array, "vrnt_genpos")
        check_ndarray_dtype_is_floating(array, "vrnt_genpos")
        check_ndarray_ndim(array, "vrnt_genpos", 1)
        check_ndarray_len_eq(array, "vrnt_genpos", len(self.vrnt_chrgrp))
        self._vrnt_genpos = array

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
    def spline_fill_value(self) -> Union[str,numpy.ndarray,None]:
        """Default spline fill value."""
        return self._spline_fill_value
    @spline_fill_value.setter
    def spline_fill_value(self, value: Union[str,numpy.ndarray,None]) -> None:
        """Set the default spline fill value"""
        if value is not None:
            check_is_str_or_ndarray(value, "spline_fill_value")
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
        return self.__copy__()

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'StandardGeneticMap':
        """
        Make a deep copy of the StandardGeneticMap.

        Parameters
        ----------
        memo : dict, None
            Dictionary of memo metadata.

        Returns
        -------
        out : StandardGeneticMap
            A deep copy of the original StandardGeneticMap.
        """
        return self.__deepcopy__(memo)

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
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : A (N,) ndarray of ints
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
        Reorder markers in-place in the GeneticMap using an array of indices.

        Parameters
        ----------
        indices : A (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.

        kwargs : dict
            Additional keyword arguments.
        """
        # sort internal self
        self.vrnt_chrgrp = self.vrnt_chrgrp[indices]
        self.vrnt_phypos = self.vrnt_phypos[indices]
        self.vrnt_genpos = self.vrnt_genpos[indices]

        # reset variant chromosome group name, stix, spix, len
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len = None

    def sort(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None
        ) -> None:
        """
        Sort slements of the GeneticMap using a sequence of keys.
        Note this modifies the GeneticMap in-place.

        Parameters
        ----------
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        # get indices for sort
        indices = self.lexsort(keys)

        # reorder internals
        self.reorder(indices)

    def group(
            self,
            **kwargs: dict
        ) -> None:
        """
        Sort the GeneticMap jointly by chromosome group and physical position, 
        then populate grouping indices.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # sort the genetic map using default keys
        self.sort()

        # get unique group names, start indices, group lengths
        uniq = numpy.unique(self._vrnt_chrgrp, return_index = True, return_counts = True)

        # make assignments
        self._vrnt_chrgrp_name, self._vrnt_chrgrp_stix, self._vrnt_chrgrp_len = uniq

        # calculate chr_grp_spix (stop indices)
        self._vrnt_chrgrp_spix = self._vrnt_chrgrp_stix + self._vrnt_chrgrp_len

    def ungroup(
            self,
            **kwargs: dict
        ) -> None:
        """
        Remove grouping metadata from the GeneticMap.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # set grouping metadat to None
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len  = None

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
        Remove indices from the GeneticMap. If the GeneticMap was grouped 
        beforehand, then re-sort and re-group internal arrays after removing 
        indices.

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
        self.vrnt_chrgrp = numpy.delete(self.vrnt_chrgrp, indices)
        self.vrnt_phypos = numpy.delete(self.vrnt_phypos, indices)
        self.vrnt_genpos = numpy.delete(self.vrnt_genpos, indices)

        # if GeneticMap was previously grouped, re-sort and re-group
        if self.is_grouped():
            self.group()

    def select(
            self, 
            indices: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Keep only selected markers, removing all others from the GeneticMap.
        If the GeneticMap was grouped beforehand, then re-sort and re-group 
        internal arrays after removing indices.

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
        self.vrnt_chrgrp = self.vrnt_chrgrp[indices]
        self.vrnt_phypos = self.vrnt_phypos[indices]
        self.vrnt_genpos = self.vrnt_genpos[indices]

        # if GeneticMap was previously grouped, re-sort and re-group
        if self.is_grouped():
            self.group()

    ################## Integrity Methods ###################
    def congruence(
            self
        ) -> numpy.ndarray:
        """
        Assess physical and genetic map site congruency. If the genetic map is 
        not grouped, it will be grouped. 

        Returns
        -------
        out : numpy.ndarray
            A boolean matrix of map concordancies where:

            - ``True`` = the current marker has a map_pos >= the previous 
              position
            - ``False`` = the current marker has a map_pos < the previous 
              position

        Notes
        -----
        This assumes high contiguity between physical and genetic maps
        (i.e. a high quality reference genome). This assumption may cause
        major issues if there are incorrect markers at the beginning of the
        chromosome. This also assumes the first marker on the chromosome is 
        placed correctly.
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
        Determine if all sites in the genetic map demonstrate congruence with
        their supposed physical and genetic positions.

        Returns
        -------
        out : bool
            Whether all genetic map loci demonstrate congruence between 
            physical and genetic positions.
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

        Notes
        -----
        This assumption may cause major issues if there are incorrect markers 
        at the beginning of the chromosome.
        """
        # get boolean mask for concordant map positions
        mask = self.congruence()

        # if all are concordant, exit this function
        if( numpy.all(mask) ):
            return

        # select (and group) only sites that are congruent
        self.select(mask)

    ################ Interpolation Methods #################
    def build_spline(
            self, 
            kind: str = 'linear', 
            fill_value: Union[str,numpy.ndarray] = 'extrapolate', 
            **kwargs: dict
        ) -> None:
        """
        Build a spline for estimating genetic map distances.

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
        self.spline = {}

        # set the spline type and fill value
        self.spline_kind = kind
        self.spline_fill_value = fill_value

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
        """
        Return whether or not the GeneticMap has a built spline.
        
        Returns
        -------
        out : bool
            Whether the GeneticMap has a spline built.
        """
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
        Interpolate genetic positions given variant physical positions.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            Chromosome/linkage group labels for each marker variant.
        
        vrnt_phypos : numpy.ndarray
            Chromosome/linkage group physical positions for each marker variant.

        Returns
        -------
        out : numpy.ndarray
            Interpolated genetic positions for each marker variant.
        """
        # raise error if no spline is found
        if not self.has_spline():
            raise ValueError("interpolation spline not built")

        # raise warning if map is not congruent
        if not self.is_congruent():
            warnings.warn("genetic map is not congruent: markers are out of order", RuntimeWarning)

        # allocate empty memory
        out = numpy.empty(vrnt_phypos.shape, dtype = float)

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
            **kwargs: dict
        ) -> 'StandardGeneticMap':
        """
        Interpolate a new genetic map from the current genetic map. Associate 
        spline of current GeneticMap with new GeneticMap.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            Chromosome/linkage group labels for each marker variant.
        
        vrnt_phypos : numpy.ndarray
            Chromosome/linkage group physical positions for each marker variant.

        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : StandardGeneticMap
            An interpolated genetic map sharing a copy of the spline from the 
            original genetic map.
        """
        # interpolate genetic positions
        vrnt_genpos = self.interp_genpos(vrnt_chrgrp, vrnt_phypos)

        # get the class of self
        cls = self.__class__

        # create a new object instance
        out = cls.__new__(cls)

        # initialize the object instance
        out.__init__(
            vrnt_chrgrp       = vrnt_chrgrp, 
            vrnt_phypos       = vrnt_phypos,
            vrnt_genpos       = vrnt_genpos,
            spline            = copy.deepcopy(self.spline),
            spline_kind       = copy.deepcopy(self.spline_kind),
            spline_fill_value = copy.deepcopy(self.spline_fill_value),
            vrnt_genpos_units = "M",
            auto_group        = False,
            auto_build_spline = False,
        )

        # copy variant metadata to the new object
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len  = copy.deepcopy(self.vrnt_chrgrp_len)

        return out

    ############### Genetic Distance Methods ###############
    def gdist1g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            ast: Optional[Integral] = None, 
            asp: Optional[Integral] = None
        ) -> numpy.ndarray:
        """
        Calculate sequential genetic distances using genetic positions.
        Requires ``vrnt_chrgrp`` and ``vrnt_genpos`` to have been sorted 
        jointly in ascending order.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
            Must be sorted in ascending order jointly with ``vrnt_genpos``.

        vrnt_genpos : numpy.ndarray
            A 1D array of variant genetic positions.
            Must be sorted in ascending order jointly with ``vrnt_chrgrp``.

        ast : Integral, None
            Optional array start index (inclusive).
            If ``None``, assume that all array elements are to be used for 
            sequential genetic distance calculations.

        asp : Integral, None
            Optional array stop index (exclusive).
            If ``None``, assume that all array elements are to be used for 
            sequential genetic distance calculations.

        Returns
        -------
        out : numpy.ndarray
            A 1D array of distances between the marker prior.

        Notes
        -----
        Sequential distance arrays will start every chromosome with numpy.inf!
        """
        # type checks
        check_is_ndarray(vrnt_chrgrp, "vrnt_chrgrp")
        check_ndarray_dtype_is_integer(vrnt_chrgrp, "vrnt_chrgrp")
        check_is_ndarray(vrnt_genpos, "vrnt_genpos")
        check_ndarray_dtype_is_floating(vrnt_genpos, "vrnt_genpos")

        # get views of only sections we'll be looking at
        view_chrgrp = vrnt_chrgrp[ast:asp]
        view_genpos = vrnt_genpos[ast:asp]

        # get chromosome groups, start index, length of region
        uniq, start, counts = numpy.unique(view_chrgrp, return_index = True, return_counts = True)

        # get stop index
        stop = start + counts

        # allocate empty array
        out = numpy.empty(view_genpos.shape, dtype = float)

        # for each chromosome group
        for st,sp in zip(start, stop):
            # infinite distance between chromosomes
            out[st] = numpy.inf
            # subtract distances
            out[st+1:sp] = view_genpos[st+1:sp] - view_genpos[st:sp-1]

        return out

    def gdist2g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            rst: Optional[Integral] = None, 
            rsp: Optional[Integral] = None, 
            cst: Optional[Integral] = None, 
            csp: Optional[Integral] = None
        ) -> numpy.ndarray:
        """
        Calculate pairwise genetic distances using genetic positions.
        Requires ``vrnt_chrgrp`` and ``vrnt_genpos`` to have been sorted 
        jointly in ascending order.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
            Must be sorted in ascending order jointly with ``vrnt_genpos``.

        vrnt_genpos : numpy.ndarray
            A 1D array of variant genetic positions.
            Must be sorted in ascending order jointly with ``vrnt_chrgrp``.

        rst : Integral, None
            Optional row start index (inclusive).
            If ``None``, assume that all rows are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        rsp : Integral, None
            Optional row stop index (exclusive).
            If ``None``, assume that all rows are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        cst : Integral, None
            Optional column start index (inclusive).
            If ``None``, assume that all columns are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        csp : Integral, None
            Optional column stop index (exclusive).
            If ``None``, assume that all columns are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        Returns
        -------
        out : numpy.ndarray
            A 2D array of distances between marker pairs.
        """
        # TODO: # OPTIMIZE: fill with inf, then add distances? - assumes sort
        # create sparse meshgrid indicating where chromosomes differ
        mi, mj = numpy.meshgrid(vrnt_chrgrp[rst:rsp], vrnt_chrgrp[cst:csp], indexing='ij', sparse=True)

        # create sparse meshgrid indicating where genetic positions are
        gi, gj = numpy.meshgrid(vrnt_genpos[rst:rsp], vrnt_genpos[cst:csp], indexing='ij', sparse=True)

        # calculate absolute distances
        out = numpy.abs(gi - gj)

        # where chromosomes are different, set to numpy.inf
        out[mi != mj] = numpy.inf

        return out

    def gdist1p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            ast: Optional[Integral] = None, 
            asp: Optional[Integral] = None
        ) -> numpy.ndarray:
        """
        Calculate sequential genetic distances using physical positions.
        Requires ``vrnt_chrgrp`` and ``vrnt_phypos`` to have been sorted 
        jointly in ascending order. Requires an interpolation spline to have 
        been built beforehand.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
            Must be sorted in ascending order jointly with ``vrnt_phypos``.

        vrnt_phypos : numpy.ndarray
            A 1D array of variant genetic positions.
            Must be sorted in ascending order jointly with ``vrnt_chrgrp``.

        ast : Integral, None
            Optional array start index (inclusive).
            If ``None``, assume that all array elements are to be used for 
            sequential genetic distance calculations.

        asp : Integral, None
            Optional array stop index (exclusive).
            If ``None``, assume that all array elements are to be used for 
            sequential genetic distance calculations.

        Returns
        -------
        out : numpy.ndarray
            A 1D array of distances between the marker prior.

        Notes
        -----
        Sequential distance arrays will start every chromosome with numpy.inf!
        """
        # type checks
        check_is_ndarray(vrnt_chrgrp, "vrnt_chrgrp")
        check_ndarray_dtype_is_integer(vrnt_chrgrp, "vrnt_chrgrp")
        check_is_ndarray(vrnt_phypos, "vrnt_phypos")
        check_ndarray_dtype_is_integer(vrnt_phypos, "vrnt_phypos")

        # TODO: # OPTIMIZE: don't interpolate all markers
        # interpolate genetic positions
        vrnt_genpos = self.interp_genpos(vrnt_chrgrp, vrnt_phypos)

        # calculate genetic distances
        out = self.gdist1g(vrnt_chrgrp, vrnt_genpos, ast, asp)

        return out

    def gdist2p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            rst: Optional[Integral] = None, 
            rsp: Optional[Integral] = None, 
            cst: Optional[Integral] = None, 
            csp: Optional[Integral] = None
        ) -> numpy.ndarray:
        """
        Calculate pairwise genetic distances using physical positions.
        Requires ``vrnt_chrgrp`` and ``vrnt_phypos`` to have been sorted 
        jointly in ascending order.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            A 1D array of variant chromosome groups.
            Must be sorted in ascending order jointly with ``vrnt_phypos``.

        vrnt_phypos : numpy.ndarray
            A 1D array of variant genetic positions.
            Must be sorted in ascending order jointly with ``vrnt_chrgrp``.

        rst : Integral, None
            Optional row start index (inclusive).
            If ``None``, assume that all rows are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        rsp : Integral, None
            Optional row stop index (exclusive).
            If ``None``, assume that all rows are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        cst : Integral, None
            Optional column start index (inclusive).
            If ``None``, assume that all columns are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        csp : Integral, None
            Optional column stop index (exclusive).
            If ``None``, assume that all columns are to be calculated in the 
            pairwise genetic distance matrix are to be calculated.

        Returns
        -------
        out : numpy.ndarray
            A 2D array of distances between marker pairs.
        """
        # TODO: # OPTIMIZE: don't interpolate all markers
        # interpolate genetic positions
        vrnt_genpos = self.interp_genpos(vrnt_chrgrp, vrnt_phypos)

        # calculate genetic distances
        out = self.gdist2g(vrnt_chrgrp, vrnt_genpos, rst, rsp, cst, csp)

        return out

    #################### Export Methods ####################
    def to_pandas(
            self, 
            vrnt_chrgrp_col: str = "chr", 
            vrnt_phypos_col: str = "pos", 
            vrnt_genpos_col: str = "cM", 
            vrnt_genpos_units: str = "cM",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a GeneticMap to a pandas.DataFrame.

        Parameters
        ----------
        vrnt_chrgrp_col : str, default = "chr"
            Name of the chromosome/linkage group name column to which to export.

        vrnt_phypos_col : str, default = "pos"
            Name of the physical position column to which to export.

        vrnt_genpos_col : str, default = "cM"
            Name of the genetic position column to which to export.

        vrnt_genpos_units : str, default = "cM"
            Units of the genetic position column to which to export.
            Options are listed below and are case-sensitive:
            
            - ``"M"`` - genetic position units are in Morgans
            - ``"Morgans"`` - genetic position units are in Morgans
            - ``"cM"`` - genetic position units are in centiMorgans
            - ``"centiMorgans"`` - genetic position units are in centiMorgans

        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            pandas.DataFrame.
        
        Returns
        -------
        out : pandas.DataFrame
            An output dataframe.
        """
        # calculate genetic positions
        vrnt_genpos_convert = self.vrnt_genpos
        if vrnt_genpos_units in ("M","Morgans"):
            pass
        elif vrnt_genpos_units in ("cM", "centiMorgans"):
            vrnt_genpos_convert = 100.0 * vrnt_genpos_convert
        else:
            raise ValueError(
                "vrnt_genpos units must be 'M', 'Morgans', 'cM', or 'centiMorgans' but received value '{0}'".format(
                    vrnt_genpos_units
                )
            )

        # construct output dataframe
        out = pandas.DataFrame({
            vrnt_chrgrp_col: self.vrnt_chrgrp,
            vrnt_phypos_col: self.vrnt_phypos,
            vrnt_genpos_col: vrnt_genpos_convert,
        })

        return out

    def to_csv(
            self, 
            filename: str,
            vrnt_chrgrp_col: str = "chr", 
            vrnt_phypos_col: str = "pos", 
            vrnt_genpos_col: str = "cM", 
            vrnt_genpos_units: str = "cM",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Write a StandardGeneticMap to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.

        vrnt_chrgrp_col : str, default = "chr"
            Name of the chromosome/linkage group name column to which to export.

        vrnt_phypos_col : str, default = "pos"
            Name of the physical position column to which to export.

        vrnt_genpos_col : str, default = "cM"
            Name of the genetic position column to which to export.

        vrnt_genpos_units : str, default = "cM"
            Units of the genetic position column to which to export.
            Options are listed below and are case-sensitive:
            
            - ``"M"`` - genetic position units are in Morgans
            - ``"Morgans"`` - genetic position units are in Morgans
            - ``"cM"`` - genetic position units are in centiMorgans
            - ``"centiMorgans"`` - genetic position units are in centiMorgans

        sep : str, default = ","
            Separator to use in the exported CSV file.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV file.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # convert GeneticMap to pandas.DataFrame
        df = self.to_pandas(
            vrnt_chrgrp_col = vrnt_chrgrp_col,
            vrnt_phypos_col = vrnt_phypos_col,
            vrnt_genpos_col = vrnt_genpos_col,
            vrnt_genpos_units = vrnt_genpos_units,
        )

        # use pandas to save to csv
        df.to_csv(
            path_or_buf = filename,
            sep = sep,
            header = header,
            index = index,
            **kwargs
        )

    ############################## Class Methods ###############################

    #################### Import Methods ####################
    @classmethod
    def from_pandas(
            cls,
            df: pandas.DataFrame,
            vrnt_chrgrp_col: Union[str,Integral] = "chr", 
            vrnt_phypos_col: Union[str,Integral] = "pos", 
            vrnt_genpos_col: Union[str,Integral] = "cM", 
            spline: dict = None,
            spline_kind: str = "linear",
            spline_fill_value: Union[str,numpy.ndarray] = "extrapolate",
            vrnt_genpos_units: str = "M",
            auto_group: bool = True,
            auto_build_spline: bool = True,
            **kwargs: dict
        ) -> 'StandardGeneticMap':
        """
        Read a StandardGeneticMap from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        vrnt_chrgrp_col : str, Integral, default = "chr"
            Name or number of the chromosome/linkage group name column from which to import.

        vrnt_phypos_col : str, Integral, default = "pos"
            Name or number of the physical position column from which to import.

        vrnt_genpos_col : str, Integral, default = "cM"
            Name or number of the genetic position column from which to import.

        spline : dict, None, default = None
            Pre-built interpolation spline to associate with the genetic map.
        
        spline_kind : str, default = "linear"
            In automatic building of splines, the spline kind to be built.
            Specifies the kind of interpolation as a string ('linear',
            'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous',
            'next', where 'zero', 'slinear', 'quadratic' and 'cubic' refer to a
            spline interpolation of zeroth, first, second or third order;
            'previous' and 'next' simply return the previous or next value of
            the point) or as an integer specifying the order of the spline
            interpolator to use.
            
        spline_fill_value : str, numpy.ndarray, default = "extrapolate"
            In automatic building of splines, the spline fill value to use.
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

        vrnt_genpos_units : str, default = "M"
            Units in which genetic positions in the ``vrnt_genpos`` array are 
            stored. Options are listed below and are case-sensitive:
            
            - ``"M"`` - genetic position units are in Morgans
            - ``"Morgans"`` - genetic position units are in Morgans
            - ``"cM"`` - genetic position units are in centiMorgans
            - ``"centiMorgans"`` - genetic position units are in centiMorgans
            
            Internally, all genetic positions are stored in Morgans. Providing 
            the units of the input  

        auto_group : bool
            Whether to automatically sort and group variants into chromosome groups.
        
        auto_build_spline : bool
            Whether to automatically construct a spline on object construction.
            If ``spline`` is provided, then this spline is overwritten.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            pandas.DataFrame.

        Returns
        -------
        out : StandardGeneticMap
            A StandardGeneticMap read from a pandas.DataFrame.
        """
        # data type checks
        check_is_pandas_DataFrame(df, "df")
        check_is_str_or_Integral(vrnt_chrgrp_col, "vrnt_chrgrp_col")
        check_is_str_or_Integral(vrnt_phypos_col, "vrnt_phypos_col")
        check_is_str_or_Integral(vrnt_genpos_col, "vrnt_genpos_col")

        # extract data (pandas.Series data type)
        vrnt_chrgrp = df[vrnt_chrgrp_col] if isinstance(vrnt_chrgrp_col, str) else df.iloc[:,vrnt_chrgrp_col]
        vrnt_phypos = df[vrnt_phypos_col] if isinstance(vrnt_phypos_col, str) else df.iloc[:,vrnt_phypos_col]
        vrnt_genpos = df[vrnt_genpos_col] if isinstance(vrnt_genpos_col, str) else df.iloc[:,vrnt_genpos_col]

        # convert to numpy
        vrnt_chrgrp = vrnt_chrgrp.to_numpy(dtype = int)
        vrnt_phypos = vrnt_phypos.to_numpy(dtype = int)
        vrnt_genpos = vrnt_genpos.to_numpy(dtype = float)
        
        # construct object
        out = cls(
            vrnt_chrgrp       = vrnt_chrgrp,
            vrnt_phypos       = vrnt_phypos,
            vrnt_genpos       = vrnt_genpos,
            spline            = spline,
            spline_kind       = spline_kind,
            spline_fill_value = spline_fill_value,
            vrnt_genpos_units = vrnt_genpos_units,
            auto_group        = auto_group,
            auto_build_spline = auto_build_spline,
            **kwargs
        )

        return out

    @classmethod
    def from_csv(
            cls, 
            filename: str,
            vrnt_chrgrp_col: Union[str,Integral] = "chr", 
            vrnt_phypos_col: Union[str,Integral] = "pos", 
            vrnt_genpos_col: Union[str,Integral] = "cM", 
            spline: dict = None,
            spline_kind: str = "linear",
            spline_fill_value: Union[str,numpy.ndarray] = "extrapolate",
            vrnt_genpos_units: str = "M",
            auto_group: bool = True,
            auto_build_spline: bool = True,
            sep: str = ',', 
            header: int = 0, 
            **kwargs: dict
        ) -> 'StandardGeneticMap':
        """
        Read a StandardGeneticMap from a CSV file.

        Parameters
        ----------
        filename : str, path object, or file-like object
            Any valid string path, including URLs. Valid URL schemes include 
            http, ftp, s3, and file. For file URLs, a host is expected (see 
            pandas docs).
        
        vrnt_chrgrp_col : str, Integral, default = "chr"
            Name or number of the chromosome/linkage group name column from 
            which to import.

        vrnt_phypos_col : str, Integral, default = "pos"
            Name or number of the physical position column from which to import.

        vrnt_genpos_col : str, Integral, default = "cM"
            Name or number of the genetic position column from which to import.

        spline : dict, None, default = None
            Pre-built interpolation spline to associate with the genetic map.
        
        spline_kind : str, default = "linear"
            In automatic building of splines, the spline kind to be built.
            Specifies the kind of interpolation as a string ('linear',
            'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous',
            'next', where 'zero', 'slinear', 'quadratic' and 'cubic' refer to a
            spline interpolation of zeroth, first, second or third order;
            'previous' and 'next' simply return the previous or next value of
            the point) or as an integer specifying the order of the spline
            interpolator to use.
            
        spline_fill_value : str, numpy.ndarray, default = "extrapolate"
            In automatic building of splines, the spline fill value to use.
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

        vrnt_genpos_units : str, default = "M"
            Units in which genetic positions in the ``vrnt_genpos`` array are 
            stored. Options are listed below and are case-sensitive:
            
            - ``"M"`` - genetic position units are in Morgans
            - ``"Morgans"`` - genetic position units are in Morgans
            - ``"cM"`` - genetic position units are in centiMorgans
            - ``"centiMorgans"`` - genetic position units are in centiMorgans
            
            Internally, all genetic positions are stored in Morgans. Providing 
            the units of the input  

        auto_group : bool
            Whether to automatically sort and group variants into chromosome groups.
        
        auto_build_spline : bool
            Whether to automatically construct a spline on object construction.
            If ``spline`` is provided, then this spline is overwritten.

        sep : str, default = ','
            CSV delimiter to use.
        
        header : int, list of int, default=0
            Row number(s) to use as the column names, and the start of the data.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : StandardGeneticMap
            A StandardGeneticMap read from a CSV file.
        """
        # read file using pandas
        df = pandas.read_csv(
            filename,
            sep = sep,
            header = header,
            **kwargs
        )

        # construct genetic map from pandas.DataFrame
        out = cls.from_pandas(
            df = df,
            vrnt_chrgrp_col = vrnt_chrgrp_col,
            vrnt_phypos_col = vrnt_phypos_col,
            vrnt_genpos_col = vrnt_genpos_col,
            spline = spline,
            spline_kind = spline_kind,
            spline_fill_value = spline_fill_value,
            vrnt_genpos_units = vrnt_genpos_units,
            auto_group = auto_group,
            auto_build_spline = auto_build_spline,
        )

        return out



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
        raise TypeError(
            "variable '{0}' must be of type '{1}' but received type '{2}'.".format(
                vname,
                StandardGeneticMap.__name__,
                type(v).__name__
            )
        )
