"""
Module defining basal genetic map interfaces and associated error checking routines.
"""

__all__ = [
    "GeneticMap",
    "check_is_GeneticMap",
]

from abc import ABCMeta, abstractmethod
from numbers import Integral
from typing import Optional, Union
import numpy
import pandas
from pybrops.core.io.CSVInputOutput import CSVInputOutput

from pybrops.core.io.PandasInputOutput import PandasInputOutput

class GeneticMap(PandasInputOutput,CSVInputOutput,metaclass=ABCMeta):
    """
    An abstract class for genetic map objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Genetic map representation.
        2) Genetic map metadata.
        3) Genetic map routines.
        4) Genetic map interpolation spline construction.
        5) Genetic map spline interpolation.
        6) Import and export of genetic maps.
    """

    ########################## Special Object Methods ##########################
    @abstractmethod
    def __len__(self):
        """Get the number of markers in the genetic map."""
        raise NotImplementedError("method is abstract")

    def __repr__(
            self
        ) -> str:
        """
        Return repr(self).
        
        Returns
        -------
        out : str
            A representation of the object.
        """
        return "<{0} of length {1} at {2}>".format(
            type(self).__name__,
            len(self),
            hex(id(self))
        )

    ################## GeneticMap copying ##################
    @abstractmethod
    def __copy__(
            self
        ) -> 'GeneticMap':
        """
        Make a shallow copy of the GeneticMap.

        Returns
        -------
        out : GeneticMap
            A shallow copy of the original GeneticMap.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'GeneticMap':
        """
        Make a deep copy of the GeneticMap.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : GeneticMap
            A deep copy of the original GeneticMap.
        """
        raise NotImplementedError("method is abstract")

    ############################ Object Properties #############################

    ################### Data Properites ####################
    @property
    @abstractmethod
    def nvrnt(self) -> Integral:
        """Number of variants in the GeneticMap."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_chrgrp(self) -> object:
        """Variant chromosome group label."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp.setter
    @abstractmethod
    def vrnt_chrgrp(self, value: object) -> None:
        """Set variant chromosome group label array"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def vrnt_phypos(self) -> object:
        """Variant physical position."""
        raise NotImplementedError("property is abstract")
    @vrnt_phypos.setter
    @abstractmethod
    def vrnt_phypos(self, value: object) -> None:
        """Set variant physical position array"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def vrnt_genpos(self) -> object:
        """Variant genetic position."""
        raise NotImplementedError("property is abstract")
    @vrnt_genpos.setter
    @abstractmethod
    def vrnt_genpos(self, value: object) -> None:
        """Set variant genetic position array"""
        raise NotImplementedError("property is abstract")
    
    ################# Metadata Properites ##################
    @property
    @abstractmethod
    def vrnt_chrgrp_name(self) -> object:
        """Variant chromosome group names."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_name.setter
    @abstractmethod
    def vrnt_chrgrp_name(self, value: object) -> None:
        """Set variant chromosome group name array"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def vrnt_chrgrp_stix(self) -> object:
        """Variant chromosome group start indices."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_stix.setter
    @abstractmethod
    def vrnt_chrgrp_stix(self, value: object) -> None:
        """Set variant chromosome group start indices array"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def vrnt_chrgrp_spix(self) -> object:
        """Variant chromosome group stop indices."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_spix.setter
    @abstractmethod
    def vrnt_chrgrp_spix(self, value: object) -> None:
        """Set variant chromosome group stop indices array"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def vrnt_chrgrp_len(self) -> object:
        """Variant chromosome group length."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_len.setter
    @abstractmethod
    def vrnt_chrgrp_len(self, value: object) -> None:
        """Set variant chromosome group length array"""
        raise NotImplementedError("property is abstract")

    ################## Spline Properites ###################
    @property
    @abstractmethod
    def spline(self) -> object:
        """Interpolation spline(s)."""
        raise NotImplementedError("property is abstract")
    @spline.setter
    @abstractmethod
    def spline(self, value: object) -> None:
        """Set interpolation spline(s)"""
        raise NotImplementedError("property is abstract")

    ############# Spline Metadata Properites ###############
    @property
    @abstractmethod
    def spline_kind(self) -> object:
        """Spline kind."""
        raise NotImplementedError("property is abstract")
    @spline_kind.setter
    @abstractmethod
    def spline_kind(self, value: object) -> None:
        """Set the spline kind"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def spline_fill_value(self) -> object:
        """Default spline fill value."""
        raise NotImplementedError("property is abstract")
    @spline_fill_value.setter
    @abstractmethod
    def spline_fill_value(self, value: object) -> None:
        """Set the default spline fill value"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ################## GeneticMap copying ##################
    @abstractmethod
    def copy(
            self
        ) -> 'GeneticMap':
        """
        Make a shallow copy of the GeneticMap.

        Returns
        -------
        out : GeneticMap
            A shallow copy of the original GeneticMap.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def deepcopy(
            self, 
            memo: dict
        ) -> 'GeneticMap':
        """
        Make a deep copy of the GeneticMap.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : GeneticMap
            A deep copy of the original GeneticMap.
        """
        raise NotImplementedError("method is abstract")

    ################### Sorting Methods ####################
    @abstractmethod
    def lexsort(
            self, 
            keys, 
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def reorder(
            self, 
            indices, 
            **kwargs: dict
        ) -> None:
        """
        Reorder markers in-place in the GeneticMap using an array of indices.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.

        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def sort(
            self, 
            keys, 
            **kwargs: dict
        ) -> None:
        """
        Sort slements of the GeneticMap using a sequence of keys.
        Note this modifies the GeneticMap in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Grouping Methods ###################
    @abstractmethod
    def group(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Sort the GeneticMap, then populate grouping indices.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def is_grouped(
            self, 
            **kwargs: dict
        ) -> bool:
        """
        Determine whether the GeneticMap has been sorted and grouped.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        raise NotImplementedError("method is abstract")

    ################ Insert/Delete Methods #################
    @abstractmethod
    def remove(
            self, 
            indices, 
            **kwargs: dict
        ) -> None:
        """
        Remove indices from the GeneticMap. Sort and group internal arrays.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape ``(a,)``, ``slice`` or ``int`` of item(s) to remove.

            Where:

            - ``a`` is the number of indices to remove.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def select(
            self, 
            indices, 
            **kwargs: dict
        ) -> None:
        """
        Keep only selected markers, removing all others from the GeneticMap.
        Sort and group internal arrays.

        Parameters
        ----------
        indices : numpy.ndarray, slice, int
            Array of shape ``(a,)``, ``slice`` or ``int`` of item(s) to remove.

            Where:

            - ``a`` is the number of indices to remove.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    # TODO: def insert(self, ...)

    ################## Integrity Methods ###################
    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    ################# Interpolation Methods ################
    @abstractmethod
    def build_spline(
            self, 
            kind: str, 
            fill_value: Union[str,numpy.ndarray], 
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def has_spline(
            self
        ) -> bool:
        """Return whether or not the GeneticMap has a built spline."""
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def interp_gmap(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            **kwargs: dict
        ) -> 'GeneticMap':
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
        out : GeneticMap
            An interpolated genetic map sharing a copy of the spline from the 
            original genetic map.
        """
        raise NotImplementedError("method is abstract")

    ############### Genetic Distance Methods ###############
    @abstractmethod
    def gdist1g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            ast: Optional[Integral], 
            asp: Optional[Integral]
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gdist2g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            rst: Optional[Integral], 
            rsp: Optional[Integral], 
            cst: Optional[Integral], 
            csp: Optional[Integral]
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gdist1p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            ast: Optional[Integral], 
            asp: Optional[Integral]
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gdist2p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            rst: Optional[Integral], 
            rsp: Optional[Integral], 
            cst: Optional[Integral], 
            csp: Optional[Integral]
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
        raise NotImplementedError("method is abstract")

    #################### Export Methods ####################
    @abstractmethod
    def to_pandas(
            self, 
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a GeneticMap object to a pandas.DataFrame.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            pandas.DataFrame.
        
        Returns
        -------
        out : pandas.DataFrame
            An output dataframe.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def to_csv(
            self, 
            filename: str,
            **kwargs: dict
        ) -> None:
        """
        Write a GeneticMap to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    #################### Import Methods ####################
    @classmethod
    @abstractmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame,
            **kwargs: dict
        ) -> 'GeneticMap':
        """
        Read an object from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            pandas.DataFrame.

        Returns
        -------
        out : GeneticMap
            A GeneticMap read from a pandas.DataFrame.
        """
        raise NotImplementedError("class method is abstract")

    @classmethod
    @abstractmethod
    def from_csv(
            cls, 
            filename: str,
            **kwargs: dict
        ) -> 'GeneticMap':
        """
        Read a GeneticMap from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : GeneticMap
            A GeneticMap read from a CSV file.
        """
        raise NotImplementedError("class method is abstract")


################################## Utilities ###################################
def check_is_GeneticMap(v: object, vname: str) -> None:
    """
    Check if object is of type GeneticMap. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GeneticMap):
        raise TypeError("variable '{0}' must be a GeneticMap".format(vname))
