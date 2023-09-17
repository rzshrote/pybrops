"""
Module defining basal genetic map interfaces and associated error checking routines.
"""

__all__ = [
    "GeneticMap",
    "check_is_GeneticMap",
]

from abc import ABCMeta, abstractmethod
from typing import Optional, Union
import numpy
import pandas

class GeneticMap(metaclass=ABCMeta):
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

    ############################ Object Properties #############################

    ################### Data Properites ####################
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
        Reorder markers in the GeneticMap using an array of indices.
        Note this modifies the GeneticMap in-place.

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

    @abstractmethod
    def prune(
            self, 
            nt: int, 
            M: float, 
            **kwargs: dict
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
            based on their genetic distance exceeds ``nt`` (if provided), the
            additional markers are sought between those regions.
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
        If not grouped, will group.
        Assess physical and genetic map site congruency.

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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def is_congruent(
            self
        ) -> bool:
        """
        Determine if the genetic map is congruent
        Determine if all sites in the genetic map demonstrate congruence with
        their supposed physical and genetic positions.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    ################# Interpolation Methods ################
    @abstractmethod
    def build_spline(
            self, 
            kind: str, 
            fill_value: object, 
            **kwargs: dict
        ) -> None:
        """
        Build a spline for estimating genetic map distances. This is built
        using the marker start indices (self.chr_start)

        Parameters
        ----------
        kind : str
            Specifies the kind of interpolation as a string.
        fill_value : obj
            Fill value for points extrapolated outside the spline.
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
        Interpolate genetic positions given variant physical positions

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
        vrnt_phypos : numpy.ndarray

        Returns
        -------
        out : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def interp_gmap(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            **kwargs: dict
        ):
        """
        Interpolate a new genetic map from the current genetic map.
        Associate spline of current GeneticMap with new GeneticMap.
        """
        raise NotImplementedError("method is abstract")

    ############### Genetic Distance Methods ###############
    @abstractmethod
    def gdist1g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            ast: Optional[int], 
            asp: Optional[int]
        ) -> numpy.ndarray:
        """
        Calculate sequential genetic distances using genetic positions.
        Requires vrnt_chrgrp and vrnt_genpos to have been sorted descending.

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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gdist2g(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            rst: Optional[int], 
            rsp: Optional[int], 
            cst: Optional[int], 
            csp: Optional[int]
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
            A 2D array of pairwise distances between markers.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gdist1p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            ast: Optional[int], 
            asp: Optional[int]
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
        ast : int, None
            Optional array start index (inclusive).
        asp : int, None
            Optional array stop index (exclusive).

        Returns
        -------
        gdist : numpy.ndarray
            A 1D array of distances between the marker prior.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gdist2p(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            rst: Optional[int], 
            rsp: Optional[int], 
            cst: Optional[int], 
            csp: Optional[int]
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
            A 2D array of pairwise distances between markers.
        """
        raise NotImplementedError("method is abstract")

    #################### Export Methods ####################
    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def to_csv(
            self, 
            fname: str, 
            sep: str, 
            header: bool, 
            index: Union[bool,int], 
            **kwargs: dict
        ) -> None:
        """
        Convert a GeneticMap object to a csv file.

        Parameters
        ----------
        fname : str
        sep : str
        header : bool
        index : bool, int
        kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")



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
