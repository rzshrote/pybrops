"""
Module defining interfaces and error checking routimes for progeny mean 
estimated breeding value matrices.
"""

__all__ = [
    "ProgenyMeanEstimatedBreedingValueMatrix",
    "check_is_ProgenyMeanEstimatedBreedingValueMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pathlib import Path
from typing import Optional, Union
import h5py
import numpy
import pandas
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.SquareTaxaTraitMatrix import SquareTaxaTraitMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix

class ProgenyMeanEstimatedBreedingValueMatrix(
        SquareTaxaTraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
        HDF5InputOutput,
        metaclass = ABCMeta,
    ):
    """
    Abstract class for progeny mean EBV matrix representation.

    The ProgenyMeanEstimatedBreedingValueMatrix class represents a Multivariate 
    Progeny Mean Estimated Breeding Value.

    Notes
    -----
    All elements within a ProgenyMeanEstimatedBreedingValueMatrix are mean-
    centered and scaled to unit variance for each trait.

    .. math::
        BV = \\frac{X - \\mu}{\\sigma}

    Where:

    - :math:`BV` is the breeding value.
    - :math:`X` is the phenotype value.
    - :math:`\\mu` is the mean (location) for :math:`X`.
    - :math:`\\sigma` is the standard deviation (scale) for :math:`X`.

    Phenotype values can be reconstituted using:

    .. math::
        X = \\sigma BV + \\mu
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ######## Expected parental genome contributions ########
    @property
    @abstractmethod
    def epgc(self) -> tuple:
        """Expected parental genome contribution to the offspring from each parent."""
        raise NotImplementedError("property is abstract")

    ############################ Object Properties #############################
    @property
    @abstractmethod
    def location(self) -> object:
        """Mean of the phenotype values used to calculate breeding values."""
        raise NotImplementedError("property is abstract")
    @location.setter
    @abstractmethod
    def location(self, value: object) -> None:
        """Set the mean of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def scale(self) -> object:
        """Standard deviation of the phenotype values used to calculate breeding values."""
        raise NotImplementedError("property is abstract")
    @scale.setter
    @abstractmethod
    def scale(self, value: object) -> None:
        """Set the standard deviation of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ############## Matrix summary statistics ###############

    @abstractmethod
    def unscale(self) -> numpy.ndarray:
        """
        Transform values within the ProgenyMeanEstimatedBreedingValueMatrix 
        back to their unscaled and de-centered values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n,n,t)`` containing unscaled and de-centered
            values.

            Where:

            - ``n`` is the number of taxa.
            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    #################### Export Methods ####################
    
    # to_pandas (inherited)
    # to_csv (inherited)
    # to_hdf5 (inherited)

    ############################## Class Methods ###############################

    #################### Import Methods ####################

    @classmethod
    @abstractmethod
    def from_numpy(
            cls, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'ProgenyMeanEstimatedBreedingValueMatrix':
        """
        Construct a ProgenyMeanEstimatedBreedingValueMatrix from a numpy.ndarray.
        Calculates mean-centering and scaling to unit variance.

        Parameters
        ----------
        mat : numpy.ndarray
            An array of shape which to mean-center and scale.
        taxa : numpy.ndarray
            An array of taxa names.
        taxa_grp : numpy.ndarray
            An array of taxa groups.
        trait : numpy.ndarray
            An array of trait names.

        Returns
        -------
        out : ProgenyMeanEstimatedBreedingValueMatrix
            Output progeny mean estimated breeding value matrix.
        """
        raise NotImplementedError("class method is abstract")

    @classmethod
    @abstractmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame, 
            **kwargs: dict
        ) -> 'ProgenyMeanEstimatedBreedingValueMatrix':
        """
        Read a ProgenyMeanEstimatedBreedingValueMatrix from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            pandas.DataFrame.

        Returns
        -------
        out : ProgenyMeanEstimatedBreedingValueMatrix
            A ProgenyMeanEstimatedBreedingValueMatrix read from a pandas.DataFrame.
        """
        raise NotImplementedError("class method is abstract")
    
    @classmethod
    @abstractmethod
    def from_csv(
            cls, 
            filename: str, 
            **kwargs: dict
        ) -> 'ProgenyMeanEstimatedBreedingValueMatrix':
        """
        Read a ProgenyMeanEstimatedBreedingValueMatrix from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : ProgenyMeanEstimatedBreedingValueMatrix
            A ProgenyMeanEstimatedBreedingValueMatrix read from a CSV file.
        """
        raise NotImplementedError("class method is abstract")
    
    @classmethod
    @abstractmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str],
        ) -> 'ProgenyMeanEstimatedBreedingValueMatrix':
        """
        Read a ProgenyMeanEstimatedBreedingValueMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, HDF5 group name under which object data is stored.
            If ``None``, object is read from base HDF5 group.

        Returns
        -------
        out : ProgenyMeanEstimatedBreedingValueMatrix
            A ProgenyMeanEstimatedBreedingValueMatrix read from an HDF5 file.
        """
        raise NotImplementedError("class method is abstract")

    ################# Construction Methods #################
    
    @classmethod
    @abstractmethod
    def from_bvmat(
            bvmat: BreedingValueMatrix,
            **kwargs: dict
        ) -> 'ProgenyMeanEstimatedBreedingValueMatrix':
        """
        Calculate progeny mean estimated breeding values from a ``BreedingValueMatrix``.

        Parameters
        ----------
        bvmat : BreedingValueMatrix
            Breeding value matrix from which to estimate progeny mean EBVs.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : ProgenyMeanEstimatedBreedingValueMatrix
            A matrix of progeny mean EBVs.
        """
        raise NotImplementedError("class method is abstract")

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_ProgenyMeanEstimatedBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``ProgenyMeanEstimatedBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, ProgenyMeanEstimatedBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                ProgenyMeanEstimatedBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
