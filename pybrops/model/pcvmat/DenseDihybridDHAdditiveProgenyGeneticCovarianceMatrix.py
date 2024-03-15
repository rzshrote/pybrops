"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genetic covariance estimates calculated using dihybrid DH
formulae.
"""

import math
from numbers import Integral, Real
from pathlib import Path
from typing import Optional, Union
import numpy
import pandas
import h5py
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_Integral, check_is_Integral_or_None, check_is_Integral_or_inf, check_is_str, check_is_str_or_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column, check_pandas_DataFrame_has_column_index
from pybrops.core.util.arrayix import flattenix
from pybrops.core.util.subroutines import srange
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.vmat.util import cov_D1s
from pybrops.model.vmat.util import cov_D2s
from pybrops.model.pcvmat.DenseAdditiveProgenyGeneticCovarianceMatrix import DenseAdditiveProgenyGeneticCovarianceMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction, check_is_GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix(DenseAdditiveProgenyGeneticCovarianceMatrix):
    """
    A concrete class for dense additive genetic covariance matrices calculated
    for dihybrid DH progenies.

    The purpose of this concrete class is to implement functionality for:
        1) Genetic covariance estimation for dihybrid DH progenies.
        2) I/O for dihybrid DH progeny covariance matrices.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @DenseAdditiveProgenyGeneticCovarianceMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 4) # (ntaxa,ntaxa,ntrait,ntrait)
        self._mat = value

    ############## Square Metadata Properties ##############
    @DenseAdditiveProgenyGeneticCovarianceMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        return (0,1) # (female, male)

    @DenseAdditiveProgenyGeneticCovarianceMatrix.square_trait_axes.getter
    def square_trait_axes(self) -> tuple:
        return (2,3) # (trait1, trait2) covariance matrix

    #################### Trait metadata ####################

    ######## Expected parental genome contributions ########
    @DenseAdditiveProgenyGeneticCovarianceMatrix.epgc.getter
    def epgc(self) -> tuple:
        """Get a tuple of the expected parental genome contributions."""
        return (0.5, 0.5)

    ################# Parental dimensions ##################
    @property
    def nfemale(self) -> Integral:
        """Number of female parents."""
        return self._mat.shape[self.female_axis]
    
    @property
    def female_axis(self) -> Integral:
        """Axis along which female parents are stored."""
        return 0
    
    @property
    def nmale(self) -> Integral:
        """Number of male parents."""
        return self._mat.shape[self.male_axis]
    
    @property
    def male_axis(self) -> Integral:
        """Axis along which male parents are stored."""
        return 1
    
    ############################## Object Methods ##############################

    ###################### Matrix I/O ######################
    # TODO: make exporting of specific female, male, trait values, not just all
    def to_pandas(
            self, 
            female_col: str = "female",
            female_grp_col: Optional[str] = "female_grp",
            male_col: str = "male",
            male_grp_col: Optional[str] = "male_grp",
            trait1_col: str = "trait1",
            trait2_col: str = "trait2",
            covariance_col: str = "covariance",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix to a pandas.DataFrame.

        Parameters
        ----------
        female_col : str, default = "female"
            Name of the column to which to write female taxa names.

        female_grp_col : str, None, default = "female_grp"
            Name of the column to which to write female taxa groups.

        male_col : str, default = "male"
            Name of the column to which to write male taxa names.

        male_grp_col : str, None, default = "male_grp"
            Name of the column to which to write male taxa groups.

        trait_col : str, default = "trait"
            Name of the column to which to write trait taxa names.

        variance_col : str, default = "covariance"
            Name of the column to which to write covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            pandas.DataFrame.
        
        Returns
        -------
        out : pandas.DataFrame
            An output dataframe.
        """
        # type checks
        check_is_str(female_col, "female_col")
        if female_grp_col is not None:
            check_is_str(female_grp_col, "female_grp_col")
        check_is_str(male_col, "male_col")
        if male_grp_col is not None:
            check_is_str(male_grp_col, "male_grp_col")
        check_is_str(trait1_col, "trait1_col")
        check_is_str(trait2_col, "trait2_col")
        check_is_str(covariance_col, "variance_col")

        # get names for taxa
        if self.taxa is None:
            nzero = math.ceil(math.log10(self.ntaxa))+1
            taxa = numpy.array(["Taxon"+str(e).zfill(nzero) for e in range(self.ntaxa)], dtype = object)
        else:
            taxa = self.taxa

        # get names for traits
        if self.trait is None:
            nzero = math.ceil(math.log10(self.ntrait))+1
            trait = numpy.array(["Trait"+str(e).zfill(nzero) for e in range(self.ntrait)], dtype = object)
        else:
            trait = self.trait

        # calculate flattened array and corresponding axis indices
        flatmat, (femaleix, maleix, trait1ix, trait2ix) = flattenix(self.mat)

        # make dictionary to store output columns in specific column ordering
        out_dict = {}
        out_dict.update({female_col: taxa[femaleix]})
        if female_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[femaleix]
            out_dict.update({female_grp_col: values})
        out_dict.update({male_col: taxa[maleix]})
        if male_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[maleix]
            out_dict.update({male_grp_col: values})
        out_dict.update({trait1_col: trait[trait1ix]})
        out_dict.update({trait2_col: trait[trait2ix]})
        out_dict.update({covariance_col: flatmat})

        # create a pandas DataFrame from the data
        out = pandas.DataFrame(out_dict, **kwargs)

        return out

    def to_csv(
            self, 
            filename: str,
            female_col: str = "female",
            female_grp_col: Optional[str] = "female_grp",
            male_col: str = "male",
            male_grp_col: Optional[str] = "male_grp",
            trait1_col: str = "trait1",
            trait2_col: str = "trait2",
            covariance_col: str = "covariance",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Write a DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        
        female_col : str, default = "female"
            Name of the column to which to write female taxa names.

        female_grp_col : str, None, default = "female_grp"
            Name of the column to which to write female taxa groups.

        male_col : str, default = "male"
            Name of the column to which to write male taxa names.

        male_grp_col : str, None, default = "male_grp"
            Name of the column to which to write male taxa groups.

        trait_col : str, default = "trait"
            Name of the column to which to write trait taxa names.

        variance_col : str, default = "covariance"
            Name of the column to which to write covariance taxa names.

        sep : str, default = ","
            Separator to use in the exported CSV file.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV file.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # convert DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix to pandas.DataFrame
        df = self.to_pandas(
            female_col = female_col,
            female_grp_col = female_grp_col,
            male_col = male_col,
            male_grp_col = male_grp_col,
            trait1_col = trait1_col,
            trait2_col = trait2_col,
            covariance_col = covariance_col,
        )

        # export using pandas
        df.to_csv(
            path_or_buf = filename,
            sep = sep,
            header = header,
            index = index,
            **kwargs
        )

    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            HDF5 file name or HDF5 file stream to which to write.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` data is stored.
            If ``None``, the ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` is written to the base HDF5 group.
        
        overwrite : bool
            Whether to overwrite data fields if they are present in the HDF5 file.
        """
        # call super function
        super(DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix, self).to_hdf5(
            filename  = filename,
            groupname = groupname,
            overwrite = overwrite,
        )
    
    ############################## Class Methods ###############################

    ###################### Matrix I/O ######################
    @classmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame, 
            female_col: Union[str,Integral] = "female",
            female_grp_col: Optional[Union[str,Integral]] = "female_grp",
            male_col: Union[str,Integral] = "male",
            male_grp_col: Optional[Union[str,Integral]] = "male_grp",
            trait1_col: Union[str,Integral] = "trait1",
            trait2_col: Union[str,Integral] = "trait2",
            covariance_col: Union[str,Integral] = "covariance",
            **kwargs: dict
        ) -> 'DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Read a ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` from a ``pandas.DataFrame``.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        female_col : str, Integral, default = "female"
            Name or index of the column from which to read female taxa names.

        female_grp_col : str, Integral, None, default = "female"
            Name or index of the column from which to read female taxa groups.

        male_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        male_grp_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        trait_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        variance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``pandas.DataFrame``.
        
        Returns
        -------
        out : DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix
            A ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` read from a ``pandas.DataFrame``.
        """
        ### type checks and get df indices

        # df
        check_is_pandas_DataFrame(df, "df")

        # female_col
        if isinstance(female_col, str):
            check_pandas_DataFrame_has_column(df, "df", female_col)
            female_colix = df.columns.get_loc(female_col)
        elif isinstance(female_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", female_col)
            female_colix = female_col
        else:
            check_is_str_or_Integral(female_col, "female_col")
        
        # female_grp_col
        if female_grp_col is not None:
            if isinstance(female_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", female_grp_col)
                female_grp_colix = df.columns.get_loc(female_grp_col)
            elif isinstance(female_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", female_grp_col)
                female_grp_colix = female_grp_col
            else:
                check_is_str_or_Integral(female_grp_col, "female_grp_col")
        
        # male_col
        if isinstance(male_col, str):
            check_pandas_DataFrame_has_column(df, "df", male_col)
            male_colix = df.columns.get_loc(male_col)
        elif isinstance(male_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", male_col)
            male_colix = male_col
        else:
            check_is_str_or_Integral(male_col, "male_col")
        
        # male_grp_col
        if male_grp_col is not None:
            if isinstance(male_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", male_grp_col)
                male_grp_colix = df.columns.get_loc(male_grp_col)
            elif isinstance(male_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", male_grp_col)
                male_grp_colix = male_grp_col
            else:
                check_is_str_or_Integral(male_grp_col, "male_grp_col")
        
        # trait1_col
        if isinstance(trait1_col, str):
            check_pandas_DataFrame_has_column(df, "df", trait1_col)
            trait1_colix = df.columns.get_loc(trait1_col)
        elif isinstance(trait1_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", trait1_col)
            trait1_colix = trait1_col
        else:
            check_is_str_or_Integral(trait1_col, "trait1_col")
        
        # trait2_col
        if isinstance(trait2_col, str):
            check_pandas_DataFrame_has_column(df, "df", trait2_col)
            trait2_colix = df.columns.get_loc(trait2_col)
        elif isinstance(trait2_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", trait2_col)
            trait2_colix = trait2_col
        else:
            check_is_str_or_Integral(trait2_col, "trait2_col")
        
        # covariance_col
        if isinstance(covariance_col, str):
            check_pandas_DataFrame_has_column(df, "df", covariance_col)
            covariance_colix = df.columns.get_loc(covariance_col)
        elif isinstance(covariance_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", covariance_col)
            covariance_colix = covariance_col
        else:
            check_is_str_or_Integral(covariance_col, "variance_col")
        
        # get required data columns (type numpy.ndarray)
        female_data     = df.iloc[:,female_colix    ].to_numpy(dtype = object)
        male_data       = df.iloc[:,male_colix      ].to_numpy(dtype = object)
        trait1_data     = df.iloc[:,trait1_colix    ].to_numpy(dtype = object)
        trait2_data     = df.iloc[:,trait2_colix    ].to_numpy(dtype = object)
        covariance_data = df.iloc[:,covariance_colix].to_numpy(dtype = float)

        # get unique female, male taxa (type numpy.ndarray)
        female_taxa, female_taxaix = numpy.unique(female_data, return_index = True)
        male_taxa, male_taxaix     = numpy.unique(male_data, return_index = True)

        # combine female and male taxa names (type numpy.ndarray)
        taxa = numpy.union1d(female_taxa, male_taxa)

        # allocate female, male indices
        femaleix = numpy.empty(len(female_data), dtype = int)
        maleix   = numpy.empty(len(male_data  ), dtype = int)

        # calculate female, male indices
        for i,taxon in enumerate(taxa):
            femaleix[female_data == taxon] = i
            maleix  [  male_data == taxon] = i
        
        # calculate unique trait values, trait indices
        trait1, trait1ix = numpy.unique(trait1_data, return_inverse = True)
        trait2, trait2ix = numpy.unique(trait2_data, return_inverse = True)

        # get trait names
        trait = numpy.union1d(trait1, trait2)

        # get optional taxa group data
        taxa_grp = None
        if female_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            female_grp_data = df.iloc[:,female_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in female_taxa:
                    ix = female_taxaix[female_taxa == taxon][0]
                    taxa_grp[i] = female_grp_data[ix]
        if male_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            male_grp_data = df.iloc[:,male_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in male_taxa:
                    ix = male_taxaix[male_taxa == taxon][0]
                    taxa_grp[i] = male_grp_data[ix]

        # get array dimensions
        nfemale = len(taxa)
        nmale = len(taxa)
        ntrait = len(trait)

        # allocate NaN array for covariance matrix
        mat = numpy.full((nfemale,nmale,ntrait,ntrait), numpy.nan, dtype = float)

        # overwrite NaN values with covariance values
        mat[femaleix,maleix,trait1ix,trait2ix] = covariance_data

        # construct an object
        out = cls(
            mat = mat, 
            taxa = taxa, 
            taxa_grp = taxa_grp, 
            trait = trait, 
        )

        return out

    @classmethod
    def from_csv(
            cls, 
            filename: str, 
            female_col: Union[str,Integral] = "female",
            female_grp_col: Optional[Union[str,Integral]] = "female_grp",
            male_col: Union[str,Integral] = "male",
            male_grp_col: Optional[Union[str,Integral]] = "male_grp",
            trait1_col: Union[str,Integral] = "trait1",
            trait2_col: Union[str,Integral] = "trait2",
            covariance_col: Union[str,Integral] = "covariance",
            sep: str = ',',
            header: int = 0,
            **kwargs: dict
        ) -> 'DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Read an object from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.

        female_col : str, Integral, default = "female"
            Name or index of the column from which to read female taxa names.

        female_grp_col : str, Integral, None, default = "female"
            Name or index of the column from which to read female taxa groups.

        male_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        male_grp_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        trait_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        variance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : CSVInputOutput
            An object read from a CSV file.
        """
        df = pandas.read_csv(
            filepath_or_buffer = filename,
            sep = sep,
            header = header,
            **kwargs
        )

        out = cls.from_pandas(
            df = df,
            female_col = female_col,
            female_grp_col = female_grp_col,
            male_col = male_col,
            male_grp_col = male_grp_col,
            trait1_col = trait1_col,
            trait2_col = trait2_col,
            covariance_col = covariance_col,
        )

        return out
    
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Read ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` data is stored.
            If ``None``, ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` is read from base HDF5 group.

        Returns
        -------
        out : DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix
            A ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix`` read from file.
        """
        # call super function
        return super(DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix, cls).from_hdf5(
            filename  = filename,
            groupname = groupname,
        )

    ############# Matrix Factory Class Methods #############
    # TODO: provide support for non-linear models
    @classmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Union[Integral,Real],
            gmapfn: GeneticMapFunction,
            **kwargs: dict
        ) -> 'DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Calculate a symmetrical tensor of progeny variances for each possible
        3-way cross between *inbred* individuals.

        Parameters
        ----------
        gmod : GenomicModel
            Genomic Model with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic covariance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            covariance.
        nself : Integral, Real
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-----------------+-------------------------+
            | Example         | Description             |
            +=================+=========================+
            | ``nself = 0``   | Derive gametes from F1  |
            +-----------------+-------------------------+
            | ``nself = 1``   | Derive gametes from F2  |
            +-----------------+-------------------------+
            | ``nself = 2``   | Derive gametes from F3  |
            +-----------------+-------------------------+
            | ``...``         | etc.                    |
            +-----------------+-------------------------+
            | ``nself = inf`` | Derive gametes from SSD |
            +-----------------+-------------------------+
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix
            A matrix of additive genetic covariance estimations.
        """
        # type checks
        check_is_GenomicModel(gmod, "gmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral_or_inf(nself, "nself")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        
        # if genomic model is an additive linear genomic model, then use specialized routine
        if isinstance(gmod, AdditiveLinearGenomicModel):
            return cls.from_algmod(
                algmod = gmod, 
                pgmat = pgmat, 
                ncross = ncross, 
                nprogeny = nprogeny, 
                nself = nself, 
                gmapfn = gmapfn, 
                **kwargs
            )
        # otherwise raise error since non-linear support hasn't been implemented yet
        else:
            raise NotImplementedError("support for non-linear models not implemented yet")

    @classmethod
    def from_algmod(
            cls, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Union[Integral,Real], 
            gmapfn: Integral, 
            mem: Union[Integral,None] = 1024
        ) -> 'DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Calculate a symmetrical matrix of DH progeny variances for each possible
        pairwise dihybrid cross between *non-inbred* individuals.
        Calculations are derived from Lehermeier et al. (2019).

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic covariance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            covariance.
        nself : Integral, Real
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-----------------+-------------------------+
            | Example         | Description             |
            +=================+=========================+
            | ``nself = 0``   | Derive gametes from F1  |
            +-----------------+-------------------------+
            | ``nself = 1``   | Derive gametes from F2  |
            +-----------------+-------------------------+
            | ``nself = 2``   | Derive gametes from F3  |
            +-----------------+-------------------------+
            | ``...``         | etc.                    |
            +-----------------+-------------------------+
            | ``nself = inf`` | Derive gametes from SSD |
            +-----------------+-------------------------+
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : Integral, default = 1024
            Memory chunk size to use during matrix operations. If ``None``,
            then memory chunk size is not limited.

            WARNING: Setting ``mem = None`` might result in memory allocation
            errors! For reference, ``mem = 1024`` refers to a matrix of size
            1024x1024, which needs about 8.5 MB of storage. Matrices of course
            need a quadratic amount of memory: :math:`O(n^2)`.

        Returns
        -------
        out : ProgenyGeneticCovarianceMatrix
            A matrix of additive genetic covariance estimations.
        """
        # type checks
        check_is_AdditiveLinearGenomicModel(algmod, "algmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral_or_inf(nself, "nself")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_Integral_or_None(mem, "mem")

        # check for chromosome grouping
        if not pgmat.is_grouped_vrnt():
            raise ValueError("pgmat must be grouped along the vrnt axis")
        
        # check for genetic positions
        if pgmat.vrnt_genpos is None:
            raise ValueError("pgmat must have genetic positions")

        # gather shapes of data input
        ntrait = algmod.ntrait                  # number ot traits (t)
        ntaxa = pgmat.ntaxa                     # number of individuals (n)

        # gather data pointers
        geno = pgmat.mat                        # (m,n,p) genotype matrix
        chrgrp_stix = pgmat.vrnt_chrgrp_stix    # (g,) chromosome group start indices
        chrgrp_spix = pgmat.vrnt_chrgrp_spix    # (g,) chromosome group stop indices
        genpos = pgmat.vrnt_genpos              # (p,) marker genetic positions
        u = algmod.u_a                          # (p,t) marker effect coefficients

        # allocate a square matrix for each pairwise covariance
        varA_mat = numpy.zeros(
            (ntaxa, ntaxa, ntrait),     # (n,n,t) covariance matrix
            dtype = float
        )

        # for each linkage group
        for lst, lsp in zip(chrgrp_stix, chrgrp_spix):
            # determine memory chunk step (for iterators)
            step = (lsp - lst) if mem is None else mem
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                for cst,csp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                    # create sparse meshgrid indicating where genetic positions are
                    gi, gj = numpy.meshgrid(
                        genpos[rst:rsp],        # row block genetic positions
                        genpos[cst:csp],        # column block genetic positions
                        indexing = 'ij',        # use ij indexing
                        sparse = True           # generate a spare array tuple for speed
                    )

                    # calculate recombination probability matrix for chunk
                    r = gmapfn.mapfn(numpy.abs(gi - gj)) # (rb,cb) covariance matrix

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1 = cov_D1s(r, nself)  # (rb,cb) covariance matrix

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2 = cov_D2s(r, nself)  # (rb,cb) covariance matrix

                    # get marker coefficients for rows and columns
                    ru = u[rst:rsp].T # (rb,t)' -> (t,rb)
                    cu = u[cst:csp].T # (cb,t)' -> (t,cb)

                    # for each mate pair (including selfs)
                    # subscript codes:
                    #   1 = female phase 2
                    #   2 = female phase 1
                    #   3 = male phase 2
                    #   4 = male phase 1
                    for female in range(0,ntaxa):       # varA row index
                        for male in range(0,female):    # varA col index
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,female,rst:rsp] - geno[1,female,rst:rsp] # (rb,)
                            cdgeno21 = geno[0,female,cst:csp] - geno[1,female,cst:csp] # (cb,)
                            rdgeno31 = geno[1,male,rst:rsp] - geno[1,female,rst:rsp] # (rb,)
                            cdgeno31 = geno[1,male,cst:csp] - geno[1,female,cst:csp] # (cb,)
                            rdgeno32 = geno[1,male,rst:rsp] - geno[0,female,rst:rsp] # (rb,)
                            cdgeno32 = geno[1,male,cst:csp] - geno[0,female,cst:csp] # (cb,)
                            rdgeno41 = geno[0,male,rst:rsp] - geno[1,female,rst:rsp] # (rb,)
                            cdgeno41 = geno[0,male,cst:csp] - geno[1,female,cst:csp] # (cb,)
                            rdgeno42 = geno[0,male,rst:rsp] - geno[0,female,rst:rsp] # (rb,)
                            cdgeno42 = geno[0,male,cst:csp] - geno[0,female,cst:csp] # (cb,)
                            rdgeno43 = geno[0,male,rst:rsp] - geno[1,male,rst:rsp] # (rb,)
                            cdgeno43 = geno[0,male,cst:csp] - geno[1,male,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect21 = rdgeno21 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect21 = cdgeno21 * cu # (cb,)*(t,cb) -> (t,cb)
                            reffect31 = rdgeno31 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect31 = cdgeno31 * cu # (cb,)*(t,cb) -> (t,cb)
                            reffect32 = rdgeno32 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect32 = cdgeno32 * cu # (cb,)*(t,cb) -> (t,cb)
                            reffect41 = rdgeno41 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect41 = cdgeno41 * cu # (cb,)*(t,cb) -> (t,cb)
                            reffect42 = rdgeno42 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect42 = cdgeno42 * cu # (cb,)*(t,cb) -> (t,cb)
                            reffect43 = rdgeno43 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect43 = cdgeno43 * cu # (cb,)*(t,cb) -> (t,cb)

                            # calculate varA part for female-male
                            # (t,rb)@(rb,cb)@(cb,t) -> (t,t)
                            varA_part21 = reffect21 @ D2 @ ceffect21.T
                            varA_part31 = reffect31 @ D1 @ ceffect31.T
                            varA_part32 = reffect32 @ D1 @ ceffect32.T
                            varA_part41 = reffect41 @ D1 @ ceffect41.T
                            varA_part42 = reffect42 @ D1 @ ceffect42.T
                            varA_part43 = reffect43 @ D2 @ ceffect43.T

                            # calculate varA part for this matrix chunk
                            # (t,t) + (t,t) + (t,t) + (t,t) + (t,t) + (t,t) -> (t,t)
                            varA_part = varA_part21 + varA_part31 + varA_part32 + varA_part41 + varA_part42 + varA_part43

                            # add this partial covariance to the lower triangle
                            # (1,1,t,t) + (t,t) -> (1,1,t,t)
                            varA_mat[female,male,:] += varA_part

        # divide entire matrix by 4 to get covariance per the equation
        # multiplication is faster computationally
        varA_mat *= 0.25

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female in range(1, ntaxa):
            for male in range(0, female):
                varA_mat[male,female,:,:] = varA_mat[female,male,:,:]

        # construct output
        out = cls(
            mat = varA_mat,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp
        )

        return out



################################## Utilities ###################################
def check_is_DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix):
        raise TypeError("'{0}' must be a DenseDihybridDHAdditiveProgenyGeneticCovarianceMatrix".format(vname))
