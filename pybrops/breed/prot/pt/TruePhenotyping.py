"""
Module implementing phenotyping protocols for extracting true breeding values.
"""

import copy
import math
from numbers import Real
from pathlib import Path
from typing import Optional, Union
import numpy
import pandas
import h5py

from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group, check_h5py_File_is_readable, check_h5py_File_is_writable
from pybrops.core.util.h5py import h5py_File_write_dict
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class TruePhenotyping(
        PhenotypingProtocol
    ):
    """
    Class implementing phenotyping protocols for extracting true breeding values.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            gpmod: GenomicModel, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class TruePhenotyping.

        Parameters
        ----------
        gpmod : GenomicModel
            Genomic prediction model to use to determine phenotypes.
        kwargs : dict
            Additional keyword arguments
        """
        self.gpmod = gpmod

    def __copy__(
            self
        ) -> 'TruePhenotyping':
        """
        Make a shallow copy of the ``TruePhenotyping`` object.

        Returns
        -------
        out : TruePhenotyping
            A shallow copy of the ``TruePhenotyping`` object.
        """
        # get the class
        cls = type(self)

        # create a new class object
        out = cls(
            gpmod = copy.copy(self.gpmod),
        )

        return out

    def __deepcopy__(
            self,
            memo: Optional[dict] = None,
        ) -> 'TruePhenotyping':
        """
        Make a deep copy of the ``TruePhenotyping`` object.

        Parameters
        ----------
        memo : dict, None
            An optional dictionary of memo metadata.

        Returns
        -------
        out : TruePhenotyping
            A deep copy of the ``TruePhenotyping`` object.
        """
        # get the class
        cls = type(self)

        # create a new class object
        out = cls(
            gpmod = copy.deepcopy(self.gpmod, memo),
        )

        return out

    ############################ Object Properties #############################

    ############### Genomic Model Properties ###############
    @property
    def gpmod(self) -> GenomicModel:
        """Genomic prediction model."""
        return self._gpmod
    @gpmod.setter
    def gpmod(self, value: GenomicModel) -> None:
        """Set genomic prediction model"""
        check_is_GenomicModel(value, "gpmod")
        self._gpmod = value

    ################ Stochastic Parameters #################
    @property
    def var_err(self) -> numpy.ndarray:
        """Error variance for each trait."""
        return numpy.repeat(0.0, self.gpmod.ntrait)
    @var_err.setter
    def var_err(self, value: numpy.ndarray) -> None:
        """Set error variance"""
        error_readonly("var_err")

    ############################## Object Methods ##############################
    def copy(
            self
        ) -> 'TruePhenotyping':
        """
        Make a shallow copy of the ``TruePhenotyping`` object.

        Returns
        -------
        out : TruePhenotyping
            A shallow copy of the ``TruePhenotyping`` object.
        """
        return copy.copy(self)

    def deepcopy(
            self,
            memo: Optional[dict] = None,
        ) -> 'TruePhenotyping':
        """
        Make a deep copy of the ``TruePhenotyping`` object.

        Parameters
        ----------
        memo : dict, None
            An optional dictionary of memo metadata.

        Returns
        -------
        out : TruePhenotyping
            A deep copy of the ``TruePhenotyping`` object.
        """
        return copy.deepcopy(self, memo)
    
    def phenotype(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        gpmod : GenomicModel, None
            Genomic prediction model to use to determine phenotypes.
            If None, use default genomic prediction model.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : pandas.DataFrame
            A pandas.DataFrame containing phenotypes for individuals.
        """
        # check argument data types
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        if miscout is not None:
            check_is_dict(miscout, "miscout")

        # calculate true breeding values
        gvmat = self.gpmod.gegv(pgmat)

        # construct dataframe labels
        labels_dict = {}
        taxazfill = math.ceil(math.log10(gvmat.ntaxa))+1
        labels_dict["taxa"] = ["Taxon"+str(i+1).zfill(taxazfill) for i in range(gvmat.ntaxa)] if gvmat.taxa is None else gvmat.taxa
        if gvmat.taxa_grp is not None:
            labels_dict["taxa_grp"] = gvmat.taxa_grp
        labels_df = pandas.DataFrame(labels_dict)

        # construct dataframe data values
        mat = gvmat.unscale()   # calculate unscaled breeding values
        traitzfill = math.ceil(math.log10(gvmat.ntrait))+1
        cols = ["Trait"+str(i+1).zfill(traitzfill) for i in range(gvmat.ntrait)] if gvmat.trait is None else gvmat.trait
        values_df = pandas.DataFrame(
            data = mat,
            columns = cols
        )

        # combine data labels and values
        out_df = pandas.concat([labels_df, values_df], axis = 1)

        return out_df

    def set_h2(
            self, 
            h2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise AttributeError("unsupported operation: heritability always set at 1.0")

    def set_H2(
            self, 
            H2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : float, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise AttributeError("unsupported operation: heritability always set at 1.0")

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write a ``TruePhenotyping`` object to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an already opened HDF5 file to which to write. File remains open after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which object data is stored.
            If ``None``, object is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            # empty string
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ############ process ``filename`` argument #############

        # h5 file object
        h5file = None
        
        # if we have a string or Path, open HDF5 file in append (``r+``) mode
        if isinstance(filename, (str,Path)):
            h5file = h5py.File(filename, "a")

        # elif we have an h5py.File, make sure mode is writable, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_writable(filename)
            h5file = filename

        # else raise TypeError
        else:
            raise TypeError(
                "filename must be of type ``str`` or ``h5py.File`` but received type {0}".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ### populate HDF5 file

        # data dictionary
        data = {
        }

        # save data
        h5py_File_write_dict(h5file, groupname, data, overwrite)

        # close the file, only if the provided filename was a string or Path and not a h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            gpmod : GenomicModel = None,
        ) -> 'TruePhenotyping':
        """
        Read an object from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which object data is stored.
            If ``None``, object is read from base HDF5 group.
        gpmod : GenomicModel
            A genomic model to bind to the ``TruePhenotyping`` protocol.
            This will be eliminated when a better storage mechanism is available.

        Returns
        -------
        out : TruePhenotyping
            An object read from an HDF5 file.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r``) mode
        if isinstance(filename, (str,Path)):
            check_file_exists(filename)
            h5file = h5py.File(filename, "r")

        # elif we have an h5py.File, make sure mode is in at least ``r`` mode, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_readable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # FIXME: errors if groupname == "" or "/"
            # if the group does not exist in the file, close and raise error
            check_h5py_File_has_group(h5file, groupname)

            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ################## process ``gpmod`` ###################

        check_is_GenomicModel(gpmod, "gpmod")

        ########################################################
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = []

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)
        
        ########################################################
        ### read data from HDF5 file and (optionally) close ####
        
        # output dictionary
        data = {
        }

        ##################################
        ### read mandatory data fields ###

        #################################
        ### read optional data fields ###

        #####################################
        ### read optional metadata fields ###

        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string or Path an not an h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

        ########################################################
        ### create object
        
        # create object from read data
        out = cls(
            gpmod = gpmod,
        )

        return out
