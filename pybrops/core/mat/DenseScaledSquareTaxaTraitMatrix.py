"""
Module implementing a dense, scaled matrix with taxa axes that are square and
a trait axis which is not square, and associated error checking routines.
"""

import copy
import h5py
from numbers import Real
from pathlib import Path
from typing import Optional, Union
import numpy
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group
from pybrops.core.error.error_value_h5py import check_h5py_File_is_readable
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable
from pybrops.core.mat.DenseScaledMatrix import DenseScaledMatrix
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.core.mat.ScaledSquareTaxaTraitMatrix import ScaledSquareTaxaTraitMatrix
from pybrops.core.util.h5py import h5py_File_write_dict


class DenseScaledSquareTaxaTraitMatrix(
        DenseSquareTaxaTraitMatrix,
        DenseScaledMatrix,
        ScaledSquareTaxaTraitMatrix,
    ):
    """
    A concrete class for dense, scaled matrices with taxa axes that are square 
    and a trait axis which is not square.

    The purpose of this abstract class is to merge the following implementations
    and interfaces:

        1. DenseSquareTaxaTraitMatrix (implementation)
        2. DenseScaledMatrix (implementation)
        3. ScaledSquareTaxaTraitMatrix (interface)
    """

    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self,
            mat: numpy.ndarray, 
            location: Union[numpy.ndarray,Real] = 0.0, 
            scale: Union[numpy.ndarray,Real] = 1.0, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseScaledSquareTaxaTraitMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            Matrix values to store.

        location : numpy.ndarray, Real
            An array of shape ``(t,)`` containing locations for each trait.
            If ``Real``, then the provided location is used for each trait.

        scale : numpy.ndarray, Real
            An array of shape ``(t,)`` containing scales for each trait.
            If ``Real``, then the provided scale is used for each trait.

        taxa : numpy.ndarray, None
            An array of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.

        taxa_grp : numpy.ndarray, None
            An array of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.

        trait : numpy.ndarray, None
            An array of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.

        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call DenseSquareTaxaTraitMatrix constructor
        super(DenseScaledSquareTaxaTraitMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )
        # set location and scale
        self.location = location
        self.scale = scale

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A copy of the DenseScaledSquareTaxaTraitMatrix.
        """
        return self.__class__(
            mat = copy.copy(self.mat),
            location = copy.copy(self.location),
            scale = copy.copy(self.scale),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            trait = copy.copy(self.trait),
        )
    
    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict, None, default = None
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A deep copy of the DenseScaledSquareTaxaTraitMatrix.
        """
        return self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            location = copy.deepcopy(self.location, memo),
            scale = copy.deepcopy(self.scale, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            trait = copy.deepcopy(self.trait, memo),
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a shallow copy of the DenseScaledSquareTaxaTraitMatrix.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A shallow copy of the original DenseScaledSquareTaxaTraitMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseScaledSquareTaxaTraitMatrix':
        """
        Make a deep copy of the DenseScaledSquareTaxaTraitMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseScaledSquareTaxaTraitMatrix
            A deep copy of the original DenseScaledSquareTaxaTraitMatrix.
        """
        return copy.deepcopy(self, memo)

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write DenseSquareTaxaTraitMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseSquareTaxaTraitMatrix`` data is stored.
            If ``None``, the ``DenseSquareTaxaTraitMatrix`` is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
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
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
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
        #### write data to HDF5 file and (optionally) close ####

        # data dictionary
        data = {
            "mat"           : self.mat,
            "location"      : self.location,
            "scale"         : self.scale,
            "taxa"          : self.taxa,
            "taxa_grp"      : self.taxa_grp,
            "trait"         : self.trait,
            # metadata
            "taxa_grp_name" : self.taxa_grp_name,
            "taxa_grp_stix" : self.taxa_grp_stix,
            "taxa_grp_spix" : self.taxa_grp_spix,
            "taxa_grp_len"  : self.taxa_grp_len,
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
            groupname: Optional[str] = None
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Read DenseSquareTaxaTraitMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which DenseSquareTaxaTraitMatrix data is stored.
            If None, DenseSquareTaxaTraitMatrix is read from base HDF5 group.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
            A dense matrix read from file.
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
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = ["mat","location","scale"]

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)

        ########################################################
        ### read data from HDF5 file and (optionally) close ####
        
        # output dictionary
        data = {
            "mat"           : None,
            "location"      : None,
            "scale"         : None,
            "taxa"          : None,
            "taxa_grp"      : None,
            "trait"         : None,
            # metadata
            "taxa_grp_name" : None,
            "taxa_grp_stix" : None,
            "taxa_grp_spix" : None,
            "taxa_grp_len"  : None,
        }
        
        ##################################
        ### read mandatory data fields ###

        # read mat array
        data["mat"] = h5file[groupname + "mat"][()]

        # read location array
        data["location"] = h5file[groupname + "location"][()]

        # read scale array
        data["scale"] = h5file[groupname + "scale"][()]
        
        #################################
        ### read optional data fields ###

        # read taxa data, if "groupname/taxa" in HDF5 file
        if groupname + "taxa" in h5file:
            data["taxa"] = numpy.array([s.decode("utf-8") if isinstance(s,bytes) else s for s in h5file[groupname+"taxa"][()]], dtype=object)

        # read taxa_grp data, if "groupname/taxa_grp" in HDF5 file
        if groupname + "taxa_grp" in h5file:
            data["taxa_grp"] = h5file[groupname + "taxa_grp"][()]
        
        # read trait data, if "groupname/trait" in HDF5 file
        if groupname + "trait" in h5file:
            data["trait"] = numpy.array([s.decode("utf-8") if isinstance(s,bytes) else s for s in h5file[groupname+"trait"][()]], dtype=object)

        #####################################
        ### read optional metadata fields ###

        # read taxa_grp_name data, if "groupname/taxa_grp_name" in HDF5 file
        if groupname + "taxa_grp_name" in h5file:
            data["taxa_grp_name"] = h5file[groupname + "taxa_grp_name"][()]

        # read taxa_grp_stix data, if "groupname/taxa_grp_stix" in HDF5 file
        if groupname + "taxa_grp_stix" in h5file:
            data["taxa_grp_stix"] = h5file[groupname + "taxa_grp_stix"][()]

        # read taxa_grp_spix data, if "groupname/taxa_grp_spix" in HDF5 file
        if groupname + "taxa_grp_spix" in h5file:
            data["taxa_grp_spix"] = h5file[groupname + "taxa_grp_spix"][()]

        # read taxa_grp_len data, if "groupname/taxa_grp_len" in HDF5 file
        if groupname + "taxa_grp_len" in h5file:
            data["taxa_grp_len"] = h5file[groupname + "taxa_grp_len"][()]

        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string or Path an not an h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

        ########################################################
        ################### Object creation ####################
        
        # create object from read data
        mat = cls(
            mat      = data["mat"],
            taxa     = data["taxa"],
            taxa_grp = data["taxa_grp"],
            trait    = data["trait"],
        )

        # assign metadata
        mat.taxa_grp_name = data["taxa_grp_name"]
        mat.taxa_grp_stix = data["taxa_grp_stix"]
        mat.taxa_grp_spix = data["taxa_grp_spix"]
        mat.taxa_grp_len  = data["taxa_grp_len"]

        return mat

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseScaledSquareTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseScaledSquareTaxaTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseScaledSquareTaxaTraitMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseScaledSquareTaxaTraitMatrix.__name__,
                type(v).__name__
            )
        )
