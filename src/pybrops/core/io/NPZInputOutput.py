"""
Module defining binary NumPy archive format I/O interfaces and assocated error 
checking routines.
"""

from abc import ABCMeta, abstractmethod
from typing import Optional, Sequence, Union


class NPZInputOutput(metaclass=ABCMeta):
    """
    Abstract class for defining NPZ input/output functionality.

    This abstract class defines two functions with the following purposes:
    
    - ``to_npz`` - write an object to a ``.npz`` file.
    - ``from_npz`` - load an object from a ``.npz`` file.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_npz(
            self, 
            filename: str,
            **kwargs: dict
        ) -> None:
        """
        Write an object to a ``.npz`` file.

        Parameters
        ----------
        filename : str
            ``.npz`` file name to which to write object field(s).
        
        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            ``.npz`` file.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_npz(
            cls, 
            filename: str,
            mmap_mode: Optional[str],
            allow_pickle: bool,
            fix_imports: bool,
            encoding: str,
            max_header_size: int,
            **kwargs: dict
        ) -> 'NPZInputOutput':
        """
        Read an object from a ``.npz`` file.

        Parameters
        ----------
        filename : str, Sequence of str
            ``.npz`` file name from which to load object field(s).
        
        mmap_mode : str, None, default = None
            Options are ``{None, "r+", "r", "w+", "c"}``
            If not ``None``, then memory-map the file, using the given mode 
            (see ``numpy.memmap`` for a detailed description of the modes). A 
            memory-mapped array is kept on disk. However, it can be accessed and 
            sliced like any ``ndarray``. Memory mapping is especially useful for 
            accessing small fragments of large files without reading the entire 
            file into memory.
            Default: ``None``

        allow_pickle : bool, default = False
            Allow loading pickled object arrays stored in ``.npz`` files. 
            Reasons for disallowing pickles include security, as loading pickled 
            data can execute arbitrary code. If pickles are disallowed, loading 
            object arrays will fail. 
            Default: ``False``

        fix_imports : bool, default = True
            Only useful when loading Python 2 generated pickled files on 
            Python 3, which includes ``.npz``/``.npz`` files containing object 
            arrays. If ``fix_imports`` is ``True``, pickle will try to map the 
            old Python 2 names to the new names used in Python 3.
            Default: ``True``
        
        encoding : str, default = "ASCII"
            What encoding to use when reading Python 2 strings. Only useful when 
            loading Python 2 generated pickled files in Python 3, which includes 
            ``.npz``/``.npz`` files containing object arrays. Values other than 
            ``"latin1"``, ``"ASCII"``, and ``"bytes"`` are not allowed, as they 
            can corrupt numerical data.
            Default: ``"ASCII"``

        max_header_size : int, optional
            Maximum allowed size of the header. Large headers may not be safe to 
            load securely and thus require explicitly passing a larger value. See 
            ``ast.literal_eval`` for details. This option is ignored when 
            ``allow_pickle`` is passed. In that case the file is by definition 
            trusted and the limit is unnecessary.
            Default: ``10000``
        
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``.npz`` file(s).

        Returns
        -------
        out : NPZInputOutput
            An object read from a ``.npz`` file.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_NPZInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type NPZInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, NPZInputOutput):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,NPZInputOutput.__name__,type(v).__name__))
