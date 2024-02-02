"""
Module defining binary NumPy format I/O interfaces and assocated error 
checking routines.
"""

from abc import ABCMeta, abstractmethod
from typing import Optional, Sequence, Union


class NPYInputOutput(metaclass=ABCMeta):
    """
    Abstract class for defining ``.npy`` file input/output functionality.

    This abstract class defines two functions with the following purposes:
    
    - ``to_npy`` - write an object to a ``.npy`` file.
    - ``from_npy`` - load an object from a ``.npy`` file.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_npy(
            self, 
            filename: Union[str,Sequence[str]],
            fieldname: Union[str,Sequence[str]],
            allow_pickle: bool,
            fix_imports: bool,
            **kwargs: dict
        ) -> None:
        """
        Write an object to a ``.npy`` file.

        Parameters
        ----------
        filename : str, Sequence of str
            ``.npy`` file name(s) to which to write field(s).
        
        fieldname : str, Sequence of str
            Name(s) of field(s) within the object for which to save.
        
        allow_pickle : bool, default = True
            Allow saving object arrays using Python pickles. Reasons for 
            disallowing pickles include security (loading pickled data can 
            execute arbitrary code) and portability (pickled objects may not be 
            loadable on different Python installations, for example if the 
            stored objects require libraries that are not available, and not 
            all pickled data is compatible between Python 2 and Python 3). 
            Default: ``True``

        fix_imports : bool, default = True
            Only useful in forcing objects in object arrays on Python 3 to be 
            pickled in a Python 2 compatible way. If fix_imports is ``True``, 
            pickle will try to map the new Python 3 names to the old module names 
            used in Python 2, so that the pickle data stream is readable with 
            Python 2.
            Default: ``True``

        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            ``.npy`` file.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_npy(
            cls, 
            filename: Union[str,Sequence[str]],
            fieldname: Union[str,Sequence[str]],
            mmap_mode: Optional[str],
            allow_pickle: bool,
            fix_imports: bool,
            encoding: str,
            max_header_size: int,
            **kwargs: dict
        ) -> 'NPYInputOutput':
        """
        Read an object from a ``.npy`` file(s).

        Parameters
        ----------
        filename : str, Sequence of str
            ``.npy`` file name(s) from which to load field(s).
        
        fieldname : str, Sequence of str
            Name(s) of field(s) within the object for which to load.

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
            Allow loading pickled object arrays stored in ``.npy`` files. 
            Reasons for disallowing pickles include security, as loading pickled 
            data can execute arbitrary code. If pickles are disallowed, loading 
            object arrays will fail. 
            Default: ``False``

        fix_imports : bool, default = True
            Only useful when loading Python 2 generated pickled files on 
            Python 3, which includes ``.npy``/``.npz`` files containing object 
            arrays. If ``fix_imports`` is ``True``, pickle will try to map the 
            old Python 2 names to the new names used in Python 3.
            Default: ``True``
        
        encoding : str, default = "ASCII"
            What encoding to use when reading Python 2 strings. Only useful when 
            loading Python 2 generated pickled files in Python 3, which includes 
            ``.npy``/``.npz`` files containing object arrays. Values other than 
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
            ``.npy`` file(s).

        Returns
        -------
        out : NPYInputOutput
            An object read from a ``.npy`` file(s).
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_NPYInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type NPYInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, NPYInputOutput):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,NPYInputOutput.__name__,type(v).__name__))
