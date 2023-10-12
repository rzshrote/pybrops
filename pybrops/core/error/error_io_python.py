"""
Module containing error subroutines related to generic file and paths.
"""

__all__ = [
    "check_path_exists",
    "check_file_exists",
    "check_directory_exists",
]

from os.path import exists, isfile

def check_path_exists(path: str) -> None:
    """
    Subroutine to check whether a given path exists.
    If the path does not exist, raise an FileNotFoundError with a custom error
    message.

    Parameters
    ----------
    path : str, path-like object
        Path to check.
    """
    if not exists(path):
        raise FileNotFoundError("{0} does not exist".format(path))

def check_file_exists(path: str) -> None:
    """
    Subroutine to check whether a given file exists.
    If the file does not exist, raise an FileNotFoundError with a custom error
    message.

    Parameters
    ----------
    path : str, path-like object
        Path to check.
    """
    if not isfile(path):
        raise FileNotFoundError("{0} does not exist".format(path))

def check_directory_exists(path: str) -> None:
    """
    Subroutine to check whether a given directory exists.
    If the directory does not exist, raise an NotADirectoryError with a custom
    error message.

    Parameters
    ----------
    path : str, path-like object
        Path to check.
    """
    if not exists(path):
        raise NotADirectoryError("{0} does not exist".format(path))
