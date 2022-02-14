import os.path

def check_path_exists(path):
    if not os.path.exists(path):
        raise FileNotFoundError("{0} does not exist".format(path))

def check_file_exists(fname):
    if not os.path.isfile(fname):
        raise FileNotFoundError("{0} does not exist".format(fname))

def check_directory_exists(path):
    if not os.path.exists(path):
        raise NotADirectoryError("{0} does not exist".format(path))
