import sys
import pathlib

print(pathlib.Path(__file__).parents[2].resolve().as_posix())
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
