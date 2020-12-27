# append paths
import sys
import os
gmap_dir = os.path.dirname(os.path.realpath(__file__))  # /pybropt/test/popgen/gmap
popgen_dir = os.path.dirname(gmap_dir)                  # /pybropt/test/popgen
test_dir = os.path.dirname(popgen_dir)                  # /pybropt/test
pybropt_dir = os.path.dirname(test_dir)                 # /pybropt
root_dir = os.path.dirname(pybropt_dir)                 # /
sys.path.append(root_dir)                               # append /

# import 3rd party modules we'll need
import unittest
import numpy

# import our libraries
from pybropt.popgen.gmap import ExtendedGeneticMap

class test_ExtendedGeneticMap(unittest.TestCase):
    def test_from_egmap(self):
        # get test data file path
        data_path = gmap_dir + "/McMullen_2009_US_NAM.M.egmap"

        genetic_map = ExtendedGeneticMap.from_egmap(data_path)

        self.assertIsInstance(genetic_map, ExtendedGeneticMap)
