# append paths
import sys
import os
breed_dir = os.path.dirname(os.path.realpath(__file__))     # get pybropt/test/breed
test_dir = os.path.dirname(breed_dir)                       # get pybropt/test
pybropt_dir = os.path.dirname(test_dir)                     # get pybropt
sys.path.append(pybropt_dir)                                # append pybropt

# import 3rd party modules we'll need
import unittest
import numpy

# import our libraries
import breed
import popgen

class test_Breeding(unittest.TestCase):
    def test_constructor(self):
        # get test data file path
        data_path = test_dir + "/example.vcf"

        # load data
        population = popgen.Population.from_vcf(data_path)

        # create cross
        cross = popgen.Cross(
            population,
            "2wayDH",
            False,
            "2wayDH",
            "ctrl",
            "equal"
        )

        # create breeding object
        breeding = breed.Breeding(
            population,
            cross,
            "Test!"
        )

        self.assertIsInstance(breeding, breed.Breeding)
