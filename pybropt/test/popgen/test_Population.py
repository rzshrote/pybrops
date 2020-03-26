# append paths
import sys
import os
popgen_dir = os.path.dirname(os.path.realpath(__file__))    # get pybropt/test/popgen
test_dir = os.path.dirname(popgen_dir)                      # get pybropt/test
pybropt_dir = os.path.dirname(test_dir)                     # get pybropt
sys.path.append(pybropt_dir)                                # append pybropt

# import 3rd party modules we'll need
import unittest
import numpy

# import our libraries
from popgen.Population import Population

class test_Population(unittest.TestCase):
    def test_from_vcf_output_type(self):
        # get test data file path
        data_path = test_dir + "/example.vcf"

        # load data
        population = Population.from_vcf(data_path)

        # assert we have a Population type object
        self.assertIsInstance(
            population,
            Population,
            "type(population) = %s" % type(population)
        )

        true_chr_grp = numpy.string_(['20','20','20','20','20'])

        self.assertTrue(
            numpy.all(population.marker_set.chr_grp == true_chr_grp)
        )
