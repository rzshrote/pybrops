# append paths
import sys
import os
model_dir = os.path.dirname(os.path.realpath(__file__))     # get pybropt/test/model
test_dir = os.path.dirname(model_dir)                       # get pybropt/test
pybropt_dir = os.path.dirname(test_dir)                     # get pybropt
sys.path.append(pybropt_dir)                                # append pybropt

# import 3rd party modules we'll need
import unittest
import numpy

# import our libraries
import model

class test_ParametricGenomicModel(unittest.TestCase):
    def test_constructor(self):
        genomic_model = test_ParametricGenomicModel._create_fake_data()

        self.assertIsInstance(genomic_model, model.GenomicModel)

    def _create_fake_data():
        trait = numpy.string_(['trait1','trait2','trait3'])
        coeff = numpy.float64([
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6],
            [0.7, 0.8, 0.9],
            [1.0, 1.1, 1.2]
        ])

        genomic_model = model.ParametricGenomicModel(
            trait,
            coeff
        )
        
        return genomic_model
