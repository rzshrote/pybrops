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
import popgen

class test_Cross(unittest.TestCase):
    
