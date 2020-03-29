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

class test_CGS(unittest.TestCase):
    def test_constructor(self):
        pass
