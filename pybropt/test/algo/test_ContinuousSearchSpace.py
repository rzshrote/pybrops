# append paths
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))))

# import 3rd party modules we'll need
import unittest

# import our libraries
import pybropt.algo.ContinuousSearchSpace
