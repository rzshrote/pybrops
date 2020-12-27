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
from popgen.MarkerSet import MarkerSet

class test_MarkerSet(unittest.TestCase):
    """docstring for test_MarkerSet."""
    def test_constructor(self):
        # make fake data
        m1 = test_MarkerSet._create_fake_data()

        sort_chr_grp = numpy.string_(["1","1","1","2","2","2","3","3","3"])
        sort_chr_start = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_chr_stop = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_mkr_name = numpy.string_(['a','b','c','f','e','d','h','i','g'])
        sort_chr_grp_stix = numpy.array([0,3,6])
        sort_chr_grp_spix = numpy.array([3,6,9])
        sort_chr_grp_len = numpy.array([3,3,3])

        # make sure we have a MarkerSet object
        self.assertIsInstance(m1, MarkerSet)

        # make sure it is sorted
        self.assertTrue(numpy.all(m1.chr_grp == sort_chr_grp))
        self.assertTrue(numpy.all(m1.chr_start == sort_chr_start))
        self.assertTrue(numpy.all(m1.chr_stop == sort_chr_stop))
        self.assertTrue(numpy.all(m1.mkr_name == sort_mkr_name))
        self.assertTrue(numpy.all(m1.chr_grp_stix == sort_chr_grp_stix))
        self.assertTrue(numpy.all(m1.chr_grp_spix == sort_chr_grp_spix))
        self.assertTrue(numpy.all(m1.chr_grp_len == sort_chr_grp_len))

    def test_len(self):
        # make fake data
        m1 = test_MarkerSet._create_fake_data()

        self.assertEqual(len(m1), 9)

    def _create_fake_data():
        # make some fake data
        chr_grp = numpy.string_(["1","1","2","2","2","3","3","1","3"])
        chr_start = numpy.array([1,3,6,5,4,9,7,2,8])
        chr_stop = numpy.array([1,3,6,5,4,9,7,2,8])
        mkr_name = numpy.string_(['a','c','d','e','f','g','h','b','i'])

        sort_chr_grp = numpy.string_(["1","1","1","2","2","2","3","3","3"])
        sort_chr_start = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_chr_stop = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_mkr_name = numpy.string_(['a','b','c','f','e','d','h','i','g'])
        sort_chr_grp_stix = numpy.array([0,3,6])
        sort_chr_grp_spix = numpy.array([3,6,9])
        sort_chr_grp_len = numpy.array([3,3,3])

        # construct object
        m1 = MarkerSet(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            auto_sort = True,
            auto_mkr_rename = False
        )

        return m1
