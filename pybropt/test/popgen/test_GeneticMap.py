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

class test_GeneticMap(unittest.TestCase):
    """docstring for test_GeneticMap."""

    def test_constructor(self):
        m1 = test_GeneticMap._create_fake_data()

        sort_chr_grp = numpy.string_(["1","1","1","2","2","2","3","3","3"])
        sort_chr_start = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_chr_stop = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_map_pos = numpy.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
        sort_mkr_name = numpy.string_(['a','b','c','f','e','d','h','i','g'])
        sort_map_fncode = numpy.string_(['h','h','h','h','h','h','h','h','h'])
        sort_chr_grp_stix = numpy.array([0,3,6])
        sort_chr_grp_spix = numpy.array([3,6,9])
        sort_chr_grp_len = numpy.array([3,3,3])

        self.assertIsInstance(m1, popgen.GeneticMap)

        # make sure it is sorted
        self.assertTrue(numpy.all(m1.chr_grp == sort_chr_grp))
        self.assertTrue(numpy.all(m1.chr_start == sort_chr_start))
        self.assertTrue(numpy.all(m1.chr_stop == sort_chr_stop))
        self.assertTrue(numpy.all(m1.map_pos == sort_map_pos))
        self.assertTrue(numpy.all(m1.mkr_name == sort_mkr_name))
        self.assertTrue(numpy.all(m1.map_fncode == sort_map_fncode))
        self.assertTrue(numpy.all(m1.chr_grp_stix == sort_chr_grp_stix))
        self.assertTrue(numpy.all(m1.chr_grp_spix == sort_chr_grp_spix))
        self.assertTrue(numpy.all(m1.chr_grp_len == sort_chr_grp_len))

    def test_from_egmap(self):
        # get test data file path
        data_path = test_dir + "/maize_genetic_map_McMullen_2009_US_NAM.egmap"

        # load data
        gmap = popgen.GeneticMap.from_egmap(data_path, "kosambi")

        self.assertIsInstance(gmap, popgen.GeneticMap)

    def test_to_egmap(self):
        # make file paths
        data_path = test_dir + "/maize_genetic_map_McMullen_2009_US_NAM.egmap"
        export_path = test_dir + "/test_to_egmap_export.egmap"

        # load data
        gmap = popgen.GeneticMap.from_egmap(data_path, "kosambi")

        # export data
        gmap.to_egmap(export_path)

        self.assertTrue(os.path.exists(export_path))

    def _create_fake_data():
        # make some fake data
        chr_grp = numpy.string_(["1","1","2","2","2","3","3","1","3"])
        chr_start = numpy.array([1,3,6,5,4,9,7,2,8])
        chr_stop = numpy.array([1,3,6,5,4,9,7,2,8])
        map_pos = numpy.array([0.1,0.3,0.6,0.5,0.4,0.9,0.7,0.2,0.8])
        mkr_name = numpy.string_(['a','c','d','e','f','g','h','b','i'])
        map_fncode = numpy.string_(['h','h','h','h','h','h','h','h','h'])
        mapfn = 'haldane'

        sort_chr_grp = numpy.string_(["1","1","1","2","2","2","3","3","3"])
        sort_chr_start = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_chr_stop = numpy.array([1,2,3,4,5,6,7,8,9])
        sort_mkr_name = numpy.string_(['a','b','c','f','e','d','h','i','g'])
        sort_chr_grp_stix = numpy.array([0,3,6])
        sort_chr_grp_spix = numpy.array([3,6,9])
        sort_chr_grp_len = numpy.array([3,3,3])

        # construct object
        m1 = popgen.GeneticMap(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            map_pos = map_pos,
            mkr_name = mkr_name,
            map_fncode = map_fncode,
            auto_sort = True,
            auto_mkr_rename = False,
            auto_fncode = False,
            auto_spline = True
        )

        return m1
