#!/usr/bin/env python

from greengenes.pick_chimeras import get_overlap, determine_taxon_conflict
from cogent.util.unit_test import TestCase, main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class PickChimerasTests(TestCase):
    def setUp(self):
        pass

    def test_get_overlap(self):
        """get overlap on ids"""
        b3 = [(1,2,3,4),(5,6,7,8),(9,10,11,12)]
        cs = [(5,6,7),(2,3,4)]
        exp = set([5])
        obs = get_overlap(b3, cs)
        self.assertEqual(obs,exp)

    def test_determine_taxon_conflict(self):
        """determine if there is a taxon conflict"""
        taxmap = {1:['a','b','c','d','e','f','g'],
                  2:['a','b','c','d','e','f','x'],
                  3:['a','b','z','q','w','t','y']}
        self.assertFalse(determine_taxon_conflict(taxmap, 1,2))
        self.assertTrue(determine_taxon_conflict(taxmap, 1,3))

if __name__ == '__main__':
    main()
