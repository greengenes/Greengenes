#!/usr/bin/env python

"""chop the genus name from the species names"""

from greengenes.chop_taxonomy import chop_name
from cogent.util.unit_test import TestCase, main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class ChopTaxonomyTests(TestCase):
    def setUp(self):
        pass

    def test_chop_name(self):
        """Take a species name and choppermafy it"""
        str1 = "s__foo bar"
        exp1 = "s__bar"
        str2 = "s__bar"
        exp2 = "s__bar"
        str3 = "s__foo bar taco"
        exp3 = "s__bar"
        str4 = "s__Candidatus bar foo"
        exp4 = "s__foo"
        str5 = "s__candidatus bar foo"
        exp5 = "s__foo"
        str6 = "s__candidatus sp. foo"
        exp6 = "s__"
        str7 = "s__candidatus foo"
        exp7 = "s__"

        obs1 = chop_name(str1)
        obs2 = chop_name(str2)
        obs3 = chop_name(str3)
        obs4 = chop_name(str4)
        obs5 = chop_name(str5)
        obs6 = chop_name(str6)
        obs7 = chop_name(str7)

        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)
        self.assertEqual(obs4, exp4)
        self.assertEqual(obs5, exp5)
        self.assertEqual(obs6, exp6)
        self.assertEqual(obs7, exp7)

if __name__ == '__main__':
    main()
