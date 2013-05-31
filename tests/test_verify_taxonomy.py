#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from greengenes.verify_taxonomy import check_parse, check_n_levels, check_gap, \
        check_prefixes, ParseError

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class VerifyTaxonomy(TestCase):
    def setUp(self):
        pass

    def test_check_parse(self):
        """returns valid parsed or raises"""
        exp = ("1",["k__a","p__b","c__c","o__d","f__e","g__f","s__g"])
        obs = check_parse(good_string)
        self.assertEqual(obs, exp)

        self.assertRaises(ParseError, check_parse, bad_string)

    def test_check_n_levels(self):
        """requires N levels, or unclassified"""
        id_, parsed = check_parse(good_string)
        self.assertTrue(check_n_levels(parsed, 7))

        self.assertFalse(check_n_levels(parsed, 8))

        id_, parsed = check_parse(good_unclassified)
        self.assertTrue(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified1)
        self.assertFalse(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified2)
        self.assertFalse(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified3)
        self.assertFalse(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified4)
        self.assertFalse(check_n_levels(parsed, 7))

    def test_check_gap(self):
        """check if a gap exists in a string"""
        id_, parsed = check_parse(good_string)
        self.assertTrue(check_gap(parsed))

        id_, parsed = check_parse(good_trailing)
        self.assertTrue(check_gap, parsed)

        id_, parsed = check_parse(gap)
        self.assertFalse(check_gap(parsed))
    
    def test_check_prefixes(self):
        """Verify the expected prefixes are present"""
        prefixes = ['k','p','c','o','f','g','s']
        id_, parsed = check_parse(good_string)
        self.assertTrue(check_prefixes(parsed, prefixes))

        id_, parsed = check_parse(bad_prefix)
        self.assertFalse(check_prefixes(parsed, prefixes))
    
        self.fail("check duplicate ranks")

good_string = "1	k__a; p__b; c__c; o__d; f__e; g__f; s__g"
bad_string = "2 k__a; p__b; c__c; o__d; f__e; g__f; s__g" # space instead of tab
#multiple_parents = """10	k__a; p__b; c__c; o__d; f__e; g__f; s__g
#11	k__a; p__b; c__c; o__X; f__e; g__f; s__g"""
gap = "20	k__a; p__b; c__; o__d; f__e; g__f; s__g"
good_trailing = "30	k__a; p__b; c__c; o__d; f__e; g__; s__"
good_unclassified = "30	k__; p__; c__; o__; f__; g__; s__"
bad_unclassified1 = "40	Unclassified"
bad_unclassified2 = "50	k__x; unclassified"
bad_unclassified3 = "60	k__x; "
bad_unclassified4 = "70	k__x"
bad_nlevels = "80	k__a; p__b; c__c"
bad_prefix = "1	k__a; p__b; c__c; q__d; f__e; g__f; s__g"

if __name__ == '__main__':
    main()
