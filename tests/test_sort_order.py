#!/usr/bin/env python

from greengenes.sort_order import parse_dec, parse_lengths, merge_dec_length, sort_order
from cogent.util.unit_test import TestCase,main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class SortOrderTests(TestCase):
    def setUp(self):
                       #acc:(len,decision,407k,no sp.)
        self.map407k = {'a':(123,'named_isolate',True,True),
                        'b':(122,'named_isolate',True,True),
                        'c':(123,'named_isolate',True,False),
                        'd':(122,'named_isolate',True,False),
                        'e':(123,'clone',True,False),
                        'f':(122,'clone',True,False),
                        'g':(123,'named_isolate',False,True),
                        'h':(122,'named_isolate',False,True),
                        'i':(123,'named_isolate',False,False),
                        'j':(122,'named_isolate',False,False),
                        'k':(123,'clone',False,False),
                        'l':(122,'clone',False,False)}

    def test_parse_dec(self):
        """parse gg dec"""
        lines = ["a\tnamed_isolate\tblah sp. foo","b\tclone","c\tcrap",
                 "d\tnamed_isolate\timgood"]
        exp = {'a':('named_isolate', False),'b':('clone', False),
                'c':('clone', False),'d':('named_isolate', True)}
        obs = parse_dec(lines)
        self.assertEqual(obs,exp)

    def test_parse_lengths(self):
        lines = ["a\t123","b\t250"]
        exp = {'a':123,'b':250}
        obs = parse_lengths(lines)
        self.assertEqual(obs,exp)

    def test_merge_dec_length(self):
        in407k = set(['a','b'])
        lengths = {'a':123,'b':250, 'd':10, 'c':100}
        decinfo = {'a':('named_isolate', False),'b':('clone', False),
                'c':('clone', False),'d':('named_isolate', True)}
        exp = {'a':(123,'named_isolate',True,False),
               'b':(250,'clone',True,False),
               'c':(100,'clone',False,False),
               'd':(10,'named_isolate',False,True)}
        obs = merge_dec_length(in407k, lengths, decinfo)
        self.assertEqual(obs,exp)

    def test_sort_order(self):
        exp = ['a','b','c','d','e','f','g','h','i','j','k','l']
        obs = sort_order(self.map407k)
        self.assertEqual(obs,exp)

if __name__ == '__main__':
    main()
