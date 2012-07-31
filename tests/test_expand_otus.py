#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from greengenes.expand_otus import expand

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class ExpandTests(TestCase):
    def setUp(self):
        pass

    def test_expand(self):
        """expands from ref seqs"""
        tax = {'a':"string does not matter",
               'b':'asdasd',
               'X123':'w000000p'}
        otus = [('otu_id',['a','other','seq']),
                ('asdasd',['b']),
                ('a123123',['X123','1','2','3'])]
        exp = [('a','string does not matter'),
               ('other','string does not matter'),
               ('seq','string does not matter'),
               ('b','asdasd'),
               ('X123','w000000p'),
               ('1','w000000p'),
               ('2','w000000p'),
               ('3','w000000p')]
        obs = expand(tax,otus)
        self.assertEqual(obs,exp)

if __name__ == '__main__':
    main()
