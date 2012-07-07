#!/usr/bin/env python

from cogent.util.unit_test import TestCase,main
from greengenes.write import write_sequence, write_gb_summary, \
    write_obs_record
from StringIO import StringIO

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class WriteTests(TestCase):
    def setup(self):
        pass

    def test_write_sequence(self):
        """writes a sequence"""
        f = StringIO()
        exp = ">abc\nAATTGG\n>efg\nATGGG\n"
        write_sequence(f, "abc", "AATTGG")
        write_sequence(f, "efg", "ATGGG")
        f.seek(0)
        obs = f.read()
        self.assertEqual(obs,exp)

    def test_write_gb_summary(self):
        """writes a gb summary line"""
        f = StringIO()
        # its just an honest method...
        exp = "a\tb\tc\tg\nx\ty\tz\tfoo\n"
        write_gb_summary(f, 'a', ['b','c','g'])
        write_gb_summary(f, 'x', ['y','z','foo'])
        f.seek(0)
        obs = f.read()
        self.assertEqual(obs, exp)

    def test_write_obs_record(self):
        """writes observed records"""
        f = StringIO()
        exp = "a\nb\nc\n"
        write_obs_record(f, 'a') 
        write_obs_record(f, 'b') 
        write_obs_record(f, 'c')
        f.seek(0)
        obs = f.read()
        self.assertEqual(obs,exp)

if __name__ == '__main__':
    main()
