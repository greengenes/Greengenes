#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from greengenes.parse import parse_existing_records
from StringIO import StringIO

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class ParseTests(TestCase):
    def setUp(self):
        pass

    def test_parse_existing_records(self):
        """Parse existing records"""
        recs = StringIO(existing_records)
        exp = set(['abc','def','xyz','foo','bar'])
        obs = parse_existing_records(recs)
        self.assertEqual(obs,exp)

existing_records = """#accession\tprokMSA_id
abc\t10
def\t20
xyz
foo\t
bar\t30
"""

if __name__ == '__main__':
    main()
