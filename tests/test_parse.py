#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from greengenes.parse import parse_column, parse_gg_summary_flat, \
        parse_invariants, parse_otus, parse_b3_chimeras, \
        parse_cs_chimeras, parse_uchime_chimeras
from greengenes.util import GreengenesRecord
from StringIO import StringIO

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"
    
class ParseTests(TestCase):
    def setUp(self):
        pass

    def test_parse_b3_chimeras(self):
        """Parse B3 chimera output"""
        inlines = ["a\t1.2\tparentfoo\tparentbar\n",
                   "123x\t1.323\txyz\tabc\n"]
        exp = [('a',1.2,'parentfoo','parentbar'),
               ('123x',1.323,'xyz','abc')]
        obs = parse_b3_chimeras(inlines)
        self.assertEqual(obs,exp)

    def test_parse_uchime_chimeras(self):
        """Parse uchime chimera output"""

        inlines = ["1.7308	GU234509.1	94_12740_108965	94_6884_155789	99.6	98.0	95.5	93.9	98.0	51	1	3	22	2	2	1.63	Y\n",
                   "0.0000	HQ316978.1	*	*	*	*	*	*	*	*	*	*	*	*	*	*	N\n",
                   "0.0190	HQ316968.1	94_6884_155789	94_28256_783727	97.8	97.5	97.4	99.3	97.5	6	0	18	4	3	8	0.33	N\n"]
        exp = [('GU234509.1', 1.7308, '94_12740_108965','94_6884_155789')]
        obs = parse_uchime_chimeras(inlines)
        self.assertEqual(obs, exp)

    def test_parse_cs_chimeras(self):
        """Parse ChimeraSlayer output"""
        inlines = ["a\tparentfoo\tparentbar\n",
                   "123x\txyz\tabc\n"]
        exp = [('a','parentfoo','parentbar'),
               ('123x','xyz','abc')]
        obs = parse_cs_chimeras(inlines)
        self.assertEqual(obs,exp)

    def test_parse_otus(self):
        """parse me some tasty otus"""
        inlines = ["1\tx\ty\tz\n", "10\tpoo\n", "20xad\t30\t20\tfoo\n"]
        exp = [("1",["x","y","z"]), ("10",["poo"]), 
               ("20xad", ["30", "20", "foo"])]
        obs = parse_otus(inlines)
        self.assertEqual(obs,exp)

    def test_parse_invariants(self):
        """Parse the invariants file"""
        exp = [('inv_a','NNNNAANNANNTTNGNANNNAAANN',10),
               ('inv_b','AANNANNANANANANNTTATTANNN',13)]
        obs = parse_invariants(StringIO(invariants))    
        self.assertEqual(obs,exp)

    def test_parse_column(self):
        """Parse existing records"""
        recs = StringIO(existing_records)
        exp = set(['abc','def','xyz','foo','bar'])
        obs = parse_column(recs)
        self.assertEqual(obs,exp)

    def test_parse_gg_summary_flat(self):
        """Parse the gg summary files from flat_files.py"""
        exp = [GreengenesRecord({'prokMSA_id':'1', 'ncbi_acc_w_ver':'xyzf'}),
               GreengenesRecord({'prokMSA_id':'25', 'ncbi_acc_w_ver':'abcd',
                                 'country':'australia'}),
               GreengenesRecord({'prokMSA_id':'50', 'ncbi_acc_w_ver':'223xx'})]
        obs = list(parse_gg_summary_flat(StringIO(gg_summary)))

        self.assertEqual(obs,exp)

invariants = """>inv_a
NNNNAANNANNTTNGNANNNAAANN
>inv_b
AANNANNANANANANNTTATTANNN
"""

existing_records = """#accession\tprokMSA_id
abc\t10
def\t20
xyz
foo\t
bar\t30
"""

gg_summary = """#prokMSA_id\tncbi_acc_w_ver\tcountry
1\txyzf\t
25\tabcd\taustralia
50\t223xx\t
"""

if __name__ == '__main__':
    main()
