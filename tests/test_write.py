#!/usr/bin/env python

from cogent.util.unit_test import TestCase,main
from greengenes.write import write_sequence,  \
    write_obs_record, write_gg_record
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

    #def test_write_gb_summary(self):
    #    """writes a gb summary line"""
    #    f = StringIO()
    #    # its just an honest method...
    #    exp = "a\tb\tc\tg\nx\ty\tz\tfoo\n"
    #    write_gb_summary(f, 'a', ['b','c','g'])
    #    write_gb_summary(f, 'x', ['y','z','foo'])
    #    f.seek(0)
    #    obs = f.read()
    #    self.assertEqual(obs, exp)

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

    def test_write_gg_record(self):
        """Writes a gg record"""
        exp = sorted(['BEGIN',
            'prokmsa_id=123',
            'gg_id=',
            'hugenholtz_tax_string=',
            'ncbi_acc_w_ver=xyz',
            'ncbi_gi=333',
            'n_pos_aligned=',
            'n_pos_unaligned=',
            'db_name=',
            'gold_id=',
            'decision=',
            'prokmsaname=',
            'isolation_source=',
            'clone=foo',
            'organism=',
            'strain=',
            'specific_host=',
            'authors=',
            'title=',
            'pubmed=123',
            'journal=',
            'study_id=',
            'submit_date=',
            'country=',
            'ncbi_tax_string=',
            'silva_tax_string=',
            'rdp_tax_string=',
            'greengenes_tax_string=',
            'non_acgt_percent=0.5',
            'perc_ident_to_invariant_core=',
            'small_gap_intrusions=',
            'bellerophon=',
            'bel3_div_ratio=',
            'chim_slyr_a=',
            'chim_slyr_b=',
            'chim_slyr_a_tax=',
            'chim_slyr_b_tax=',
            'aligned_seq=',
            'unaligned_seq=',
            'END',''])
        ggrec = GreengenesRecord({'prokmsa_id':123,'ncbi_acc_w_ver':'xyz',
                                'ncbi_gi':'333','pubmed':123,'clone':'foo',
                                  'non_acgt_percent':'0.5'})
        f = StringIO()
        write_gg_record(f, ggrec)
        f.seek(0)
        obs = sorted(f.read().splitlines())
        self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
