#!/usr/bin/env python

from greengenes.db import GreengenesMySQL
from cogent.util.unit_test import TestCase,main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class GGMySQLTests(TestCase):
    def setUp(self):
        self.db = GreengenesMySQL(debug=True)

    def tearDown(self):
        del self.db

    def test_get_max_seqid(self):
        """get the max seq id"""
        # as of gg_13_5... 
        self.assertEqual(self.db._get_max_seqid(), 2796094)

    def test_locking(self):
        """Lock some tables"""
        self.assertRaises(ValueError, self.db._lock, [['a','b','c']])
        self.assertEqual(self.db._have_locks, False)
        
        self.db._lock([['greengenes','g','read']])
        self.assertTrue(self.db._have_locks)
        self.assertEqual(self.db._lock_aliases, {'greengenes':'g'})
        
        self.db._unlock()
        self.assertFalse(self.db._have_locks)
        self.assertEqual(self.db._lock_aliases, {})

    def test_updatePyNASTSequences(self):
        """Update the seqs"""
        self.db.cursor.execute('select gg_id from greengenes where pynast_aligned_seq_id is null')
        ggids = [str(i[0]) for i in self.db.cursor.fetchall()]
        self.assertGreaterThan(len(ggids), 0)

        seqs = dict([(i,'ATGC_%s' % i) for i in ggids])
        self.db.updatePyNASTSequences(seqs)

        exp = seqs
        obs_cur = self.db.cursor.execute("""select g.gg_id, s.sequence 
                                from greengenes g inner join sequence s
                                    on g.pynast_aligned_seq_id=s.seq_id
                                where g.gg_id in (%s)""" % ','.join(exp.keys()))
        obs = dict([(str(id_), seq) for id_, seq in self.db.cursor.fetchall()])
        self.assertEqual(obs,exp)

    def test_updateSSUAlignSequences(self):
        """Update the seqs"""
        self.db.cursor.execute('select gg_id from greengenes where aligned_seq_id is null')
        ggids = [str(i[0]) for i in self.db.cursor.fetchall()]
        self.assertGreaterThan(len(ggids), 0)

        seqs = dict([(i,'ATGC_%s' % i) for i in ggids])
        self.db.updateSSUAlignSequences(seqs)

        exp = seqs
        obs_cur = self.db.cursor.execute("""select g.gg_id, s.sequence 
                                from greengenes g inner join sequence s
                                    on g.aligned_seq_id=s.seq_id
                                where g.gg_id in (%s)""" % ','.join(exp.keys()))
        obs = dict([(str(id_), seq) for id_, seq in self.db.cursor.fetchall()])
        self.assertEqual(obs,exp)

    def test_updateUnlignedSequences(self):
        """Update the seqs"""
        self.db.cursor.execute('select gg_id from greengenes where unaligned_seq_id is null')
        ggids = [str(i[0]) for i in self.db.cursor.fetchall()]
        self.assertGreaterThan(len(ggids), 0)

        seqs = dict([(i,'ATGC_%s' % i) for i in ggids])
        self.db.updateUnalignedSequences(seqs)

        exp = seqs
        obs_cur = self.db.cursor.execute("""select g.gg_id, s.sequence 
                                from greengenes g inner join sequence s
                                    on g.unaligned_seq_id=s.seq_id
                                where g.gg_id in (%s)""" % ','.join(exp.keys()))
        obs = dict([(str(id_), seq) for id_, seq in self.db.cursor.fetchall()])
        self.assertEqual(obs,exp)

    def test_update_seq_field(self):
        """tested implicitly by other sequence field update methods"""
        pass

if __name__ == '__main__':
    main()

