#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser
from cogent.seqsim.tree import RangeNode
from greengenes.util import GreengenesRecord, prune_tree, make_tree_arb_safe

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.9-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class UtilTests(TestCase):
    def setUp(self):
        pass

    def test_light_getsubtree(self):
        """get a subtree"""
        t = DndParser("(((a,b)c,(d,e)f)g,(h,i)j)k;", constructor=RangeNode)
        exp = "(((a,b)c)g,(h,i)j)k;"
        obs = prune_tree(t, ['a','b','h','i']).getNewick()
        self.assertEqual(obs,exp)

        t = DndParser("((((a,b)c,(d,e)f)g,(h,i)j,(k,l,m)n)o,(p,q)r)s;",
                        constructor=RangeNode)
        exp = "((((a,b)c,(d,e)f)g,(h,i)j,(k,l)n)o)s;"
        obs = prune_tree(t, ['a','b','d','e','h','i','k','l'])
        self.assertEqual(obs.getNewick(),exp)

        t = DndParser("((((a,b)c,(d,e)f)g,(h,i)j,(k,l,m)n)o,(p,q)r)s;",
                        constructor=RangeNode)
        exp = "((((a,b)c,(d,e)f)g,(h,i)j,(k)n)o)s;"
        obs = prune_tree(t, ['a','b','d','e','h','i','k'])
        self.assertEqual(obs.getNewick(),exp)

    def test_make_tree_arb_safe(self):
        """ARB can't handle single descendent nodes"""
        t = DndParser("(((((a)b)c)d)e)f;", constructor=RangeNode)
        exp = "(((((a,X)b,X)c,X)d,X)e,X)f;"
        make_tree_arb_safe(t)
        self.assertEqual(t.getNewick(), exp)

class GreengenesRecordTests(TestCase):
    def setUp(self):
        self.ggrecord = GreengenesRecord({'prokMSA_id':123})

    def test_init(self):
        """test initialization"""
        exp = {'prokMSA_id':123,
                'ncbi_acc_w_ver':None,
                'ncbi_gi':None,
                'db_name':None,
                'gold_id':None,
                'decision':None,
                'prokMSAname':None,
                'isolation_source':None,
                'clone':None,
                'organism':None,
                'strain':None,
                'specific_host':None,
                'authors':None,
                'title':None,
                'pubmed':None,
                'journal':None,                 
                'study_id':None,   
                'submit_date':None, 
                'country':None,
                'ncbi_tax_string':None,   
                'Silva_tax_string':None,
                'RDP_tax_string':None,
                'greengenes_tax_string':None,
                'non_ACGT_percent':None,
                'perc_ident_to_invariant_core':None,
                'small_gap_intrusions':None,
                'bellerophon':None,
                'bel3_div_ratio':None,
                'chim_slyr_a':None,
                'chim_slyr_b':None,
                'chim_slyr_a_tax':None,
                'chim_slyr_b_tax':None,
                'aligned_seq':None,
                'unaligned_seq':None
                } 
        obs = self.ggrecord
        self.assertEqual(obs,exp)
    
    def test_setTypes(self):
        """Sets types GG fields"""
        self.ggrecord['ncbi_acc_w_ver'] = 'asd'
        self.ggrecord.setTypes()
        self.assertEqual(self.ggrecord['prokMSA_id'], 123)
        self.assertEqual(self.ggrecord['ncbi_acc_w_ver'], 'asd')

    def test_getARBRules(self):
        """pull out arb rules right"""
        obs = sorted(self.ggrecord.getARBRules().split("\n\n"))
        exp = sorted(arbrules.split("\n\n")[:-1])
        for a,b in zip(obs,exp):
            self.assertEqual(a,b)
    
    def test_toGreengenesFormat(self):
        """Stringamify self"""
        obs = sorted(self.ggrecord.toGreengenesFormat().splitlines())
        exp = sorted(exp_testrecord.splitlines())

        self.assertEqual(obs,exp)
    
    def test_sanityCheck(self):
        """verify types are right"""
        self.assertEqual(self.ggrecord.sanityCheck(), None)
        self.ggrecord['prokMSA_id'] = "bad"
        self.assertRaises(ValueError, self.ggrecord.sanityCheck)
exp_testrecord = """BEGIN
prokMSA_id=123
ncbi_acc_w_ver=
ncbi_gi=
db_name=
gold_id=
decision=
prokMSAname=
isolation_source=
clone=
organism=
strain=
specific_host=
authors=
title=
pubmed=
journal=
study_id=
submit_date=
country=
ncbi_tax_string=
Silva_tax_string=
RDP_tax_string=
greengenes_tax_string=
non_ACGT_percent=
perc_ident_to_invariant_core=
small_gap_intrusions=
bellerophon=
bel3_div_ratio=
chim_slyr_a=
chim_slyr_b=
chim_slyr_a_tax=
chim_slyr_b_tax=
aligned_seq=
unaligned_seq=
END

"""
arbrules = """MATCH   "prokMSA_id\=*"
	SRT "*\=="
	WRITE "name"

MATCH   "ncbi_acc_w_ver\=*"
	SRT "*\=="
	WRITE "acc"

MATCH   "ncbi_gi\=*"
	SRT "*\=="
	WRITE "ncbi_gi"

MATCH   "db_name\=*"
	SRT "*\=="
	WRITE "IMG_taxon_OID"

MATCH   "gold_id\=*"
	SRT "*\=="
	WRITE "GOLD_stamp"

MATCH   "decision\=*"
	SRT "*\=="
	WRITE "sequence_type"

MATCH   "prokMSAname\=*"
	SRT "*\=="
	WRITE "full_name"

MATCH   "isolation_source\=*"
	SRT "*\=="
	WRITE "isolation_source"

MATCH   "clone\=*"
	SRT "*\=="
	WRITE "clone"

MATCH   "organism\=*"
	SRT "*\=="
	WRITE "organism"

MATCH   "strain\=*"
	SRT "*\=="
	WRITE "strain"

MATCH   "specific_host\=*"
	SRT "*\=="
	WRITE "specific_host"

MATCH   "authors\=*"
	SRT "*\=="
	WRITE "author"

MATCH   "title\=*"
	SRT "*\=="
	WRITE "title"

MATCH   "pubmed\=*"
	SRT "*\=="
	WRITE "pubmed"

MATCH   "journal\=*"
	SRT "*\=="
	WRITE "journal"

MATCH   "study_id\=*"
	SRT "*\=="
	WRITE "study_id"

MATCH   "submit_date\=*"
	SRT "*\=="
	WRITE "submit_date"

MATCH   "country\=*"
	SRT "*\=="
	WRITE "country"

MATCH   "ncbi_tax_string\=*"
	SRT "*\=="
	WRITE "ncbi_tax"

MATCH   "Silva_tax_string\=*"
	SRT "*\=="
	WRITE "Silva_tax"

MATCH   "RDP_tax_string\=*"
	SRT "*\=="
	WRITE "RDP_tax"

MATCH   "greengenes_tax_string\=*"
	SRT "*\=="
	WRITE "greengenes_tax"

MATCH   "non_ACGT_percent\=*"
	SRT "*\=="
	WRITE "percent_non_ACGT"

MATCH   "perc_ident_to_invariant_core\=*"
	SRT "*\=="
	WRITE "percent_invariant_match"

MATCH   "small_gap_intrusions\=*"
	SRT "*\=="
	WRITE "small_gap_intrusions"

MATCH   "bellerophon\=*"
	SRT "*\=="
	WRITE "bellerophon"

MATCH   "bel3_div_ratio\=*"
	SRT "*\=="
	WRITE "bel_div_ratio"

MATCH   "chim_slyr_a\=*"
	SRT "*\=="
	WRITE "chim_slyr_a"

MATCH   "chim_slyr_b\=*"
	SRT "*\=="
	WRITE "chim_slyr_b"

MATCH   "chim_slyr_a_tax\=*"
	SRT "*\=="
	WRITE "chim_slyr_a_tax"

MATCH   "chim_slyr_b_tax\=*"
	SRT "*\=="
	WRITE "chim_slyr_b_tax"

MATCH   "aligned_seq\=*"
	SRT "*\=="
	WRITE "aligned_seq"

MATCH   "unaligned_seq\=*"
	SRT "*\=="
	WRITE "unaligned_seq"
"""

if __name__ == '__main__':
    main()

