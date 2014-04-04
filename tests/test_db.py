#!/usr/bin/env python

from greengenes.db import GreengenesDB
from unittest import TestCase,main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class GGDBTests(TestCase):
    def setUp(self):
        self.db = GreengenesDB(debug=True)

    def tearDown(self):
        del self.db

    def test_update_greengenes_tax(self):
        """Update Greengenes taxonomy strings"""
        taxmap = {86:"foo bar", 3:'yay number 3', '4':'slightly malformed id',
                  57:None}
        exp = dict([(int(k), v) for k,v in taxmap.items()])
        self.db.update_greengenes_tax(taxmap, "TESTING")
        obs = self.db.get_greengenes_tax_multiple([86,3,4,57])
        self.assertEqual(obs,exp)

    def test_update_ncbi_tax(self):
        """Update NCBI taxonomy strings"""
        taxmap = {86:"foo bar", 3:'yay number 3', '4':'slightly malformed id',
                  57:None}
        exp = dict([(int(k), v) for k,v in taxmap.items()])
        self.db.update_ncbi_tax(taxmap, "TESTING")
        obs = self.db.get_ncbi_tax_multiple([86,3,4,57])
        self.assertEqual(obs,exp)

    def test_update_tax(self):
        """implicitly tested by other taxonomy update methods"""
        pass

    def test_get_greengenes_tax(self):
        """Get a Greengenes taxonomy string by GGID"""
        exp_86 = "k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobrevibacter; s__"
        exp_3 = None
        exp_4 = "k__Archaea; p__[Parvarchaeota]; c__[Parvarchaea]; o__YLA114; f__; g__; s__"

        obs_86 = self.db.get_greengenes_tax(86)
        obs_3 = self.db.get_greengenes_tax(3)
        obs_4 = self.db.get_greengenes_tax(4)

        self.assertEqual(obs_86, exp_86)
        self.assertEqual(obs_3, exp_3)
        self.assertEqual(obs_4, exp_4)

    def test_get_greengenes_tax_multiple(self):
        """Get multiple taxonomy strings by GGIDs"""
        exp = {3:None,4:"k__Archaea; p__[Parvarchaeota]; c__[Parvarchaea]; o__YLA114; f__; g__; s__",
               86:"k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobrevibacter; s__"}
        obs = self.db.get_greengenes_tax_multiple([3,4,86])
        self.assertEqual(obs, exp)

    def test_get_ncbi_tax(self):
        """Get a NCBI taxonomy string by GGID"""
        exp_3 = "Archaea; Euryarchaeota; Methanococci; Methanococcales; Methanocaldococcaceae; Methanocaldococcus"
        exp_4 = "Archaea; environmental samples"
        exp_86 = "Archaea; Euryarchaeota; Methanobacteria; Methanobacteriales; Methanobacteriaceae; Methanobrevibacter"

        obs_3 = self.db.get_ncbi_tax(3)
        obs_4 = self.db.get_ncbi_tax(4)
        obs_86 = self.db.get_ncbi_tax(86)

        self.assertEqual(obs_3, exp_3)
        self.assertEqual(obs_4, exp_4)
        self.assertEqual(obs_86, exp_86)

    def test_get_ncbi_tax_multiple(self):
        """Get multiple NCBI taxonomy strings by GGID"""
        exp = {3:"Archaea; Euryarchaeota; Methanococci; Methanococcales; Methanocaldococcaceae; Methanocaldococcus",
               4:"Archaea; environmental samples",
               86:"Archaea; Euryarchaeota; Methanobacteria; Methanobacteriales; Methanobacteriaceae; Methanobrevibacter"}
        obs = self.db.get_ncbi_tax_multiple([3,4,86])
        self.assertEqual(obs, exp)

    def test_get_multiple_taxonomy_string_ggid(self):
        """implicitly tested by other taxonomy getters"""
        pass

    def test_get_single_taxonomy_string_ggid(self):
        """implicitly tested by other taxonomy getters"""
        pass

    def test_to_arb(self):
        """Bulk fetch records for ARB import"""
        ids = [86,79,50]
        exp = [e for e in example_arb_recs.splitlines() if e]
        obs = self.db.to_arb(ids, "aligned_seq_id")
        obs = [o.strip() for o in obs if o.strip()]

        self.assertEqual(sorted(obs), sorted(exp))

    def test_get_pynast_seq(self):
        """Get a single PyNAST sequence by GG_ID"""
        exp = gg_id_86_pynast_seq
        obs = self.db.get_pynast_seq(86)
        self.assertEqual(obs,exp)

    def test_get_ssualign_seq(self):
        """Get a single SSU Align sequence by GG_ID"""
        exp = gg_id_4_ssualigned_seq
        obs = self.db.get_ssualign_seq(4)
        self.assertEqual(obs, exp)

    def test_get_unaligned_seq(self):
        exp = gg_id_86_unaligned_seq
        obs = self.db.get_unaligned_seq(86)
        self.assertEqual(obs, exp)

    def test_insert_sequence(self):
        exp_seq = "AATTGGCC"
        exp_id = self.db._get_max_seqid() + 1
        obs_id = self.db.insert_sequence(exp_seq)
        self.db.cursor.execute("select sequence from sequence where seq_id=%d"\
                % exp_id)
        obs_seq = self.db.cursor.fetchone()[0]
        self.assertEqual(obs_id, exp_id)
        self.assertEqual(obs_seq, exp_seq)

    def test_insert_taxonomy(self):
        exp_tax = "a; b; c"
        exp_name = "TEST"
        exp_id = self.db._get_max_taxid() + 1
        obs_id = self.db.insert_taxonomy(exp_tax, exp_name)
        self.db.cursor.execute("""select tax_string, tax_version
                                  from taxonomy where tax_id=%d""" % exp_id)
        obs_tax, obs_name = self.db.cursor.fetchone()
        self.assertEqual(obs_id, exp_id)
        self.assertEqual(obs_tax, exp_tax)
        self.assertEqual(obs_name, exp_name)

    def test_insert_record(self):
        exp_id = self.db._get_max_ggid() + 1
        exp_rec = {'ncbi_acc_w_ver': 'test',
                   'decision':'test_dec'}
        exp_name = "test_name"
        obs_id = self.db.insert_record(exp_rec.copy(), exp_name)

        self.db.cursor.execute("""
                select ncbi_acc_w_ver, decision
                from record
                where gg_id=%d""" % exp_id)
        obs_rec = self.db.cursor.fetchone()
        obs_rec = {k:v for k,v in zip(["ncbi_acc_w_ver","decision"], obs_rec)}

        self.db.cursor.execute("select name from gg_release where gg_id=%d" \
                              % exp_id)
        obs_name = self.db.cursor.fetchone()[0]

        self.assertEqual(obs_id, exp_id)
        self.assertEqual(obs_rec, exp_rec)
        self.assertEqual(obs_name, exp_name)

    def test_insert_otu(self):
        self.db.insert_otu(49, [13,7,32], 'test', 0.123, '13_5')
        self.db.cursor.execute("select * from otu_cluster")
        exp = (1, 49, 2150456, 0.123, 'test')
        obs = self.db.cursor.fetchall()[0]
        self.assertEqual(obs, exp)

        self.db.cursor.execute("select * from otu")
        exp = [(1, 1, 13), (2, 1, 7), (3, 1, 32), (4, 1, 49)]
        obs = self.db.cursor.fetchall()
        self.assertEqual(sorted(obs), exp)

    def test_get_sequence_ggid(self):
        """Implicitly tested by other sequence obtaining methods"""
        pass

    def test_get_max_seqid(self):
        """get the max seq id"""
        # as of gg_13_5...
        self.assertEqual(self.db._get_max_seqid(), 2811319)

    def test_get_max_taxid(self):
        """get the max seq id"""
        # as of gg_13_5...
        self.assertEqual(self.db._get_max_taxid(), 7558683)

    def test_update_pynast_seq(self):
        """Update the seqs"""
        self.db.cursor.execute('select gg_id from record where pynast_aligned_seq_id is null')
        ggids = [str(i[0]) for i in self.db.cursor.fetchall()]
        self.assertNotEqual(len(ggids), 0)

        seqs = dict([(i,'ATGC_%s' % i) for i in ggids])
        self.db.update_pynast_seq(seqs)

        exp = seqs
        obs_cur = self.db.cursor.execute("""select g.gg_id, s.sequence
                                from record g inner join sequence s
                                    on g.pynast_aligned_seq_id=s.seq_id
                                where g.gg_id in (%s)""" % ','.join(exp.keys()))
        obs = dict([(str(id_), seq) for id_, seq in self.db.cursor.fetchall()])
        self.assertEqual(obs,exp)

    def test_update_ssualign_seq(self):
        """Update the seqs"""
        self.db.cursor.execute('select gg_id from record where aligned_seq_id is null')
        ggids = [str(i[0]) for i in self.db.cursor.fetchall()]
        self.assertNotEqual(len(ggids), 0)

        seqs = dict([(i,'ATGC_%s' % i) for i in ggids])
        self.db.update_ssualign_seq(seqs)

        exp = seqs
        obs_cur = self.db.cursor.execute("""select g.gg_id, s.sequence
                                from record g inner join sequence s
                                    on g.aligned_seq_id=s.seq_id
                                where g.gg_id in (%s)""" % ','.join(exp.keys()))
        obs = dict([(str(id_), seq) for id_, seq in self.db.cursor.fetchall()])
        self.assertEqual(obs,exp)

    def test_update_unaligned_seq(self):
        """Update the seqs"""
        self.db.cursor.execute('select gg_id from record where unaligned_seq_id is null')
        ggids = [str(i[0]) for i in self.db.cursor.fetchall()]
        self.assertNotEqual(len(ggids), 0)

        seqs = dict([(i,'ATGC_%s' % i) for i in ggids])
        self.db.update_unaligned_seq(seqs)

        exp = seqs
        obs_cur = self.db.cursor.execute("""select g.gg_id, s.sequence
                                from record g inner join sequence s
                                    on g.unaligned_seq_id=s.seq_id
                                where g.gg_id in (%s)""" % ','.join(exp.keys()))
        obs = dict([(str(id_), seq) for id_, seq in self.db.cursor.fetchall()])
        self.assertEqual(obs,exp)

    def test_update_seq_field(self):
        """tested implicitly by other sequence field update methods"""
        pass

gg_id_86_unaligned_seq = "GCTACTGCTATTGGGATTCGATTAAGCCATNCAAGTTGAACGAATTTAGATTCGTGGCGTACGGCTCAGTAACACGTGGATAACCTACCCTTAGGACTGGGATAACTCTGGGAAACTGGGGATAATACCGGATAGGCAATTTTTCCTGTAATGGTTTTTTGTTTAAATGTTTTTTTTCGCCTAAGGATGGGTCTGCGGCAGATTAGGTAGTTGGTTAGGTAATGGCTTACCAAGCCGTTGATCTGTACGGGTTGTGAGAGCAAGAGCCCGGAGATGGAACCTGAGACAAGGTTCCAGGCCCTACGGGGCGCAGCAGGCGCGAAACCTCCGCAATGTGAGAAATCGCGACGGGGGGATCCCAAGTGCCATTCTTAACGGGATGGCTTTTCATTAGTGTAAAGAGCTTTTGGAATAAGAGCTGGGCAAGACCGTTGCCAACCGCCGCGGTAACACCGTCAGCTCTAGTGGTAGCAGTTTTTATTGGGCCTAAAGCGTCCGTAGCCGGTTTATTAAGTCTCTGGTGAAATCCTGTAGCTTAACTGTGGGAATTGCTGGAGATACTAGTAGACTTGAGATCGGGAGAGGTTAGAGGTACTCCCAGGGTAGAGGTGAAATTCTGTAATCCTGGGAGGACCGCCTGTGGCGAAGGCGTCTAGCTGGAACGATTCTGCCGGTGAGGGACGAAAGCTAGGGGCGCGAACCGGATTAGATACCCGGGTAGTCCTAGCTGTAAACGATGCGGACTTGGTGTTGGGATGGCTTTGAGCTGCTCCAGTGCCGAAGGGAAGCTGTTAAGTCCGCCGCCTGGGAAGTACGGTCGCAAGACTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAACGCGTGGAGCCTGCGGTTTAATTGGATTCAACGCCGGACATCTCACCAGAGGCGACAGCTGTATGATGACCAGTTTGATGAGCTTGTTTGACTAGCTGAGAGGAGGTGCATGGCCGCCGTCAGCTCGTACCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCACGCCCTTAGTTACCAGCGGATTCTAGTGAATGCCGGGCACACTAGGGGGACCGCCTGTGATAAATAGGAGGAAGGAGTGGACGACGGTAGGTCCGTATGCCCCGAATCCTCTGGGCAACACGCGGGCTACAATGGATGAGACAATGGGTTCCGACGCCGAAAGGTGGAGGTAATCCTCTAAACTTATTCGTAGTTCGGATTGAGGACTGTAACTCGTTCTCATGAAGCTGGAATGCGTAGTAATCGCGTGTCATAATCGCGCGGTGAATACGTCCCTGCTCCTTGGACACACCGCCCGTCACG"
gg_id_86_pynast_seq = "-----------------------------------------------------------------------------------------------------------------------------------------GC-T-AC-T-GC--TAT-T--G-GG-ATTC--GA---T-T--AAGCCA-T-NC-A-AGT-TGA-A-CGA---------A-T------------------------------------------------------------------------------------------------TTA-GA---------------------------------------------------------------------------------------------------------------------TT--CG-T-GG-C-GT-A--C-------------GGC-TCAGT-A--AC-AC-G-T-G-GA---TAA--C-CT-A--C-C-CTT--AG-G------------------------------------------------------------------A-CT----GGG-AT-AA-CTC-------------------------T-G-G-----------------------GAA-A---CTG-GGG-ATAA-TA---CC-G--G-AT-A---------------------------------G--G-C-A--A--T-----------------TT-TTCC-T-----------------------------------------------------------------------------------------------------------------------G-TA-A--------------------------------------------------------------------------------------------------------------------------------------T-G-GT-T-T---------------T--T-T-G-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAATGTTTT---------------------------------------------------------------------------------------------------------------------------------------TTT-T----------------------------------------------------------------------------------------------------------------------------------C-----G--------------C----C-T---A-AG-G---AT---G-G-----G-TCT-GCG--G-CAG--A------TT--A--G-GT-A----G---TTGG-T-T-AG-G-T----AAT-GG-C-T-T-ACCA--A-GC-C-G--T-TG-A------------TCT-G-T------AC-GG-G-T-TGT-G-AG----A--GC-AA--G-AG-C-CC-GGAG-A-TGGAA--C-C-TG-A-GA-C-AA-G-G-TTCCAG-GCCC-TAC-G--G-G-G-C-GC-A-GC-A-G-GC---GC-G-A-AAC-CTCCG-C-AA-T-GT--GA-GA-A----A-T-CG-C-GA-CG-GG-GGGA-TCCC-A-AG-T---G-C-C--A--T----------T-C-T--------TA-AC-------------G-G-G--------A--T-GGC--------TT-TT-C-A--T-TAG----T------------------------------G--T--AA-A---G----A------------------------------G-C-TT-T-TG-G---------AA-----------TAAGA-GCTGGG-C--AA---G-AC-CGTT--GCCA--A-C---C--GCCG---C-GG--TA-AC--AC---CG-TC-AGC-TCT-A-G-TG-GTAG-C-AGT-TT-TT-A--T-T--GGGC-CTA----AA-GCGT-CC--G-TA-G-C-C-G------------G--T-TT-A-T-T-AA----G-T-C-T---C-TGG-TG-A-AA-TC--CT-GTA-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CT-G-T-GG-GA-AT---T-G-C-T-G-G--------A--GA-T-A-C-T-A-GTA--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-A-T-----C-GG--GA-G-A------------G-GT-T-AG-A----GG--TACT-CCC-A-GG--GT-A-GAG-GTGAAA-TT-CTG-TAAT-C-CT-G-GGA--GG-A-CC-G-CC-TG--T--G--GC-GAA-G--G-C---G----T--C-T-AGCTG------G-AA-C---------------------------------------------------------------GATT-C-T--GC--CG-----GT-GA-GG--G-A-CGA--AA-G-C--------------T-AGGG-GCG-C-G-AACC--GG-ATTA-G-ATA-C-----CC-G-G-GTA-G-T----C-CT--A-G-CTG-T-AAA--C-GATG-CG--GA-CT---------T-GG--T--G-T-TG-G-GA-T--G--GC----------------------------------------------------------------------------------TTT-GA---------------------------------------------------------------------------------------------------------------------------------------------GC---T-G-C-TC--C-A-G-T-GC-C------GA--A----GG-GAA--GC-T-G-T--T--AA-GT--C----C-GCC-GCC-T-G-GG-AAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTT-AAA---------GGAA-TTG-GCGGG-G-G-AGCA----CCA--C-A-A-CGC-GT-G--G--AG-CC-T--GC-GGT-TT-AATT-G-G-ATT-CAAC-G-CC-G-GA-C-A-TC-TC-A-CC-AGAGG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CGACAGC--TGTATGATGACCAGTTTGATGAGCTTGTT-TGA-CTAGCTGAG-A-G-G-A-GGTG-CA-TGG-CC--GCC-GTC-A-GC-TC---G-TA-CC-G--TGA-GG-CGT-C-CT-G-TT-AA-GT-CAGGC-AA--------C-GAG-CGA-G-ACC-C-A-CG--CC--C-TTAG--T-T-A-C-C---AG-C-G--G--AT-TC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGA---A---TG---C----C-G------------G----G---C-A--CA---------------C-T-A-G-G-GG-G--AC-C-G-CCT--G-T------------------------------------G-A---TAA----------------------------------A-T-A-G--G-A-GG-A--AGG-A--GTGG-A-CGAC-GGT--AGGT-C---CGT-A-T-G-C-C-C-CGA----AT-C--CT-C-T-GG-GC-AA-CAC-GCGGG-C--TA--CAATG---G-ATGA-G-A--C-AAT-GG-GT--------------------------------------------------------------------------------------------------T-C-C-G-A--C-GCCG-A--A---------------------------------------A-GG-T-G-----------G--A-G-GT---A----------A--TCC-T------C-T-AAACT-TA-T-T-C-G-TAG-TTC--------GGA-T-TGAGG-AC--T-GTAA-CT-C-------------------------------------------------------------------------------------------------G-TTCTC-A-T-G-AA-G-CT-GGAAT-GC-G-TA--G-TA-AT-C-G-C----GTG-TC-A-T--A------AT--CGC-GC-G-GT-G-AAT-ACGT-C-CCTGCTCCT-TGGA----CACACCG-CCC-GTC-----A---CG---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

gg_id_4_ssualigned_seq = "-----------------------------------------------------------------------------------------------------------------------------------------------C-T-GC--TAT-T--G-GG-ATCC--GA---C-T--AAGCCA-T-GC-A-AGC-CCA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GG-C-AG-A--A-------------GGC-TCAGT-A--AC-AC-G-T-C-GC---TAA--C-CT-G--C-C-CTA--AG-G------------------------------------------------------------------T-CC----AGG-AT-AA-CCT-------------------------C-G-G-----------------------GAA-A---CTG-AGG-ATAA-TA---CT-G--G-AT-G---------------------------------G--G-G-A--A--A-----------------AG-ATAC-T-----------------------------------------------------------------------------------------------------------------------G-GA-A--------------------------------------------------------------------------------------------------------------------------------------T-G-TT-C-T---------------C--T-T-C-C-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G--------------C----C-T---T-AG-G---AT---G-G-----G-GCG-GCG--G-CCG--A------TT--A--G-GT-T----G---TTGG-C-G-GT-G-T----AAT-GG-A-C-C-ACCA--A-GC-C-T--A-TG-A------------TCG-G-T------AC-GG-G-C-CAT-G-AG----A--GT-GG--T-AG-C-CC-GGAG-A-AGGGT--A-C-TG-A-GA-T-AC-G-G-ACCCTA-GCCC-TAC-G--G-G-G-T-GC-A-GC-A-G-GC---GC-G-A-AAC-CTCTG-C-AA-T-GC--AC-GA-A----A-G-TG-T-GA-CA-GG-GGGA-TCCC-A-AG-T---G-----------------------------------------------------------------------C--------TT-TT-G-T--C-AAG----G------------------------------G--T--AA-G---T----A------------------------------C-C-TT-G-AC-----------AA-----------TAAGC-GGTGGG-T--AA---G--C-TGGT--GCCA--G-C---C--GCCG---C-GG--TA-AC--AC---CA-GC-GCC-GCA-A-G-TG-GGGA-C-CGC-TA-TT-A--T-T--GGGC-CTA----AA-GCAT-CC--G-TA-G-C-C-G------------G--T-CA-A-A-T-AA----A-T-C-T---T-CTG-TG-A-AA-TC--GT-TCA-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CT-G-A-AC-GG------T-G-C-A-G-A--------A--GA-C-A-C-T-G-TTT--G-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-G-G-A-C-----C-GG--GA-G-G------------C-GT-A-AG-A----GG--TATT-CGT-A-GG--GT-A-GCG-GTAAAA-TG-CGA-TAAT-C-CT-A-CGA--AG-A-CC-A-CC-TG--T--G--GC-GAA-G--G-C---G----T--C-T-TACGA------G-AA-C---------------------------------------------------------------GGCT-C-C--GA--CG-----GT-GA-GG--G-A-TGA--AG-G-C--------------C-AGGG-GAG-C-G-AATG--GG-ATTA-G-ATA-C-----CC-C-A-GTA-G-T----C-CT--G-G-CAG-T-AAA--C-CCTG-CG--AG-CT---------A-GG--T--G-T-CG-C-GC-A--T--CC----------------------------------------------------------------------------------TCC------------------------------------------------------------------------------------------------------------------------------------------------GG---G-T-G-TG--C-G-G-T-GT-C------GT--A----GA-GAA--GT-C-G-T--T--AA-GC--T----C-GCC-GCC-T-G-GG-AAG-TA---CGG-----T-C--G-C-A-A-GGC-T--GAA-ACTT-AAA---------GGAA-TTG-GCGGG-G-G-AGCA----CTA--C-A-A-GGG-GT-G--C--GG-CG-T--GC-GGT-TT-AATT-G-G-ATT-CAAC-G-CC-G-AA-A-A-TC-TC-A-CC-AGGGG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AGACGGC--AGAATGAAGGTCAGTCTAAAGGGCTTACC-TGA-CGAGCCGAG-A-G-G-T-GGTG-CA-TGG-CC--GTC-GTC-A-GC-TC---G-TG-CC-G--TGA-GG-TGT-C-CT-G-TT-AA-GT-CAGGC-AA--------C-GAG-CGA-G-ACC-T-G-TG--CC--T-ACAC--T-T-G-C-C---A-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T---------------G-T-G-T-A-GG-G--AC-T-G-CCT--G-G------------------------------------G-----AAA----------------------------------C-C-A-G--G-A-GG-A--AGG-T--ACAG-G-CAAC-GGT--AGGT-C---TGT-A-T-G-C-C-C-CGA----AT-C--CA-C-T-GG-GC-TA-CAC-GCGCG-C--AA--CAATG---G-ACGA-G-A--C-AAT-G--GC--------------------------------------------------------------------------------------------------T-G-C-A-A--C-ACCG-A--A---------------------------------------A-GG-T-G-----------A--A-G-CT---A----------A--TC-------------AAACT-CG-T-C-C-T-CAG-TAC--------GGA-T-TGAGG-CT--T-GTAA-CT-C-------------------------------------------------------------------------------------------------A-GCCTC-A-T-G-AC-G-CC-GGAAT-CC-C-TA--G-TA-AT-C-G-G----AGT-TC-A-T-C-------AT--ACT-CC-G-GT-G-AAC-ATGT-C-CCTGTTCCT-TGTT----CACACCG-CCC-GTC-----A---AA--CCA-CT-CG-A--G---CAG-G-GT-CT-GGA--T-GA-------G--G--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T--CG--AAT-C----C-AGG-CT-CAG------------------------TG--AG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
example_arb_recs = """BEGIN
gg_id=50
ncbi_acc_w_ver=AJ133622.1
ncbi_gi=5042319
db_name=
gold_id=
decision=clone
prokmsaname=
isolation_source=
clone=KTK 9A
organism=uncultured archeon 'KTK 9A'
strain=
specific_host=Uncultured archaeon 'KTK 9A' 16S rRNA gene.
authors=Eder,W., Ludwig,W. and Huber,R.
title=Novel 16S rRNA gene sequences retrieved from highly saline brine sediments of kebrit deep, red Sea
journal=Arch. Microbiol. 172 (4), 213-218 (1999)
pubmed=10525737
submit_date=24-MAR-2000
country=
ncbi_tax_string=Archaea; environmental samples
silva_tax_string=Archaea;Euryarchaeota;Thermoplasmata;Thermoplasmatales;Marine Benthic Group D and DHVEG-1
greengenes_tax_string=k__Archaea; p__Euryarchaeota; c__Thermoplasmata; o__E2; f__DHVEG-1; g__; s__
hugenholtz_tax_id=
non_acgt_percent=0.0772799998522
perc_ident_to_invariant_core=99.6440963745
warning=
aligned_seq=----------------------------------------------------------------------------------------------------------------------------------------------AC-C-GC--TCT-T--G-GG-ATTC--GA---C-T--AAGCCA-T-GC-G-AGT-CGA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GG-C-CG-A--C-------------TGC-TCAGT-A--AC-AC-G-T-G-GA---TAA--T-CT-A--C-C-CTT--AG-G------------------------------------------------------------------T-GG----AGG-AT-AA-CCT-------------------------C-G-G-----------------------GAA-A---CTG-AGG-ATAA-TA---CT-C--C-AT-A---------------------------------G--A-T-C--T--A-----------------AG-AGGC-T-----------------------------------------------------------------------------------------------------------------------G-GA-A--------------------------------------------------------------------------------------------------------------------------------------T-G-CA-C-T---------------T--A-G-G-T-C-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G--------------C----C-T---A-AG-G---AT---G-A-----G-TCT-GCG--G-ACT--A------TC--A--G-GT-T----G---TAGC-T-A-GG-G-T----AAC-GT-C-C-T-AGCT--A-GC-C-T--A-CA-A------------CGG-A-T------AC-GG-G-T-TGT-G-AG----A--GC-AA--T-AA-C-CC-GGAG-A-CGGAC--T-C-TG-A-GA-C-AC-G-A-GTCCGG-GCCC-TAC-G--G-G-G-C-GC-A-GC-A-G-GC---GC-G-A-AAC-CTTCA-C-AA-T-GT--GC-GA-A----A-G-CG-C-GA-TG-AG-GGGA-TCCC-A-AG-T---G-----------------------------------------------------------------------C--------TG-TT-T-C--T-CTG----T------------------------------C--C--AA-A---A----C------------------------------G-C-AG-A-GA-----------AG-----------TAAGG-GCTGGG-T--AA---G--C-GGGT--GCCA--G-C---C--GCCG---C-GG--TA-AT--AC---CT-GC-AGC-CCG-A-G-TG-GTGA-C-CGC-TA-TT-A--T-T--GAGC-CTA----AA-ACGT-CC--G-TA-G-C-C-G------------G--C-TT-A-G-T-AA----A-T-G-C---C-TGG-GT-A-AA-TC--GT-ATA-G--------------------------------------------------------------------CC-T-AA-------------------------------------------------------------------------CT-A-T-AC-GA------T-T-C-C-G-G--------G--CA-G-A-C-T-G-CTA--A-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-A-C-----C-GG--GA-G-A------------G-GT-T-AG-A----GG--TACT-CCT-G-GG--GT-A-GGG-GTGAAA-TC-CTG-TAAT-C-CT-G-GGG--GG-A-CC-A-CC-AG--T--G--GC-GAA-A--G-C---G----T--C-T-AACCA------T-AA-C---------------------------------------------------------------GGGT-C-T--GA--CG-----GT-GA-GG--G-A-CGA--AG-C-C--------------C-TGGG-GAG-C-A-AACC--GG-ATTA-G-ATA-C-----CC-G-G-GTA-G-T----C-CA--G-G-GTG-T-AAA--C-GCTG-CT--TG-CT---------T-NA--T--G-T-TA-G-TT-G--G--GC----------------------------------------------------------------------------------TAC------------------------------------------------------------------------------------------------------------------------------------------------GC---C-C-A-GC--T-A-G-T-GT-C------GG--A----GA-GAA--GT-T-G-T--T--AA-GC--A----A-GTT-GCT-T-G-GG-AAT-TA---CGG-----C-C--G-C-A-A-GGC-T--GAA-ACTT-AAA---------GGAA-TTG-GCGGG-G-G-AGCA----CAG--C-A-A-CGG-GT-G--G--AG-CG-T--GC-GGT-TT-AATT-G-G-ATT-CAAC-G-CC-G-GA-A-A-AC-TC-A-CC-GGAGG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CGACGGT--TATATGAGGGCCAGGCTGATGACCTTGCC-TAA-TTTTCCGAG-A-G-G-T-GGTG-CA-TGG-CC--GCC-GTC-A-GT-TC---G-TA-CC-G--TAA-GG-CGT-T-CT-G-TT-AA-GT-CAGAC-AA--------C-GAA-CGA-G-ACC-C-T-TG--CC--C-CTAG--T-T-G-C-C---C-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------A---------------C-T-A-G-G-GG-A--AC-C-G-CTG--G-C------------------------------------G-C---TAA----------------------------------G-T-C-A--G-A-GG-A--AGG-T--GAGG-G-CAAC-GGT--AGGT-C---AGT-A-T-G-C-C-C-CGA----AT-C--CT-C-C-GG-GC-TA-CAC-GCGCG-C--TA--CAAAG---G-CCAG-G-A--C-AAT-G--GC--------------------------------------------------------------------------------------------------T-C-C-G-A--C-ACCG-A--A---------------------------------------A-GG-T-G-----------A--A-G-GT---A----------A--TC-------------AAACC-TG-G-T-C-G-TAG-TTC--------GGA-T-TGAGG-GC--T-GAAA-CT-C-------------------------------------------------------------------------------------------------G-CTCTC-A-T-G-AA-G-TT-GGAAT-CC-G-TA--G-TA-AT-C-G-C----GGG-TC-A-A-C-------AA--CCC-GC-G-GT-G-AAT-A-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
END

BEGIN
gg_id=79
ncbi_acc_w_ver=AF028688.1
ncbi_gi=2599406
db_name=
gold_id=
decision=named_isolate
prokmsaname=Methanobacterium bryantii
isolation_source=
clone=
organism=Methanobacterium bryantii
strain=RiH2
specific_host=Methanobacterium bryantii strain RiH2 16S ribosomal RNA gene, complete sequence.
authors=Joulian,C., Ollivier,B., Patel,B.K.C. and Roger,P.A.
title=Phenotypic and phylogenetic characterization of dominant culturable methanogens isolated from ricefield soils
journal=FEMS Microbiol. Ecol. 25 (2), 135-145 (1998)
pubmed=
submit_date=16-MAR-2000
country=
ncbi_tax_string=Archaea; Euryarchaeota; Methanobacteria; Methanobacteriales; Methanobacteriaceae; Methanobacterium
silva_tax_string=Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobacterium
greengenes_tax_string=k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium; s__bryantii
hugenholtz_tax_id=Archaea; Methanobacteria_Eury; Methanobrevibacter; Methanobacterium_subterraneum
non_acgt_percent=0.0
perc_ident_to_invariant_core=99.4792022705
warning=
aligned_seq=-----------------------------------------------------------------------------------------------------------TTCCGGTT-GA--T-CC-C-G-CCGG-AG-GC-C-AC-T-GC--TTT-T--G-GG-GTTC--GA---T-T--AAGCCA-T-GC-A-AGT-TGA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GG-C-GA-A--C-------------GGC-TCAGT-A--AC-AC-G-T-G-GA---TAA--C-CT-G--C-C-CTT--AG-G------------------------------------------------------------------A-CT----GGG-AT-AA-CCC-------------------------T-G-G-----------------------GAA-A---CTG-GGG-ATAA-TA---CC-G--G-AT-A---------------------------------T--T-T-A--G--G-----------------AG-TTCC-T-----------------------------------------------------------------------------------------------------------------------G-GA-A--------------------------------------------------------------------------------------------------------------------------------------T-G-GT-C-T---------------T--C-T-A-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G--------------C----C-T---A-AG-G---AT---G-G-----A-TCT-GCG--G-CAG--A------TT--A--G-GT-A----G---TTGG-T-G-GG-G-T----AAT-GG-C-C-C-ACCA--A-GC-C-T--T-TG-A------------TCT-G-T------AC-GG-G-T-TGT-G-AG----A--GC-AA--G-AG-C-CC-GGAG-A-TGGAA--C-C-TG-A-GA-C-AA-G-G-TTCCAG-GCCC-TAC-G--G-G-G-C-GC-A-GC-A-G-GC---GC-G-A-AAC-CTCCG-C-AA-T-GC--AC-GA-A----A-G-TG-C-GA-CG-GG-GGGA-CCCC-A-AG-T---G-----------------------------------------------------------------------C--------TT-TT-T-T--T-GAG----T------------------------------G--T--AA-A---A----A------------------------------G-C-TT-T-GA-----------AA-----------TAAGA-GCTGGG-C--AA---G--C-CGGT--GCCA--C-C---C--GCCG---C-GG--TA-AC--AC---CG-GC-AGC-TCA-A-G-TG-GTGG-C-CAT-TT-TT-A--T-T--GGGC-CTA----AA-GCGT-TC--G-TA-G-C-C-G------------G--C-TT-G-A-T-AA----G-T-C-T---C-TGG-TG-A-AA-TC--CC-ATA-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CT-G-T-GG-GA------T-G-C-T-G-G--------A--GA-T-A-C-T-A-TTA--G-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-G-C-----C-GG--GA-G-A------------G-GT-T-AG-G----GG--TACT-CCC-A-GG--GT-A-GGG-GTGAAA-TC-CTA-TAAT-C-CT-G-GGA--GG-A-CC-A-CC-TG--T--G--GC-GAA-G--G-C---G----C--C-T-AACTG------G-AA-C---------------------------------------------------------------GGAC-C-T--GA--CG-----GT-GA-GT--A-A-CGA--AA-G-C--------------C-AGGG-GCG-C-G-AACC--GG-ATTA-G-ATA-C-----CC-G-G-GTA-G-T----C-CT--G-G-CCG-T-AAA--C-GATG-CG--GA-CT---------T-GG--T--G-T-TG-G-AA-T--G--GC----------------------------------------------------------------------------------TTC------------------------------------------------------------------------------------------------------------------------------------------------GC---T-G-C-TC--C-A-G-T-GC-C------GA--A----GG-GAA--GC-T-G-T--T--AA-GT--C----C-GCC-GCC-T-G-GG-AAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTT-AAA---------GGAA-TTG-GCGGG-G-G-AGCA----CCA--C-A-A-CGC-GT-G--G--AG-CC-T--GC-GGT-TT-AATT-G-G-ATT-CAAC-G-CC-G-GA-C-A-TC-TC-A-CC-AGGGG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CGACAGC--AGAATGATGGCCAGGTTGACGACCTTGCT-TGA-CAAGCTGAG-A-G-G-A-GGTG-CA-TGG-CC--GCC-GTC-A-GC-TC---G-TA-CC-G--TGA-GG-CGT-C-CT-G-TT-AA-GT-CAGGC-AA--------C-GAG-CGA-G-ACC-C-A-CG--CC--C-TTAG--T-T-A-C-C---A-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------A---------------C-T-A-A-G-GG-G--AC-C-G-CCA--G-T------------------------------------G-A---TAA----------------------------------A-C-T-G--G-A-GG-A--AGG-A--GTGG-A-CGAC-GGT--AGGT-C---CGT-A-T-G-C-C-C-CGA----AT-C--CC-C-T-GG-GC-TA-CAC-GCGGG-C--TA--CAATG---G-CTAT-G-A--C-AAT-G--GT--------------------------------------------------------------------------------------------------T-C-C-G-A--C-ACTG-A--A---------------------------------------A-GG-T-G-----------A--A-G-GT---A----------A--TC-------------AAACA-TA-G-T-C-T-TAG-TTC--------GGA-T-CGAGG-GC--T-GTAA-CC-T-------------------------------------------------------------------------------------------------G-CCCTC-G-T-G-AA-G-CT-GGAAT-GC-G-TA--G-TA-AT-C-G-C----GTG-TC-A-T-A-------AT--CGC-GC-G-GT-G-AAT-ACGT-C-CCTGCTCCT-TGCA----CACACCG-CCC-GTC-----A---CG--CCA-CC-CA-A--A---AAG-G-GT-TT-GGA--T-GA-------G--G--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T--CG--AAT-C----T-GGG-TT-CTT------------------------TG--AGG-AGGG-CG-AAG-TCGTAACAA-GGTAG-CCGT-AGGGGAA-CCTG-CGGC-TGGATCACCTCCT------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
END

BEGIN
gg_id=86
ncbi_acc_w_ver=U55237.1
ncbi_gi=3201648
db_name=
gold_id=
decision=named_isolate
prokmsaname=Methanobrevibacter woesei
isolation_source=
clone=
organism=Methanobrevibacter woesei
strain=GS
specific_host=Methanobrevibacter woesei strain GS 16S ribosomal RNA gene, complete sequence.
authors=Lin,C. and Miller,T.L.
title=Phylogenetic analysis of Methanobrevibacter isolated from feces of humans and other animals
journal=Arch. Microbiol. 169 (5), 397-403 (1998)
pubmed=9560420
submit_date=06-JUN-2002
country=
ncbi_tax_string=Archaea; Euryarchaeota; Methanobacteria; Methanobacteriales; Methanobacteriaceae; Methanobrevibacter
silva_tax_string=Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter
greengenes_tax_string=k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobrevibacter; s__
hugenholtz_tax_id=Archaea; Methanobacteria_Eury; Methanobrevibacter; Methanobrevibacter_cuticularis
non_acgt_percent=0.074963003397
perc_ident_to_invariant_core=99.2188034058
warning=
aligned_seq=-----------------------------------------------------------------------------------------------------------------------------------------GC-T-AC-T-GC--TAT-T--G-GG-ATTC--GA---T-T--AAGCCA-T-NC-A-AGT-TGA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GG-C-GT-A--C-------------GGC-TCAGT-A--AC-AC-G-T-G-GA---TAA--C-CT-A--C-C-CTT--AG-G------------------------------------------------------------------A-CT----GGG-AT-AA-CTC-------------------------T-G-G-----------------------GAA-A---CTG-GGG-ATAA-TA---CC-G--G-AT-A---------------------------------G--G-C-A--A--T-----------------TT-TTCC-T-----------------------------------------------------------------------------------------------------------------------G-TA-A--------------------------------------------------------------------------------------------------------------------------------------T-G-GT-T-T---------------T--T-T-G-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G--------------C----C-T---A-AG-G---AT---G-G-----G-TCT-GCG--G-CAG--A------TT--A--G-GT-A----G---TTGG-T-T-AG-G-T----AAT-GG-C-T-T-ACCA--A-GC-C-G--T-TG-A------------TCT-G-T------AC-GG-G-T-TGT-G-AG----A--GC-AA--G-AG-C-CC-GGAG-A-TGGAA--C-C-TG-A-GA-C-AA-G-G-TTCCAG-GCCC-TAC-G--G-G-G-C-GC-A-GC-A-G-GC---GC-G-A-AAC-CTCCG-C-AA-T-GT--GA-GA-A----A-T-CG-C-GA-CG-GG-GGGA-TCCC-A-AG-T---G-----------------------------------------------------------------------C--------TT-TT-C-A--T-TAG----T------------------------------G--T--AA-A---G----A------------------------------G-C-TT-T-TG-----------AA-----------TAAGA-GCTGGG-C--AA---G--C-CGTT--GCCA--A-C---C--GCCG---C-GG--TA-AC--AC---CG-TC-AGC-TCT-A-G-TG-GTAG-C-AGT-TT-TT-A--T-T--GGGC-CTA----AA-GCGT-CC--G-TA-G-C-C-G------------G--T-TT-A-T-T-AA----G-T-C-T---C-TGG-TG-A-AA-TC--CT-GTA-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CT-G-T-GG-GA------T-G-C-T-G-G--------A--GA-T-A-C-T-A-GTA--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-A-T-----C-GG--GA-G-A------------G-GT-T-AG-A----GG--TACT-CCC-A-GG--GT-A-GAG-GTGAAA-TT-CTG-TAAT-C-CT-G-GGA--GG-A-CC-G-CC-TG--T--G--GC-GAA-G--G-C---G----T--C-T-AGCTG------G-AA-C---------------------------------------------------------------GATT-C-T--GC--CG-----GT-GA-GG--G-A-CGA--AA-G-C--------------T-AGGG-GCG-C-G-AACC--GG-ATTA-G-ATA-C-----CC-G-G-GTA-G-T----C-CT--A-G-CTG-T-AAA--C-GATG-CG--GA-CT---------T-GG--T--G-T-TG-G-GA-T--G--GC----------------------------------------------------------------------------------TTT------------------------------------------------------------------------------------------------------------------------------------------------GC---T-G-C-TC--C-A-G-T-GC-C------GA--A----GG-GAA--GC-T-G-T--T--AA-GT--C----C-GCC-GCC-T-G-GG-AAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTT-AAA---------GGAA-TTG-GCGGG-G-G-AGCA----CCA--C-A-A-CGC-GT-G--G--AG-CC-T--GC-GGT-TT-AATT-G-G-ATT-CAAC-G-CC-G-GA-C-A-TC-TC-A-CC-AGAGG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CGACAGC--TGTATGATGACCAGTTTGATGAGCTTGTT-TGA-CTAGCTGAG-A-G-G-A-GGTG-CA-TGG-CC--GCC-GTC-A-GC-TC---G-TA-CC-G--TGA-GG-CGT-C-CT-G-TT-AA-GT-CAGGC-AA--------C-GAG-CGA-G-ACC-C-A-CG--CC--C-TTAG--T-T-A-C-C---A-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------A---------------C-T-A-G-G-GG-G--AC-C-G-CCT--G-T------------------------------------G-A---TAA----------------------------------A-T-A-G--G-A-GG-A--AGG-A--GTGG-A-CGAC-GGT--AGGT-C---CGT-A-T-G-C-C-C-CGA----AT-C--CT-C-T-GG-GC-AA-CAC-GCGGG-C--TA--CAATG---G-ATGA-G-A--C-AAT-G--GT--------------------------------------------------------------------------------------------------T-C-C-G-A--C-GCCG-A--A---------------------------------------A-GG-T-G-----------G--A-G-GT---A----------A--TC-------------AAACT-TA-T-T-C-G-TAG-TTC--------GGA-T-TGAGG-AC--T-GTAA-CT-C-------------------------------------------------------------------------------------------------G-TTCTC-A-T-G-AA-G-CT-GGAAT-GC-G-TA--G-TA-AT-C-G-C----GTG-TC-A-T-A-------AT--CGC-GC-G-GT-G-AAT-ACGT-C-CCTGCTCCT-TGGA----CACACCG-CCC-GTC-----A---CG---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
END

"""

if __name__ == '__main__':
    main()
