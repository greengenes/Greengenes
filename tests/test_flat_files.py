#!/usr/bin/env

from cogent.util.unit_test import TestCase,main
from cogent.parse.genbank import MinimalGenbankParser
from greengenes.flat_files import get_accession, get_gi, get_sequence, \
        get_decision, get_isolation_source, get_organism, get_taxon, \
        get_country, get_ncbi_taxonomy, get_gold_id, _parse_migs_poorly, \
        get_title, get_journal, get_authors, get_pubmed, get_taxon, \
        get_ncbi_taxonomy, get_country, get_genbank_summary, get_strain
        #get_gb_summary_header, field_order
from StringIO import StringIO
from greengenes.util import GreengenesRecord

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class FlatFilesTests(TestCase):
    def setUp(self):
        self.gb1 = MinimalGenbankParser(StringIO(id_AGIY01000001_1_gb)).next()
        self.gb2 = MinimalGenbankParser(StringIO(id_FO117587_1_gb)).next()
        self.multi = MinimalGenbankParser(StringIO(combined))

    def test_get_accession(self):
        """Get an accession"""
        exp_gb1 = "AGIY01000001.1"
        exp_gb2 = "FO117587.1"

        obs_gb1 = get_accession(self.gb1)
        obs_gb2 = get_accession(self.gb2)

        self.assertEqual(obs_gb1, exp_gb1)
        self.assertEqual(obs_gb2, exp_gb2)

    def test_get_gi(self):
        """Get the GI"""
        exp_gb1 = "354825968"
        exp_gb2 = "376316277"
   
        obs_gb1 = get_gi(self.gb1)
        obs_gb2 = get_gi(self.gb2)

        self.assertEqual(obs_gb1, exp_gb1)
        self.assertEqual(obs_gb2, exp_gb2)
    
    def test_get_sequence(self):
        """get the sequence"""
        exp_gb1_start = "CCATGATTCGACCATTTTCAGAGAG"
        exp_gb1_end =   "CAACGGTCAGGCCAG"
        exp_gb1_len = 452612
        
        obs = get_sequence(self.gb1)
        obs_gb1_start = obs[:25]
        obs_gb1_end = obs[-15:]
        obs_gb1_len = len(obs)

        self.assertEqual(obs_gb1_start, exp_gb1_start)
        self.assertEqual(obs_gb1_end, exp_gb1_end)
        self.assertEqual(obs_gb1_len, exp_gb1_len)

    def test_get_decision(self):
        """get the decision"""
        self.fail("not sure yet how to get this")

    def test_get_isolation_source(self):
        """Try to get the isolation source"""
        exp = "anaerobic digested sludge"
        obs = get_isolation_source(self.gb1)
        self.assertEqual(obs, exp)
    
    def test_get_gold_id(self):
        """Get a gold id!"""
        exp = "Gi05850"
        obs = get_gold_id(self.gb1)
        self.assertEqual(obs,exp)

    def test_get_organism(self):
        """get the organism"""
        exp_gb1 = "Methanolinea tarda NOBI-1"
        obs_gb1 = get_organism(self.gb1)
        self.assertEqual(obs_gb1,exp_gb1)

    def test_gzopen(self):
        """not testing... explicit filesystem op"""
        pass

    def test_parse_migs_poorly(self):
        """Parse MIGS data... assuming its store the same way in GB recs"""
        f = lambda x: x.strip().split()[0]
        exp_it = "bacteria_archaea"
        exp_missing = None
        obs_it = _parse_migs_poorly(self.gb1['comment'], "investigation_type", f)
        obs_missing = _parse_migs_poorly(self.gb1['comment'],"w00p not here")
        self.assertEqual(obs_it, exp_it)
        self.assertEqual(obs_missing, exp_missing)

    def test_get_title(self):
        """Gets a title"""
        exp = "The draft genome of Methanolinea tarda NOBI-1"
        obs = get_title(self.gb1)
        self.assertEqual(obs,exp)

    def test_get_pubmed(self):
        """Gets pubmed if available"""
        exp = None
        obs = get_pubmed(self.gb1)
        self.assertEqual(obs,exp)

        exp = "21895912"
        obs = get_pubmed(self.gb2)
        self.assertEqual(obs,exp)

    def test_get_strain(self):
        """gets the strain"""
        exp = "NOBI-1"
        obs = get_strain(self.gb1)
        self.assertEqual(obs,exp)

    def test_get_genbank_summary(self):
        """Get the summary!!"""
        exp = GreengenesRecord({'ncbi_acc_w_ver':'AGIY01000001.1',
                'ncbi_gi':'354825968',
                'db_name':'NOT SURE YET',
                'gold_id':'Gi05850',
                'decision':'NOT SURE YET',
                'isolation_source':'anaerobic digested sludge',
                'clone':'NOT SURE YET',
                'organism':'Methanolinea tarda NOBI-1',
                'strain':'NOBI-1',
                'specific_host':'NOT SURE YET',
               'authors':'Lucas,S., Han,J., Lapidus,A., Cheng,J.-F., Goodwin,L., Pitluck,S., Peters,L., Land,M.L., Hauser,L., Imachi,H., Sekiguchi,Y., Kamagata,Y., Cadillo-Quiroz,H., Zinder,S., Liu,W.T., Tamaki,H. and Woyke,T.J.',
               'title':'The draft genome of Methanolinea tarda NOBI-1',
               'submit_date':'NOT SURE YET',
               'country':'Japan: Nagaoka',
               #'NCBI_tax_id':'882090',
               'ncbi_tax_string':'Archaea; Euryarchaeota; Methanomicrobia; Methanomicrobiales; Genera incertae sedis; Methanolinea'})
        obs = get_genbank_summary(self.gb1)
        
        for k in obs:
            if obs[k] != exp[k]:
                print k
                print obs[k]
                print exp[k]
        self.assertEqual(obs,exp)

    def test_get_authors(self):
        """Gets the authors"""
        exp = "Lucas,S., Han,J., Lapidus,A., Cheng,J.-F., Goodwin,L., Pitluck,S., Peters,L., Land,M.L., Hauser,L., Imachi,H., Sekiguchi,Y., Kamagata,Y., Cadillo-Quiroz,H., Zinder,S., Liu,W.T., Tamaki,H. and Woyke,T.J."
        obs = get_authors(self.gb1)
        self.assertEqual(obs, exp)

    def test_get_journal(self):
        """Get the journal"""
        exp = None
        obs = get_journal(self.gb1)
        self.assertEqual(obs,exp)

        exp = "Environ. Microbiol. 14 (1), 52-66 (2012)"
        obs = get_journal(self.gb2)
        self.assertEqual(obs,exp)

    def test_get_taxon(self):
        """get the taxon"""
        exp = "882090"
        obs = get_taxon(self.gb1)
        self.assertEqual(obs, exp)

        self.gb1['features'][0]['db_xref'] = ""
        exp = ""
        obs = get_taxon(self.gb1)
        self.assertEqual(obs, exp)

    def test_get_country(self):
        """get the country"""
        exp = "Japan: Nagaoka"
        obs = get_country(self.gb1)
        self.assertEqual(obs,exp)
    
    def test_get_ncbi_taxonomy(self):
        """Get the taxonomy string"""
        exp = "Archaea; Euryarchaeota; Methanomicrobia; Methanomicrobiales; Genera incertae sedis; Methanolinea"
        obs = get_ncbi_taxonomy(self.gb1)
        self.assertEqual(obs,exp)

id_AGIY01000001_1_gb = open('test_data/AGIY01000001.1.gb').read()
id_FO117587_1_gb = open('test_data/FO117587.1.gb').read()
combined = open('test_data/combined.gb').read()

if __name__ == '__main__':
    main()

