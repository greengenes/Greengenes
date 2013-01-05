#!/usr/bin/env python

"""Tests for methods that pull down likely 16S from NCBI"""

from StringIO import StringIO
from cogent.parse.genbank import RichGenbankParser
from greengenes.ncbi import esearch, parse_esearch, bulk_efetch, parse_gi_from_gb
from cogent.util.unit_test import TestCase, main

class NCBITests(TestCase):
    def setUp(self):
        pass

    def test_esearchs(self):
        """Get the first batch of 16S GIs
        
        this test isn't reliable as NCBI does not always return the same GIs as
        more records are added. Just make sure we get back 10 GI like ids
        """
        # only pulling 10 for testing simplicity
        obs = esearch('16S', retmax=10, binsize=3)
        try:
            foo = map(int, obs)
        except:
            self.fail()

    def test_parse_esearch(self):
        """hate xml"""
        exp = set(["381140089", "378760275", "376316475", "373455715",
                   "373454887", "373454086", "365275081", "359254912",
                   "354828114", "354825968",])
        obs = parse_esearch(example_xml)
        self.assertEqual(obs,exp)

    def test_bulk_efetch(self):
        """efetch a set of ids"""
        exp = [''.join([id_AGIY01000001_1_gb, id_FO117587_1_gb])]
        obs = list(bulk_efetch(['354825968','FO117587.1']))

        # ncbi records are dynamic even with the same accession versions. so,
        # a test here is not reliable.
        self.assertEqual(exp[0][:100], obs[0][:100])
        self.assertEqual(obs[0].count('//\n'), 2)

    def test_parse_gi_from_gb(self):
        """parse out gi numbers from gb records"""
        exp = set(['354825968', '376316277'])
        obs = parse_gi_from_gb(combined)
        self.assertEqual(obs,exp)

example_xml = """<?xml version="1.0" ?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult><Count>4895587</Count><RetMax>10</RetMax><RetStart>0</RetStart><QueryKey>1</QueryKey><WebEnv>NCID_1_14507967_130.14.22.101_9001_1332798701_136474237</WebEnv><IdList>
        <Id>381140089</Id>
        <Id>378760275</Id>
        <Id>376316475</Id>
        <Id>373455715</Id>
        <Id>373454887</Id>
        <Id>373454086</Id>
        <Id>365275081</Id>
        <Id>359254912</Id>
        <Id>354828114</Id>
        <Id>354825968</Id>
    </IdList><TranslationSet/><TranslationStack>   <TermSet>    <Term>rrna[fkey]</Term>    <Field>fkey</Field>    <Count>0</Count>    <Explode>Y</Explode>   </TermSet>   <OP>GROUP</OP>  </TranslationStack><QueryTranslation>rrna[fkey]</QueryTranslation></eSearchResult>"""

id_AGIY01000001_1_gb = open('test_data/AGIY01000001.1.gb').read()
id_FO117587_1_gb = open('test_data/FO117587.1.gb').read()
combined = open('test_data/combined.gb').readlines()

if __name__ == '__main__':
    main()
