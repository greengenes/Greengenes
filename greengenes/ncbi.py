#!/usr/bin/env python

"""Pull down records from NCBI that may contain 16S or 18S"""

from StringIO import StringIO
from cogent.parse.genbank import GbFinder
from cogent.db.ncbi import ESearch, EUtils
from cogent.db.util import make_lists_of_accessions_of_set_size
from xml.dom.minidom import parseString

def parse_esearch(s):
    """Returns a set of ids in ESearch results or the empty set if no results"""
    ids = []

    for record in parseString(s).getElementsByTagName('Id'):
        ids.append(str(record.childNodes[0].data))

    return set(ids)

def esearch(query, retmax, binsize=10000):
    """Wraps ESearch call"""
    ids = set([])
    bins = (retmax / binsize) + 1 # handle remainder, just overshoot it
    for i in range(bins):
        handle = ESearch(term=query, db='nucleotide', retmax=binsize, 
                         retstart=binsize * i)
        ids.update(parse_esearch(handle.read()))
    return ids

def bulk_efetch(query_ids):
    """Wraps EUtils call"""
    # pre-bin calls because we're potentially obtaining millions of records
    # and EUtils would normally store all the records in memory
    bins = make_lists_of_accessions_of_set_size(list(query_ids))
    handle = EUtils(db='nucleotide',rettype='gb')

    for queries in bins:
        data = handle[queries].read()
        yield data 

def parse_gi_from_gb(lines, verbose=False):
    """Returns a set of GIs from GB records"""
    GIs = []
    count = 0
    for l in lines:
        if not l.startswith('VERSION'):
            continue

        if verbose and (count % 1000 == 0):
            print "...%i existing records..." % count
        count += 1

        fields = l.strip().split()
        
        assert fields[0] == 'VERSION'
        assert fields[2].startswith('GI:')
    
        GIs.append(fields[2].split(':')[1])
    return set(GIs)
