#!/usr/bin/env python

from greengenes.write import write_sequence
from greengenes.util import NoSequenceError, GreengenesRecord

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def get_accession(r):
    """Get the accession"""
    return r['version'].split()[1]

def get_sequence(r):
    """Get the sequence"""
    seq = r.get('sequence', "").upper()
    if len(seq) == 0:
        raise NoSequenceError
    else:
        return seq

def get_gi(r):
    """Get the GI"""
    return r['version'].split(':')[-1]

def get_decision(r):
    # maybe r['definition'
    return "NOT SURE YET"

def get_organism(r):
    """Get the organism"""
    res = r['features'][0].get('organism', None)
    if res:
        res = res[0]
    return res

def get_taxon(r):
    """get the ncbi taxon id"""
    res = r['features'][0].get('db_xref',None)
    if res and res[0]:
        res = res[0].split(':')[1]
    return res

def get_country(r):
    res = r['features'][0].get('country', None)
    if res:
        res = res[0]
    return res

def _parse_migs_poorly(line, term, parse_f=lambda x: x.strip()):
    """Parse MIGS"""
    if '##MIGS-Data-START##' not in line:
        return None

    chunk = line.split("##MIGS-Data-START##",1)[1].split('##MIGS-Data-END##')[0]
    chunks = chunk.split('::')
    found = -1
    for idx,v in enumerate(chunks):
        v = v.strip()
        if v == term:
            found = idx
            break
    if found != -1:
        return parse_f(chunks[found + 1])
    else:
        return None
    
def get_gold_id(r):
    """Attempt to get the GOLD id"""
    f = lambda x: x.strip().split()[0]
    if 'comment' in r:
        return _parse_migs_poorly(r['comment'], 'Finished GOLD Stamp ID', f)
    else:
        return None

def get_isolation_source(r):
    """Get the Isolation Site"""
    #return _parse_migs_poorly(r['comment'], 'Isoloation Site')
    res = r['features'][0].get('isolation_source', None)
    if res:
        res = res[0]
    return res

def get_dbname(r):
    return "NOT SURE YET"

def get_clone(r):
    return "NOT SURE YET"

def get_strain(r):
    """Get the strain"""
    res = r['features'][0].get('strain', None)
    if res:
        res = res[0]
    return res

def get_specific_host(r):
    return "NOT SURE YET"

def get_title(r):
    """Get pub title"""
    # if the genbank record doesn't have a reference section, can I just call
    # it a fail?
    if 'references' not in r:
        return None
    
    return r['references'][0].get('title',None)
    

def get_authors(r):
    """Get all the authors"""
    if 'references' not in r:
        return None
    
    return r['references'][0].get('authors', None)

def get_journal(r):
    """Get journal or empty"""
    if 'references' not in r:
        return None

    res = r['references'][0].get('journal', None)
    if res.lower() == 'unpublished':
        return None
    else:
        return res

def get_pubmed(r):
    """Get PUBMED if available"""
    if 'references' not in r:
        return None
    
    return r['references'][0].get('pubmed',None)

def get_submit_date(r):
    return "NOT SURE YET"

def get_ncbi_taxonomy(r):
    return '; '.join(r['taxonomy'])

def get_genbank_summary(r):
    """Get the gb summary data"""
    rec = GreengenesRecord()

    for f,m in parse_funs:
        rec[f] = m(r)

    return rec
#    return [f(r) for field, f in parse_funs]

#def get_gb_summary_header():
#    return field_order

parse_funs = [('ncbi_acc_w_ver', get_accession),
              ('ncbi_gi', get_gi),
              ('db_name', get_dbname),
              ('gold_id', get_gold_id),
              ('decision', get_decision),
              ('isolation_source', get_isolation_source),
              ('clone', get_clone),
              ('organism', get_organism), 
              ('strain', get_strain),
              ('specific_host', get_specific_host),
              ('authors', get_authors),
              ('title', get_title),
              ('journal', get_journal),
              ('pubmed', get_pubmed),
              ('submit_date', get_submit_date),
              ('country', get_country),
              #('NCBI_tax_id', get_taxon),
              ('ncbi_tax_string', get_ncbi_taxonomy)]

# functions requiring information external 
#external_funs = [('Silva_tax_string', get_silva_taxonomy),
#                 ('RDP_tax_string', get_rdp_taxonomy)]
field_order = [a for a,b in parse_funs]

for f,m in parse_funs:
    if f not in GreengenesRecord._field:
        raise KeyError, "%s is not a valid field" % f
