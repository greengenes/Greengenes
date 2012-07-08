#!/usr/bin/env python

from cogent.seqsim.tree import RangeNode
from gzip import open as gzopen
from datetime import datetime
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class NoSequenceError(Exception):
    pass

def prune_tree(tree, ids):
    """Prunes a tree to just the ids specified"""
    tcopy = tree.deepcopy()
    all_tips = set([n.Name for n in tcopy.tips()])
    ids = set(ids)

    if not ids.issubset(all_tips):
        raise ValueError, "ids are not a subset of the tree!"
    
    while len(tcopy.tips()) != len(ids):
        for n in tcopy.tips():
            if n.Name not in ids:
                n.Parent.removeNode(n)
        tcopy.prune() 
    return tcopy

def make_tree_arb_safe(t):
    """Gives a second child to all single descendent nodes"""
    for n in t.nontips(include_self=True):
        if len(n.Children) == 1:
            n.append(RangeNode(Name="X",Length=0.0))

def log_f(line):
    start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
    return "%s\t%s\n" % (start_time, line)

def greengenes_open(file_fp, permission='U'):
    """Read or write the contents of a file
    
    file_fp : file path
    permission : either 'U','r','w','a'
    
    NOTE: univeral line breaks are always used, so 'r' is automatically changed
    into 'U'
    """
    if permission not in ['U','r','w','a']:
        raise IOError, "Unknown permission: %s" % permission

    if file_fp.endswith('gz'):
        # gzip doesn't support Ub
        if permission == 'U':
            permission = 'r'
        return gzopen(file_fp, permission)
    else:
        if permission == 'r':
            permission = 'U'
        return open(file_fp, permission)

def generate_log_fp(output_dir, basefile_name='log', suffix='txt',
                    timestamp_pattern='%Y%m%d%H%M%S'):
    """Originally from the QIIME project under qiime.workflow"""
    timestamp = datetime.now().strftime(timestamp_pattern)
    filename = '%s_%s.%s' % (basefile_name,timestamp,suffix)
    return os.path.join(output_dir,filename)

class WorkflowLogger(object):
    """Originally from the QIIME project under qiime.workflow"""
    def __init__(self,log_fp=None,open_mode='w',script_name="Not specified"):
        if log_fp:
            self._f = open(log_fp,open_mode)
        else:
            self._f = None 
        start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('Logging started at %s\n' % start_time)
        self.write('Script: %s\n' % script_name)
        self.write('Greengenes workflow version: %s\n\n' % __version__)

    def write(self,s):
        if self._f:
            self._f.write(s)
            # Flush here so users can see what step they're
            # on after each write, since some steps can take
            # a long time, and a relatively small amount of 
            # data is being written to the log files.
            self._f.flush()
        else:
            pass 
    
    def close(self):
        end_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('\nLogging stopped at %s\n' % end_time)
        if self._f:
            self._f.close()
        else:
            pass

    def __del__(self):
        """Destructor"""
        self.close()

class GreengenesRecord(dict):
    """Represent a full Greengenes record"""

    _field = {
        'prokMSA_id':{'type':int, 'desc':"unique Greengenes identifier",
                    'arb_rule':'name', 'required':True},           
        'ncbi_acc_w_ver':{'type':str, 'desc':"Genbank accession with version",
                    'arb_rule':"acc", 'required':True},
        'ncbi_gi':{'type':str, 'desc':"Genbank GI",'arb_rule':'ncbi_gi',
                    'required':True},
        'db_name':{'type':str, 
                    'desc':"Integrated Microbial Genomes database identifier",
                    'arb_rule':'IMG_taxon_OID','required':False},
        'gold_id':{'type':str, 'desc':"Genomes Online Database ID",
                    'arb_rule':'GOLD_stamp','required':False},
        'decision':{'type':str, 'desc':"Whether a clone or an isolate",
                    'arb_rule':'sequence_type','required':True},
        'prokMSAname':{'type':str, 
                    'desc':"[isolation_source] clone [clone] [organism]",
                    'arb_rule':'full_name','required':False},
        'isolation_source':{'type':str, 
                    'desc':"Genbank source features, isolation_source",
                    'arb_rule':'isolation_source','required':True},
        'clone':{'type':str, 'desc':"??????",'arb_rule':'clone',
                    'required':True},
        'organism':{'type':str, 'desc':"Genbank source features, organism",
                    'arb_rule':'organism','required':False},
        'strain':{'type':str, 'desc':"Genbank source features, strain",
                    'arb_rule':'strain','required':False},
        'specific_host':{'type':str, 'desc':"?????????",
                    'arb_rule':'specific_host','required':False},
        'authors':{'type':str, 'desc':"Associated authors",
                    'arb_rule':'author','required':True},
        'title':{'type':str, 'desc':"Associated study title",
                    'arb_rule':'title','required':True},
        'journal':{'type':str, 'desc':"Associated Journal",
                    'arb_rule':'journal','required':True},
        'pubmed':{'type':int, 'desc':"Associated PUBMED ID",
                    'arb_rule':'pubmed','required':False},
        'study_id':{'type':str, 'desc':"Unique Greengenes study identifier",
                    'arb_rule':'study_id','required':True},
        'submit_date':{'type':str, 'desc':"Submission date to Genbank",
                    'arb_rule':'submit_date','required':True},
        'country':{'type':str, 'desc':"Source country",'arb_rule':'country',
                    'required':True},
        'ncbi_tax_string':{'type':str, 
                    'desc':"NCBI Taxonomy string present in Genbank record",
                    'arb_rule':'ncbi_tax','required':False},
        'Silva_tax_string':{'type':str, 
                    'desc':"Versioned Silva taxonomy string",
                    'arb_rule':'Silva_tax','required':False},
        'RDP_tax_string':{'type':str, 'desc':"Versioned RDP taxonomy string",
                    'arb_rule':'RDP_tax','required':False},
        'greengenes_tax_string':{'type':str, 
                    'desc':"Versioned Greengenes taxonomy string",
                    'arb_rule':'greengenes_tax','required':False},
        'non_ACGT_percent':{'type':float, 
                    'desc':"Percent non-ACGT bases in sequence",
                    'arb_rule':'percent_non_ACGT','required':True},
        'perc_ident_to_invariant_core':{'type':float, 
                    'desc':"Percent not aligning to invariant bases",
                    'arb_rule':'percent_invariant_match','required':True},
        'small_gap_intrusions':{'type':float, 'desc':"????????????????",
                    'arb_rule':'small_gap_intrusions','required':True},
        'bellerophon':{'type':str, 'desc':"Bellerophon determination",
                    'arb_rule':'bellerophon','required':False},
        'bel3_div_ratio':{'type':float, 'desc':"Bellerophon divergence ratio",
                    'arb_rule':'bel_div_ratio','required':False},
        'chim_slyr_a':{'type':str, 
                    'desc':"ChimeraSlayer parent sequence identifier",
                    'arb_rule':'chim_slyr_a','required':False},
        'chim_slyr_b':{'type':str, 
                    'desc':"ChimeraSlayer parent sequence identifier",
                    'arb_rule':'chim_slyr_b','required':False},
        'chim_slyr_a_tax':{'type':str, 
                    'desc':"Greengenes consensus of parent sequence",
                    'arb_rule':'chim_slyr_a_tax','required':False},
        'chim_slyr_b_tax':{'type':str, 
                    'desc':"Greengenes consensus of parent sequence",
                    'arb_rule':'chim_slyr_b_tax','required':False},
        'aligned_seq':{'type':str, 'desc':"Full unaligned sequence",
                    'arb_rule':'aligned_seq','required':True},
        'unaligned_seq':{'type':str, 
                    'desc':"SSU-Align'ed and masked sequence in NAST width",
                    'arb_rule':'unaligned_seq','required':True}
        }

 
    _arb_rule_fmt = '\n'.join(['MATCH   "%s\=*"',
                                   '\tSRT "*\=="',
                                   '\tWRITE "%s"'])

    def __init__(self, *args, **kwargs):
        for k in self._field:
            self[k] = None
        super(GreengenesRecord, self).__init__(*args, **kwargs)

    def getARBRules(self):
        """Get the ARB rules"""
        return '\n\n'.join([self._arb_rule_fmt % (f, self._field[f]['arb_rule']) \
                for f in self._field])

    def __str__(self):
        """Print out a Greengenes record"""
        return self.toGreengenesFormat()

    def setTypes(self):
        """Coerce values to expected types"""
        for k in self._field:
            if self[k] is not None:
                self[k] = self._field[k]['type'](self[k])

    def toGreengenesFormat(self):
        """Print a Greengenes record
        
        All records are forced to str. Types are not verified.
        """
        out = ["BEGIN\n"]
        for f in self._field:
            v = self[f]
            if v is None:
                v = ''
            else:
                v = str(v)

            out.append("%s=%s\n" % (f,v))
        out.append("END\n\n")
        return ''.join(out)

    def sanityCheck(self):
        """Make sure the types are as expected"""
        for k,v in self.items():
            if v is None or '':
                continue
            if not isinstance(v, self._field[k]['type']):
                raise ValueError, "Field %s in prok %s has a bad type" % \
                                        (k, str(self['prokMSA_id']))
    
