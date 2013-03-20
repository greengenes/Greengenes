#!/usr/bin/env python

from greengenes.pick_chimeras import get_overlap, determine_taxon_conflict
from greengenes.parse import parse_uchime_chimeras
from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

script_info={}
script_info['brief_description']="""Determine chimeric accessions"""
script_info['script_description']="""Consume chimeras from UCHIME"""
script_info['script_usage']=[]
script_info['required_options'] = [\
        make_option('-c','--input-uchime',type='str',
            help="UCHIME results"),
        make_option('-o','--output',type='str',
            help="Output directory"),
        make_option('-m','--ref-taxonomy-map',type='str',
            help='Taxonomy strings for the reference database')]
script_info['version'] = __version__

def load_consensus_map(lines):
    """returns {id_:[taxon,names,...]}"""
    res = {}
    for l in lines:
        fields = l.strip().split('\t')
        taxons = fields[1].split('; ')
        
        assert len(fields) == 2
        assert len(taxons) == 7
        
        res[fields[0]] = taxons
    return res


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    taxlookup = load_consensus_map(open(opts.ref_taxonomy_map))
    uchime_results = parse_uchime_chimeras(open(opts.input_uchime))
    
    output = open(opts.output,'w')
    output.write("#accession\treason\tnote\tnote\n")

    for id_, score, parent_a, parent_b in uchime_results:
        if determine_taxon_conflict(taxlookup, parent_a, parent_b):
            o = [id_,"Class conflict found by UCHIME"]
            o.append("%s: %s" % (parent_a, '; '.join(taxlookup[parent_a])))
            o.append("%s: %s" % (parent_b, '; '.join(taxlookup[parent_b])))
            output.write('\t'.join(o))
            output.write('\n')

if __name__ == '__main__':
    main()
