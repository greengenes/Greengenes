#!/usr/bin/env python

from greengenes.pick_chimeras import get_overlap, determine_taxon_conflict
from greengenes.parse import parse_b3_chimeras, parse_cs_chimeras
from t2t.nlevel import load_consensus_map
from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

script_info={}
script_info['brief_description']="""Determine chimeric accessions"""
script_info['script_description']="""Consume chimeras from bellerophon and ChimeraSlayer. A chimera is determined to be all sequences that were found chimeric by both methods or if the sequence introduces class level variation in the reference database"""
script_info['script_usage']=[]
script_info['required_options'] = [\
        make_option('-c','--input-cs',type='str',
            help="ChimeraSlayer results"),
        make_option('-b','--input-bellerophon',type='str',
            help="Bellerophon results"),
        make_option('-o','--output',type='str',
            help="Output directory"),
        make_option('-m','--ref-taxonomy-map',type='str',
            help='Taxonomy strings for the reference database')]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    taxlookup = load_consensus_map(open(opts.ref_taxonomy_map), False)
    cs_results = parse_cs_chimeras(open(opts.input_cs))
    b3_results = parse_b3_chimeras(open(opts.input_bellerophon))
    
    output = open(opts.output,'w')
    output.write("#accession\treason\tnote\tnote\n")
    overlap = get_overlap(b3_results, cs_results)
    for id_ in overlap:
        output.write("%s\tFound by both Bellerophon and ChimeraSlayer\n" % id_)

    for id_, score, parent_a, parent_b in b3_results:
        if id_ in overlap:
            continue
        if determine_taxon_conflict(taxlookup, parent_a, parent_b):
            o = [id_,"Class conflict found by Bellerophon"]
            o.append("%s: %s" % (parent_a, '; '.join(taxlookup[parent_a])))
            o.append("%s: %s" % (parent_b, '; '.join(taxlookup[parent_b])))
            output.write('\t'.join(o))
            output.write('\n')

    for id_, parent_a, parent_b in cs_results:
        if id_ in overlap:
            continue
        if determine_taxon_conflict(taxlookup, parent_a, parent_b):
            o = [id_,"Class conflict found by ChimeraSlayer"]
            o.append("%s: %s" % (parent_a, '; '.join(taxlookup[parent_a])))
            o.append("%s: %s" % (parent_b, '; '.join(taxlookup[parent_b])))
            output.write('\t'.join(o))
            output.write('\n')

if __name__ == '__main__':
    main()
