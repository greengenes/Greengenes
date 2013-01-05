#!/usr/bin/env python

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters
from greengenes.ncbi import esearch, parse_esearch, bulk_efetch, parse_gi_from_gb
from gzip import open as open_gz
import time

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

script_info = {}
script_info['script_description'] = "Collect 16S records from NCBI"
script_info['brief_description'] = "Collecting 16S records from NCBI"
script_info['script_usage'] = [("","","")]
script_info['required_options'] = [make_option('--possible-new-gb-output',
                                    dest='possible_new_gb_out', 
                                    help='Destination for POSSIBLE new rRNA containing gb records')]
script_info['optional_options'] = [make_option('--existing-gb',dest='existing_gb', 
                                    help='A file containing existing Genbank records that to not need to be re-discovered'),
        make_option('--use-gz',dest='use_gz',action='store_true',default=False, help='Expects inputs to be .gz and writes output as .gz'),
        make_option('--cached-ids-to-obtain', dest='cached_ids', help="Use these IDs for records instead of querying NCBI")]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # if we already have these records, then we do not need to reobtain them
    if opts.existing_gb:
        existing_gis = set([l.strip() for l in open(opts.existing_gb)])
    else:
        existing_gis = set([])

    if opts.verbose:
        print "Number of existing GIs: %d" % len(existing_gis)

    if opts.possible_new_gb_out is None:
        option_parser.error("Need to specify --possible-new-gb-output")

    if opts.cached_ids:
        possible_gis = set([l.strip() for l in open(opts.cached_ids)])
    else:
        #ncbi_record_queries = ['16S','18S','small subunit','rrna[fkey]','ribosomal']
        ncbi_record_queries = ['16S AND tm7']
        # grab all the ids
        possible_gis = set([])
        for query in ncbi_record_queries:
            if opts.verbose:
                cur_size = len(possible_gis)
            possible_gis.update(esearch(query, retmax=10000000))

            if opts.verbose:
                print "Query %s added %d to set" % (query, len(possible_gis) - cur_size)

    # drop out any existing ids
    possible_gis = possible_gis - existing_gis

    if opts.verbose:
        print "Total number of GIs to query: %d" % len(possible_gis)
   
    chunk_count = 0
    total_bytes = 0
    if opts.use_gz:
        poss_output = open_gz(opts.possible_new_gb_out,'w')
    else:
        poss_output = open(opts.possible_new_gb_out,'w')
    
    collected = set([])

    retries = 0
    while possible_gis and retries < 100:
        try:
            for chunk in bulk_efetch(possible_gis):
                chunk_count += 1
                total_bytes += len(chunk)

                # Occasionally, and silently, NCBI corrupts records. 
                if '<html>' in chunk:
                    if verbose:
                        print "Erroneous record in chunk, disregarding full chunk"
                        continue

                # pullout the GIs
                records = [] 
                for l in chunk.splitlines():
                    if l.startswith('VERSION'):
                        records.append(l.split(':')[1])

                if opts.verbose:
                    print "%s - retry: %d, Chunk %d, covering %d records, writing %d bytes, %d written in total" % \
                        (time.strftime("%m-%d-%y %H:%M:%S"), retries, chunk_count, len(records), len(chunk), total_bytes)
                poss_output.write(chunk)
                collected.update(set(records))
        except Exception, e:
            retries += 1
            print "Caught exception: ", e
        possible_gis = possible_gis - collected
        collected = set([])
        
        possible_gis_at_retry = open('possible_retries_at_retry_%d.txt.gz' % retries, 'w')
        possible_gis_at_retry.write('\n'.join(possible_gis))
        possible_gis_at_retry.close()

if __name__ == '__main__':
    main()

