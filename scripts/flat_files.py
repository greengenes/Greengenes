#!/usr/bin/env python

from cogent.parse.genbank import MinimalGenbankParser, PartialRecordError
from cogent.util.misc import parse_command_line_parameters
from optparse import make_option
from greengenes.flat_files import get_sequence, get_genbank_summary, \
        get_accession
from greengenes.write import write_sequence, write_gg_record, \
        write_obs_record
from greengenes.parse import parse_column
from greengenes.util import greengenes_open as open, NoSequenceError, \
        WorkflowLogger, generate_log_fp, log_f
from sys import stdout, stderr, argv
from os import makedirs
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

script_info={}
script_info['brief_description']="""Convert possible new Greengenes records to easily consumable files"""
script_info['script_description']="""This script consumes Genbank records and dumps out summary information in a tab-delimited format and sequences with accessions for IDs for all parsable records. Failures are recorded. Records are checked by Genbank accession to determine if they've already been indexed by Greengenes."""
script_info['script_usage']=[]
script_info['required_options'] = [\
        make_option('-i','--input-gbs',type='str',
            help="Files containing Genbank records"),
        make_option('-o','--output-dir',type='str',
            help="Output directory"),
        make_option('-t','--tag',type='str',help='Output filename prefix'),
        make_option('-e','--existing',type='str',
            help="File containing previously observed accessions in the first column")]
script_info['optional_options'] = [\
        make_option('--max-failures',type='int', default=10000,
            help='Maximum parse errors per genbank file')]      
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    input_gbs = opts.input_gbs.split(',')
    output_dir = opts.output_dir
    verbose = opts.verbose
    tag = opts.tag
    existing_fp = opts.existing
    max_failures = opts.max_failures
    
    makedirs(output_dir)
    logger = WorkflowLogger(generate_log_fp(output_dir), script_name=argv[0])

    observed_records = parse_column(open(existing_fp))

    sequences_fp = os.path.join(output_dir, '%s_sequences.fasta.gz' % tag)
    gg_records_fp = os.path.join(output_dir, '%s_ggrecords.txt.gz' % tag)
    obs_records_fp = os.path.join(output_dir, '%s_obsrecords.txt.gz' % tag)
    
    sequences = open(sequences_fp,'w')
    gg_records = open(gg_records_fp, 'w')
    obs_records = open(obs_records_fp, 'w')
    
    seen = set([])
    for gb_fp in input_gbs:
        logline = log_f("Start parsing of %s..." % gb_fp)
        logger.write(logline)

        if verbose:
            stdout.write(logline)

        records = MinimalGenbankParser(open(gb_fp))
        
        failure_count = 0
        alpha = set(['A','T','G','C',
                     'a','t','g','c',
                     'N','n',
                     'R','Y','S','M',
                     'r','y','s','m',
                     'K','k','W','w',
                     'V','v','H','h','B','b','D','d'])

        while True and (failure_count < max_failures):
            # gracefully handle parser errors to a limit
            try:
                next_record = records.next()
            except PartialRecordError, e:
                failure_count += 1
                continue
            except StopIteration:
                break
            except Exception, e:
                logline = log_f("Caught: %s, previous accession: %s" % (e, accession))
                logger.write(logline)
                if verbose:
                    stdout.write(logline)
                failure_count += 1

            # accession is str including version
            try:
                accession = get_accession(next_record)
            except:
                failure_count += 1
                continue
            if accession in observed_records:
                continue
            if accession in seen:
                continue
                # accession added at the end of this while loop

            # sequence is just a str of sequence
            try:
                sequence = get_sequence(next_record)
            except NoSequenceError:
                # this isn't a failure, so no point in continuing but record
                # the accession so it isn't hit again
                write_obs_record(obs_records, accession)
                continue
            except:
                failure_count += 1
                continue

            # verify the sequence is DNA. NCBI silently corrupts records 
            # every once and a while leading to crap in the sequence. 
            # Thanks NCBI.
            seq_chars = set(sequence)
            if not seq_chars.issubset(alpha):
                logline = log_f("Corrupt sequence, accession: %s" % (accession))
                logger.write(logline)
                if verbose:
                    stdout.write(logline)
                failure_count += 1
                continue

            # gg_record contains gb summary data  
            try:
                gg_record = get_genbank_summary(next_record)
            except KeyError, e:
                failure_count += 1
                continue

            seen.add(accession)
            write_sequence(sequences, accession, sequence)
            write_gg_record(gg_records, gg_record)
            write_obs_record(obs_records, accession)
            
        if failure_count >= max_failures:
            logline = log_f("MAX FAILURES OF %d REACHED IN %s" % (max_failures, \
                                                                  gb_fp))
            logger.write(logline)
            stderr.write(logline)
        else:
            logline = log_f("Parsed %s, %d failures observed." % (gb_fp, \
                                                                  failure_count))
            logger.write(logline)

            if verbose:
                stdout.write(logline)

    sequences.close()
    gg_records.close()
    obs_records.close()

if __name__ == '__main__':
    main()
