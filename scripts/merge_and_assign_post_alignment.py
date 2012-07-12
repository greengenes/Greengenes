#!/usr/bin/env python

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option
from greengenes.util import greengenes_open as open, NoSequenceError, \
        WorkflowLogger, generate_log_fp, log_f, GreengenesRecord
from greengenes.write import write_gg_record
from greengenes.parse import parse_column, parse_invariants
from cogent.parse.greengenes import MinimalGreengenesParser
from os import makedirs
from sys import stderr, stdout, argv
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

script_info={}
script_info['brief_description']="""Incorporate sequences into GG records computing metrics where appropriate"""
script_info['script_description']="""This script consumes the unaligned and NAST width aligned sequences from SSU-Align and merges the into the Greengenes records. Any previously unobserved Genbank accession will be assigned a new unique gg_id (and identical prokMSA_id if archaea/bacteria). Additionally, metrics that can be computed, such as the percent invariant score, are computed and merged in as well. The output are updated GG records and a 1-1 mapping of the new gg_ids to Genbank accessions. WARNING: This is a memory expensive operation as records are held in memory."""
script_info['script_usage']=[]
script_info['required_options'] = [\
        make_option('-r','--gg_records',type='str',
            help="Files containing Greengenes records"),
        make_option('-o','--output-dir',type='str',
            help="Output directory"),
        make_option('--unaligned',type='str',help='paths to unaligned sequences'),
        make_option('--aligned',type='str',
            help="Paths to aligned sequences"),
        make_option('-e','--existing',type='str',
            help="Path to a file containing previously observed genbank accessions"),
        make_option('--tag',type='str',help="Output tag, gg.YY.M"),
        make_option('--starting_gg_id',type='int',help="staring gg_id assign value"),
        make_option('--invariants',type='str',\
            help="Path to the invariant maps")]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    gg_records_fp = opts.gg_records
    output_dir = opts.output_dir
    verbose = opts.verbose
    existing_fp = opts.existing
    tag = opts.tag
    gg_id = opts.starting_gg_id


    invariants = parse_invariants(open(opts.invariants))

    makedirs(output_dir)
    logger = WorkflowLogger(generate_log_fp(output_dir), script_name=argv[0])

    # gg records are not going out as gzip as python's gzip is slow relative
    # to native linux gzip and doesn't compress as well out the door (latter 
    # probably fixable)
    output_gg_fp = os.path.join(output_dir, "%s.records.txt" % tag)
    output_map_fp = os.path.join(output_dir, "%s.mapping.txt.gz" % tag)
    output_gg_noggid_fp = os.path.join(output_dir, "%s.records.noggid.txt" \
                                                    % tag)
    
    existing_records = parse_column(open(existing_fp))
    
    records = dict([(r['ncbi_acc_w_ver'], GreengenesRecord(r)) \
                    for r in MinimalGreengenesParser(open(gg_records_fp))])
    
    for f in opts.aligned.split(','):
        logline = log_f("Parsing %s..." % f)
        logger.write(logline)
        if verbose:
            stdout.write(logline)

        domain = get_domain(f)

        for aln_id, aln_seq in MinimalFastaParser(open(f)):
            id_ = aln_id.split()[0] # strip of any comments
            record = records.get(id_, None)

            if record is None:
                logline = log_f("Aligned seq %s does not have a GG record" % id_)
                logger.write(logline)
                if verbose:
                    stdout.write(logline)
                continue

            if id_ in existing_records:
                logline = log_f("%s has previously been observed!" % id_)
                logger.write(logline)
                if verbose:
                    stdout.write(logline)
                continue

            if record['gg_id'] is not None:
                logline = log_f("%s already has gg_id %d!" %\
                                    (id_,record['gg_id']))
                logger.write(logline)
                if verbose:
                    stdout.write(logline)
                continue
        
            record['gg_id'] = gg_id
            if domain != 'eukarya':
                record['prokMSA_id'] = gg_id
            gg_id += 1

            inv_score = calc_invariant(seq, invariants)
            non_ACGT = calc_nonACGT(seq)

            record['perc_ident_to_invariant_core'] = inv_score
            record['non_ACGT_percent'] = non_ACGT
            record['aligned_seq'] = seq
            record['n_pos_aligned'] = len(seq) - seq.count('-')

    for f in opts.unaligned.split(','):
        logline = log_f("Parsing %s..." % f)
        logger.write(logline)
        if verbose:
            stdout.write(logline)

        domain = get_domain(f)

        for unaln_id, unaln_seq in MinimalFastaParser(open(f)):
            id_ = unaln_id.split()[0] # strip off any comments
            record = records.get(id_, None)

            if record is None:
                logline = log_f("Unaligned seq %s does not have a GG record" %\
                                 id_)
                logger.write(logline)
                if verbose:
                    stdout.write(logline)
                continue
    
            # a gg_id should be assigned while trolling the alignment seqs
            if record['gg_id'] is None:
                logline = log_f("%s should have a gg_id by now!" % (id_))
                logger.write(logline)
                if verbose:
                    stdout.write(logline)
                continue

            record['unaligned_seq'] = seq
            record['n_pos_unaligned'] = len(seq)
    
    logline = log_f("Beginning output...")
    logger.write(logline)
    if verbose:
        stdout.write(logline)

    output_map = open(output_map_fp,'w')
    output_gg = open(output_gg_fp,'w')
    output_gg_noggid = open(output_gg_noggid_fp, 'w')
    output_gg_broken = open(output_gg_broken_fp, 'w')

    for record in records.items():
        if record['gg_id'] is None:
            write_gg_record(output_gg_noggid, record)
        else:
            try:
                record.sanityCheck()
            except:
                write_gg_record(output_gg_broken, record)
            else:
                write_gg_record(output_gg, record)
                output_map.write("%s\t%s\n" % (record['gg_id'], 
                                               record['ncbi_acc_w_ver']))
    output_gg.close()


if __name__ == '__main__':
    main()
