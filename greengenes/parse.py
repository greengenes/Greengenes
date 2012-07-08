#!/usr/bin/env python

from greengenes.util import GreengenesRecord

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def parse_existing_records(open_file):
    """Load a column of data"""
    return set([l.strip().split('\t')[0] for l in open_file \
                                         if not l.startswith('#')])

def parse_gg_summary_flat(open_file):
    """Parse a flat greengenes summary file from flat_files"""
    header_line = open_file.readline()
    if not header_line.startswith('#'):
        raise ValueError, "Missing the header!"

    header = header_line[1:].strip().split('\t')

    print "WARNING: NOT SETTING TYPES CURRENTLY"
    for line in open_file:
        record = GreengenesRecord()

        for key,value in zip(header, line.strip().split('\t')):
            record[key] = value

        #record.setTypes()
        yield record

