#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def parse_existing_records(open_file):
    """Writes out the header line for genbank summaries"""
    return set([l.strip().split('\t')[0] for l in open_file \
                                         if not l.startswith('#')])

