#!/usr/bin/env python

from cogent.parse.fasta import MinimalFastaParser
from greengenes.util import GreengenesRecord

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def parse_b3_chimeras(lines):
    """Returns [(id_, score, parent_a, parent_b)]"""
    print "FORCING  > 1.2 divergence!"
    res = []
    for line in lines:
        fields = line.strip().split('\t')
        if len(fields) != 4:
            continue
        if float(fields[1]) < 1.2:
            continue
        res.append((fields[0], float(fields[1]), fields[2], fields[3]))
    return res

def parse_cs_chimeras(lines):
    """Returns [(id_, parent_a, parent_b)]"""
    res = []
    for line in lines:
        fields = line.strip().split('\t')
        if len(fields) != 3:
            continue
        res.append((fields[0], fields[1], fields[2]))
    return res

def parse_otus(lines):
    """parses an otu mapping file

    returns [(otu_id, [ids])]
    """
    result = []
    for l in lines:
        fields = l.strip().split('\t')
        result.append((fields[0], fields[1:]))
    return result

def parse_invariants(lines):
    """Parse the invariants

    Returns [(invariant_id, invariant_map, invariant_length)]
    """
    invs = list(MinimalFastaParser(lines))
    invs_lengths = [len(s) - s.count('N') for i,s in invs]

    result = []
    for (inv_id, inv_map),inv_len in zip(invs, invs_lengths):
        result.append((inv_id, inv_map, inv_len))

    return result

def parse_column(open_file):
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

