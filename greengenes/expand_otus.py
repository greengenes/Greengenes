#!/usr/bin/env python

from greengenes.parse import parse_otus

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def expand(tax,otus):
    """Expands the taxonomy for representative seqs over the clusters"""
    result = []
    for otu_id,cluster in otus:
        rep_id = cluster[0]
        rep_tax = tax[rep_id]
        result.extend([(seqid, rep_tax) for seqid in cluster])
    return result
