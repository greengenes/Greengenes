#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def write_sequence(fp, id_, seq):
    """Write a sequence to a file like object"""
    fp.write(">%s\n%s\n" % (id_, seq))

def write_gb_summary(fp, id_, vals):
    """Writes a gb summary line"""
    fp.write("%s\t%s\n" % (id_, '\t'.join(map(str, vals))))

def write_obs_record(fp, id_):
    """Write an id to a file"""
    fp.write("%s\n" % id_)

