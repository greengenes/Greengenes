#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class ParseError(Exception):
    pass

def check_parse(line):
    """Make sure the line can parse into (id_, [some,thing])"""
    try:
        id_, initial = line.strip().split('\t')
    except ValueError:
        raise ParseError, "Unable to split in tab"

    parsed = initial.split("; ")

    if not len(parsed):
        raise ParseError, "Line appears to not have a taxonomy"

    return (id_, parsed)

def check_n_levels(parsed, n):
    """Make sure there are n levels in parsed"""
    if len(parsed) == n:
        return True
    else:
        return False

def check_gap(parsed):
    """Check that no gaps exist in parsed"""
    reduced = [p.split('__')[-1] for p in parsed]

    end_found = False
    for n in reduced:
        if n and end_found:
            return False

        if not n:
            end_found = True
    
    return True

def check_prefixes(parsed, expected_prefixes):
    """Make sure each rank has the expected prefix"""
    for name, prefix in zip(parsed, expected_prefixes):
        obs, level_name = name.split('__',1)
        if obs != prefix:
            return False

    return True

if __name__ == '__main__':
    from sys import argv, exit, stderr
  
    nlevels = int(argv[2])
    prefixes = argv[3].split(',')

    counter = 1
    for l in open(argv[1]):
        try:
            id_, parsed = check_parse(l)
        except ParseError, e:
            stderr.write("Error on line: %d\n" % counter)
            raise e
            exit(1)

        if not check_prefixes(parsed, prefixes):
            print "Line %d appears to have incorrect prefixes" % counter
            stderr.write(l)
        if not check_n_levels(parsed, nlevels):
            print "Line %d appears to have the wrong number of levels" %counter
            stderr.write(l)
        if not check_gap(parsed):
            print "Line %d appears to have a gap in the taxonomy" % counter
            stderr.write(l)
        
        counter += 1
