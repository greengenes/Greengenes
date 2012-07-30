#!/usr/bin/env python

"""chop the genus name from the species names"""

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def chop_name(species_name):
    """takes a species name (s__foo) and guts the genus"""
    if not species_name.startswith('s__'):
        raise ValueError, "Does not appear like a species name"

    fullname = species_name.split('__',1)[1]
    splitnames = fullname.split()
    
    if 'sp.' in splitnames:
        return "s__"

    if len(splitnames) < 2:
        cleaned = species_name

    if len(splitnames) == 2:
        if splitnames[0].lower() == 'candidatus':
            cleaned = "s__"
        else:
            cleaned = "s__%s" % splitnames[1]

    if len(splitnames) == 3:
        if splitnames[0].lower() == 'candidatus':
            cleaned = "s__%s" % splitnames[2]
        else:
            cleaned = "s__%s" % splitnames[1]

    if len(splitnames) > 3:
        if splitnames[0].lower() == 'candidatus':
            cleaned = "s__%s" % ' '.join(splitnames[1:])
        else:
            cleaned = species_name
    return cleaned

if __name__ == '__main__':
    from sys import argv

    conlines = [l.strip().split('\t') for l in open(argv[1])]
    tax = [(a,b.split('; ')) for a,b in conlines]

    output = []
    for id_,(k,p,c,o,f,g,s) in tax:
        new_s = chop_name(s)
        output.append('\t'.join([id_, '; '.join([k,p,c,o,f,g,new_s])]))
    f = open(argv[2],'w')
    f.write('\n'.join(output))
    f.close()


