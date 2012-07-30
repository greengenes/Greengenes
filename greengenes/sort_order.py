#!/usr/bin/env python

from cogent.parse.tree import DndParser

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def parse_dec(lines):
    """Returns {a:(decision, has_sp)}

    expects: id\tdecision\torganism
    """
    res = {}
    for line in lines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        
        if len(fields) == 2:
            id_, dec = fields
            organism = ""
        else:
            id_, dec, organism = fields

        if dec != 'named_isolate':
            res[id_] = ('clone',False)
        else:
            if 'sp.' in organism:
                res[id_] = ('named_isolate',False)
            else:
                res[id_] = ('named_isolate',True)
    return res

def parse_lengths(lines):
    """Returns {id:length}"""
    res = {}
    for l in lines:
        if l.startswith('#'):
            continue
        id_, length = l.strip().split('\t',1)
        res[id_] = int(length)
    return res

def merge_dec_length(in_pref, lengths, decinfo):
    """Returns {id_:(length, decision, in_pref, has_sp)}"""
    res = {}
    for k, (dec, has_sp) in decinfo.items():
        if k not in lengths:
            continue

        length = lengths[k]

        if k in in_pref:
            res[k] = (length, dec, True, has_sp)
        else:
            res[k] = (length, dec, False, has_sp)
    return res

def sort_order(records):
    """returns the sort order by id"""
    tree = DndParser("(((nosp,sp)named,notnamed)inpref,\
                       ((nosp,sp)named,notnamed)outpref);")
    for n in tree.tips():
        n.LengthsAndIds = []
    lookup = {}
    lookup[('named_isolate',True,True)] = \
            tree.Children[0].Children[0].Children[0]
    lookup[('named_isolate',True,False)] = \
            tree.Children[0].Children[0].Children[1]
    lookup[('clone',True,False)] = \
            tree.Children[0].Children[1]
    lookup[('named_isolate',False,True)] = \
            tree.Children[1].Children[0].Children[0]
    lookup[('named_isolate',False,False)] = \
            tree.Children[1].Children[0].Children[1]
    lookup[('clone',False,False)] = \
            tree.Children[1].Children[1]
                       
    for k,v in records.items():
        to_lookup = tuple(v[1:])
        lookup[to_lookup].LengthsAndIds.append((v[0],k))

    order = []
    # tips go left->right
    for n in tree.tips():
        order.extend([i for l,i in sorted(n.LengthsAndIds)[::-1]])

    return order

if __name__ == '__main__':
    from sys import argv
    dec = parse_dec(open(argv[1]))
    lengths = parse_lengths(open(argv[2]))
    preference = set([l.strip() for l in open(argv[3])])
    merged = merge_dec_length(preference, lengths,dec)

    order = sort_order(merged)
    f = open(argv[4],'w')
    f.write('\n'.join(order))
    f.close()
