#!/usr/bin/env python

def get_overlap(b3_chimeras, cs_chimeras):
    """Return the set intersection of the ids"""
    b3_ids = set([id_ for id_,a,b,c in b3_chimeras])
    cs_ids = set([id_ for id_,a,b in cs_chimeras])

    return b3_ids.intersection(cs_ids)

def determine_taxon_conflict(taxmap, parent_a, parent_b, rank=2):
    """Determine if there is a taxon conflict at a particular level

    taxmap : {id_:[taxon,name,...]}
    parent_a : id_ for parent_a
    parent_b : id_ for parent_b
    rank : level in the taxnomy string. rank 2 is class

    Returns true if there is a conflict
    """
    tax_a = taxmap[parent_a]
    tax_b = taxmap[parent_b]

    if tax_a[2] != tax_b[2]:
        return True
    else:
        return False
