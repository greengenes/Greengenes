#!/usr/bin/env python

def calc_invariant(seq, invariants):
    """Return the percent of bases that don't agree with invariants"""
    hits = [] 
    for id_, inv_map, inv_length in invariants:
        hit = 0.0
        for a,b in zip(inv_map, seq):
            if a == 'N':
                continue
            if a == b:
                hit += 1
        hits.append(hit / inv_length)
    return max(hits)

def calc_nonACGT(seq):
    """Calculate the percentage of non-ACGT characters"""
    aln_count = 0
    non_acgt = 0
    for c in seq:
        if c == '-':
            continue
        aln_count += 1
        if c not in ['A','G','C','T']:
            non_acgt += 1
    return float(non_acgt) / aln_count


