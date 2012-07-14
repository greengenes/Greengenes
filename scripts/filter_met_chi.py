#!/usr/bin/env python

"""filter by invariant, non_acgt and chimera"""

from sys import argv
from cogent.parse.fasta import MinimalFastaParser

if __name__ == '__main__':

    metric_lines = [l.strip().split('\t') for l in open(argv[1]) if not l.startswith('#')]
    chi_lines = [l.strip().split('\t') for l in open(argv[2]) if not l.startswith('#')]

    chi = set([r[0] for r in chi_lines])
    metrics = dict([(id_, (inv, nonacgt)) \
                        for id_,inv,nonacgt in metric_lines])

    seqs = MinimalFastaParser(open(argv[3]))
    output_seqs = open(argv[4], 'w')
    output_reasons = open(argv[4] + '.filtered', 'w')
    output_reasons.write("#ncbi_acc_w_ver\treason\n")

    for i,s in seqs:
        if i in chi:
            output_reasons.write("%s\tchimeric\n" % i)
            continue

        inv, nonacgt = metrics[i]

        if inv < 0.9:
            output_reasons.write("%s\tinvariants < 90%\n")
            continue
        if nonacgt > 0.01:
            output_reasons.write("%s\tnon ACGT > 1%\n")
            continue
        
        output_seqs.write(">%s\n%s\n" % (i,s))
