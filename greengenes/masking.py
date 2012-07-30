#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def inflate_by_mask(seq, mask):
    """Inflate seq to the mask

    complain if the number of positions in seq is != to sum(mask)
    """
    s = []
    idx = 0 
    for m in mask:
        if m:
            s.append(seq[idx])
            idx += 1
        else:
            s.append('-')

    if len(seq) != sum(mask):
        raise ValueError, "seq does fit mask!"

    return ''.join(s)

def img_best_16s_per_genome(masked, unmasked):
    """Keeps the longest masked sequence per genome
    
    expects IDs to be "accession|genome_id"
    """
    to_keep = {}
    len_genome = {}

    for seqid, seq in masked.items():
        if '|' in seqid:
            accession, genome_id = seqid.split('|')
        else:
            genome_id = seqid
            accession = None

        length = len(seq.replace('-','').replace('.',''))

        if genome_id not in len_genome:
            len_genome[genome_id] =[]

        len_genome[genome_id].append((length, accession))

    for genome_id, lengths in len_genome.items():
        masked_length, accession = sorted(lengths)[-1] # take the longest one
        
        if accession:
            seq_id = '|'.join([accession,genome_id])
        else:
            seq_id = genome_id

        unmasked_length = len(unmasked[seq_id].replace('-','').replace('.',''))

        seq_info = '\t'.join(['masked_length=%d' % masked_length,
                              'unmasked_length=%d' % unmasked_length])

        seq_header = '\t'.join([seq_id, seq_info])
        
        to_keep[seq_header] = masked[seq_id]

    return to_keep


if __name__ == '__main__':
    from cogent.parse.fasta import MinimalFastaParser
    from sys import argv

    masked = dict(MinimalFastaParser(open(argv[1])))
    unmasked = dict(MinimalFastaParser(open(argv[2])))
    mask = [int(c) for c in open(argv[3]).read().strip()]
    best = img_best_16s_per_genome(masked,unmasked)

    f = open(argv[4],'w')
    for seqid, seq in best.items():
        inflated = inflate_by_mask(seq, mask)
        f.write('>%s\n%s\n' % (seqid, inflated))

    f.close()

