#!/usr/bin/env python
"""Application controller for bellerophon3"""

from cogent.app.parameters import ValuedParameter, FlagParameter, \
       MixedParameter
from cogent.app.util import CommandLineApplication, FilePath, system, \
       CommandLineAppResult, ResultPath, remove, ApplicationError, \
       get_tmp_filename
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.core.alignment import SequenceCollection
from cogent.app.blast import blastn
from math import ceil

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Daniel McDonald", "Thomas Huber"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class Bel3(CommandLineApplication):
    """FastTree application Controller"""

    _command = 'b3'
    _input_handler = '_input_as_multiline_string'
    _parameters = {'-i':ValuedParameter('-',Delimiter=' ',Name='i', 
                        IsPath=True),
                   '-o':ValuedParameter('-',Delimiter=' ',Name='o', 
                        IsPath=True),
                   '-t':ValuedParameter('-',Delimiter=' ',Name='t',
                        IsPath=True),
                   '-c':ValuedParameter('-',Delimiter=' ',Name='c'),
                   '-w':ValuedParameter('-',Delimiter=' ',Name='w'),
                   '-f':ValuedParameter('-',Delimiter=' ',Name='f')}

    def _get_result_paths(self, data):
        result = {}

        if self.Parameters['-o'].isOn():
            opath = self._absolute(str(self.Parameters['-o'].Value))
            result['B3out'] = ResultPath(Path=opath, IsWritten=True)
        return result

    def _input_as_seqs(self, data):
        lines = []
        for i,s in data.items():
            lines.append(">%s\n%s\n" % (i,s))
        return self._input_as_lines(lines)
    
    def _input_as_lines(self, data):
        if data:
            self.Parameters['-i'].on(super(Bel3, self)._input_as_lines(data))
        return ''

def fractionate_sequence(seq, n=4):
    """Equally partition up a sequence"""
    size = int(ceil(len(seq) / float(n)))
    return [seq[i*size:(i*size)+size] for i in range(n)]

def pick_parents(blastdb, ref_seqs,target_sequences,nhits=10,params=None):
    """Pick the most likely parents for a chimeric sequence

    blastdb : pre-formatted db of all parents
    ref_seqs : aligned parent sequences
    target_sequences : a dict of sequences
    nhits : number of "best" parents to take
    """
    for seq_id, seq in target_sequences:
        query_seqs = [("%d" % i, s) for i,s in \
                                    enumerate(fractionate_sequence(seq))]
        result = blastn(query_seqs, blastdb)
        ref_db = {}
        for q, best_hits in result.bestHitsByQuery(n=nhits):
            for rec in best_hits:
                id_ = rec['SUBJECT ID']
                ref_db[rec['SUBJECT ID']] = ref_seqs[id_]
        
        yield (ref_db,seq_id, seq)

def check_chimera(refseqs, target_id, target_seq):
    """Check if target is a chimera

    refseqs : something like a dict {id:seq}
    target_id : target sequence id, string
    target_seq : the actual target sequence

    expects the refseqs and target seq to both be aligned against same ref
    """
    assert target_id not in refseqs

    inputseqs = refseqs.copy()
    inputseqs[target_id] = target_seq

    params = {'-o':get_tmp_filename(), '-w':400, '-t':target_id, '-f':'full',
              '-c':'Huber-Hugenholtz'}
    app = Bel3(InputHandler='_input_as_seqs', params=params, HALT_EXEC=False)
    res = app(inputseqs)

    how_chimeric = parse_bel3_result(res['B3out'])

    res.cleanUp()

    return how_chimeric

def parse_bel3_result(lines):
    """parses bel3 full output

    returns None if it does not appear chimeric at all, otherwise returns
    the 'preference score'
    """
    has_result = False
    parent_a = None
    parent_b = None

    for l in lines:
        if l.startswith('divergence ratio:'):
            has_result = float(l.strip().split()[-1])
            continue
        if l.startswith('parent 1 '):
            parent_a = l.strip().split('=>')[1]
            continue
        if l.startswith('parent 2 '):
            parent_b = l.strip().split('=>')[1]
            continue
        
        if l.startswith('>') and 'appears to be clean' in l:
            has_result = None
            break
        if l.startswith('Local surrounding'):
            # nothing more useful here
            break

    if has_result is False:
        raise ApplicationError, "Cannot interpret result"

    return (has_result, parent_a, parent_b)
