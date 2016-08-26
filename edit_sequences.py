#!/usr/bin/env python

from sys import stdout, stderr, exit
from optparse import OptionParser
from Bio import SeqIO, SeqRecord
from os.path import basename
import re

PAT_RENAME_CMD = re.compile('^\s*RENAME\s+(\S+)\s+TO\s+(\S+)\s*$')
PAT_REMOVE_CMD = re.compile('^\s*REMOVE\s+(\S+)\s*$')
PAT_JOIN_CMD = re.compile('^\s*JOIN\s+(\S+)\s+(t|h)\s+(\S+)\s+(t|h)\s+INTO\s+(\S+)\s*$')

if __name__ == '__main__':

    usage = 'usage: %prog [options] <FASTA FILE> <INSTRUCTION FILE>'
    parser = OptionParser(usage=usage)

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        exit(1)

    recs = dict()
    for rec in SeqIO.parse(args[0], 'fasta'):
        recs[rec.id] = rec


    c = 0
    for line in open(args[1]):
        c += 1

        if not line.strip():
            continue

        #
        # RENAME COMMAND
        #
        m = PAT_RENAME_CMD.match(line)
        if m:
            old_name, new_name = m.groups()
            
            if not recs.has_key(old_name):
                print >> stderr, ('ERROR in line %s of %s: No sequence with ' + \
                        'ID %s in current set of sequence records found. ' + \
                        'Exiting') %(c, basename(args[1]), old_name)
                exit(1)

            if recs.has_key(new_name):
                print >> stderr, ('ERROR in line %s of %s: New ID \'%s\' ' + \
                        'for is already taken by another sequence. ' + \
                        'Exiting') %(c, basename(args[1]), new_name)
                exit(1)

            rec = recs.pop(old_name)
            rec.id = new_name
            recs[new_name] = rec
            continue 
            
        #
        # REMOVE COMMAND
        #
        m = PAT_REMOVE_CMD.match(line)
        if m:
            name = m.group(1)
            
            if not recs.has_key(name):
                print >> stderr, ('ERROR in line %s of %s: No sequence with ' + \
                        'ID %s in current set of sequence records found. ' + \
                        'Exiting') %(c, basename(args[1]), name)
                exit(1)

            rec = recs.pop(name)
            continue 

        #
        # JOIN COMMAND
        #
        m = PAT_JOIN_CMD.match(line)
        if m:
            seq1, o1, seq2, o2, new_name = m.groups()
            
            if not recs.has_key(seq1) or not recs.has_key(seq2):
                print >> stderr, ('ERROR in line %s of %s: %s and/or %s are not ' + \
                        'associated with current set of seqeunces. Exiting') %(c,
                                basename(args[1]), seq1, seq2)
                exit(1)

            if new_name not in [seq1, seq2] and recs.has_key(new_name):
                print >> stderr, ('ERROR in line %s of %s: ID \'%s\' for new ' + \
                        '(joined) sequence record is already used by another ' + \
                        'sequence. Exiting') %(c, basename(args[1]), new_name)
                exit(1)

            rec1 = recs.pop(seq1)
            rec2 = recs.pop(seq2)

            new_rec = None

            if o1 == 't' and o2 == 'h':
                new_rec = SeqRecord.SeqRecord(rec1.seq + rec2.seq, id=new_name,
                        description='')
            elif o1 == 't' and o2 == 't':
                new_rec = SeqRecord.SeqRecord(rec1.seq +
                        rec2.seq.reverse_complement(), id=new_name, description='')
            elif o1 == 'h' and o2 == 'h':
                new_rec = SeqRecord.SeqRecord(rec2.seq.reverse_complement() +
                        rec1.seq, id=new_name, description='')
            elif o1 == 'h' and o2 == 't':
                new_rec = SeqRecord.SeqRecord(rec2.seq + rec1.seq, id=new_name,
                        description='')

            recs[new_name] = new_rec
        
        else:
            print >> stderr, ('ERROR in line %s of %s: Unable to parse command ' + \
                    '\'%s\'') %(c, basename(args[1]), line.strip())
            continue
        
    for k in sorted(recs.keys()):
        SeqIO.write(recs[k], stdout, 'fasta')

