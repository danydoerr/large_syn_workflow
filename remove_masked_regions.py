#!/usr/bin/env python

from sys import stdout, stderr, exit
from optparse import OptionParser
from Bio import SeqIO, Seq, Alphabet
from cStringIO import StringIO

MINLEN_DEFAULT=10

if __name__ == '__main__':

    usage = 'usage: %prog [options] <FASTA FILE>'
    parser = OptionParser(usage=usage)
    parser.add_option('-l', '--min_len', dest='minLength', type=int,
            default=MINLEN_DEFAULT, 
            help='Minimum number of consecutive X/N characters that are ' + \
                    'considered for removal. [default=%default]')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    for rec in SeqIO.parse(args[0], 'fasta'):
        p = 0
        c_count = 0
        new_seq = StringIO()
        for c in rec.seq:
            if c in {'X', 'x', 'N', 'n'}:
                c_count += 1
            elif c_count >= options.minLength:
                p -= c_count
                new_seq.seek(p)
                c_count = 0
            new_seq.write(c)
            p += 1

        rec.seq = Seq.Seq(new_seq.getvalue(), Alphabet.generic_dna)
        SeqIO.write(rec, stdout, 'fasta')
        

