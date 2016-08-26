#!/usr/bin/env python

from sys import stdout, stderr, exit
from optparse import OptionParser
from Bio import SeqIO

if __name__ == '__main__':

    usage = 'usage: %prog [options] <MIN LENGTH> <FASTA FILE>'
    parser = OptionParser(usage=usage)

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        exit(1)

    minLen = int(args[0])
    for rec in SeqIO.parse(args[1], 'fasta'):
        if len(rec) >= minLen:
            SeqIO.write(rec, stdout, 'fasta')
