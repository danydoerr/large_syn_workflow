#!/usr/bin/env python

from sys import stdout, stderr, exit
from optparse import OptionParser
from Bio import SeqIO

if __name__ == '__main__':

    usage = 'usage: %prog [options] <FASTA FILE>'
    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--cumulative', dest='cumulative',
            help='Print sum length of fasta records', action='store_true',
            default=False) 

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    c = 0
    for rec in SeqIO.parse(args[0], 'fasta'):
        if options.cumulative:
            c += len(rec)
        else:
            print '%s\t%s' %(rec.id, len(rec))
    if options.cumulative:
        print c
