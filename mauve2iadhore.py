#!/usr/bin/env python

from sys import stdout,stderr,exit
from optparse import OptionParser
from os.path import basename, exists, isdir
from itertools import combinations
from multiprocessing import cpu_count
from Bio import SeqIO
from bisect import bisect
import logging
import csv
import os

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

MIN_LENGTH_DEFAULT = 100
MAX_GAPSIZE_DEFAULT = 30
QVALUE_DEFAULT = 0.5

def readChrLengths(seqFiles, genomeNames):
    chr_locations = dict()

    for f in seqFiles:
        # assume that fasta file contains a single entry
        ident = basename(f).rsplit('.', 1)[0]

        LOG.info('reading %s with identifier %s' %(f, ident))
        cur_len = 0
        for record in SeqIO.parse(open(f), 'fasta'):
            if ident not in chr_locations:
                chr_locations[ident] = list() 
            chr_locations[ident].append((cur_len, record.id))
            cur_len += len(record)
    
    # put list into order indicated by mauve backbone filename
    chr_loc_list  = list()
    for ident in genomeNames:
        if ident not in chr_locations:
            LOG.error(('No corresponding fasta file found for identifier %s'
                + ' of backbone \'%s\'. Exiting.') %(ident,
                    '_'.join(genomeNames)))
            exit(1)

        chr_loc_list.append(tuple(chr_locations[ident]))

    return chr_loc_list 

def writeIadhoreFiles(genomes, gNames, orthologies, clusterType, gapSize, q,
        outDir):

    blast_out = open('%s/blast_table.txt' %outDir, 'w')
    for orth in orthologies:
        print >> blast_out, '%s:%s:%s\t%s:%s:%s' %(orth[0][0], orth[0][2],
                orth[0][3], orth[1][0], orth[1][2], orth[1][3])
    blast_out.flush()
    blast_out.close()

    configOut = open('%s/dataset.ini' %(outDir), 'w')
    
    for i in xrange(len(genomes)):
        gDir = '%s/%s' %(outDir, gNames[i])
        if not exists(gDir) or not isdir(gDir):
            os.mkdir(gDir)
        
        print >> configOut, 'genome= %s' %gNames[i]

        cur_out = prev_chr = None
        for gid, chr_name, start, end, strand in genomes[i]:
            if prev_chr != chr_name:
                print >> configOut, '%s %s/%s.lst' %(chr_name, gDir, chr_name)
                if cur_out != None:
                    cur_out.flush()
                    cur_out.close()
                cur_out = open('%s/%s.lst'%(gDir, chr_name), 'w')
                prev_chr = chr_name
            print >> cur_out, '%s:%s:%s%s' %(gid, start, end, strand)

        configOut.write('\n')
        if cur_out != None:
            cur_out.flush()
            cur_out.close()

    print >> configOut, 'blast_table= blast_table.txt'
    print >> configOut, 'output_path= output'
    print >> configOut, 'q_value= ' %q
    print >> configOut, 'tandem_gap=%s' %min(gapSize/2, 5)
    print >> configOut, 'gap_size= %s' %gapSize
    print >> configOut, 'cloud_gap_size= %s' %gapSize
    print >> configOut, 'cluster_gap= %s' %gapSize
    print >> configOut, 'cloud_cluster_gap= %s' %gapSize
    print >> configOut, 'cloud_filter_method= binomial_corr'
    print >> configOut, 'anchor_points= 5'
    print >> configOut, 'prob_cutoff= 0.001'
    print >> configOut, 'cluster_type= %s' %clusterType
    print >> configOut, 'number_of_threads= %s' %cpu_count()
    print >> configOut, 'alignment_method= gg2' 
    print >> configOut, 'visualizeGHM= false'
    print >> configOut, 'write_stats= true'

    configOut.close()


def parseMauveBackbone(chr_locations, mauveFile, minLength):

    genomes = [list() for _ in chr_locations]
    orthologies = list()

    isHeader = True
    outCount = 0

    for line in csv.reader(open(mauveFile), delimiter='\t'):
        if isHeader:
            isHeader = False
            continue

        outCount += 1

        cur_orth = list()
        for i in range(len(line)/2):
            if line[i*2] == '0':
                continue
            
            orient = int(line[i*2]) >= 0 and '+' or '-'
            start, end = abs(int(line[i*2]))-1, abs(int(line[i*2+1]))

            if end-start <= 0:
                LOG.warning(('sequence length of segment %s[%s:%s] is %s in' + \
                        ' line: \n\t\t%s') %(genomeNames[i], line[i*2], line[i*2+1],
                            end-start-1, '\t'.join(line)))

            if end-start-1 < minLength:
                continue

            x = bisect(chr_locations[i], (start, chr(255)))
            chr_start, chr_name = chr_locations[i][x-1]
            gid = '%s_%s' %(genomeNames[i], outCount)
            g1i = (gid, chr_name, start+1-chr_start, end-chr_start, orient)
            genomes[i].append(g1i)
            cur_orth.append(g1i)
        
        if len(cur_orth) > 1:
            orthologies.extend(combinations(cur_orth, 2))

    for i in xrange(len(genomes)):
        genomes[i].sort(key=lambda x: x[1:3])
    
    return genomes, orthologies

if __name__ == '__main__':

    usage = 'usage: %prog [options] <MAUVE BACKBONE> <FASTA FILE 1> ... <FASTA FILE N>'
    parser = OptionParser(usage=usage)

    parser.add_option('-l', '--minlength', dest='minLength',
            help='Miminimum length a segment must have to be included in ' \
                    + 'the output. [default=%default]',
                    type=int, default=MIN_LENGTH_DEFAULT, metavar='INT')

    parser.add_option('-g', '--max_gap_size', dest='gapSize', help='Maximum' + \
            ' distance between markers in a cloud or collinear cluster. ' + \
            '[default=%default]', type=int, default=MAX_GAPSIZE_DEFAULT,
            metavar='INT')

    parser.add_option('-c', '--cluster_type', dest='clusterType',
            default='hybrid', type='str', help='i-AdHoRe cluster type, ' + \
                    'must be any of (colinear|hybrid|cloud) [default: ' + \
                    '%default]')

    parser.add_option('-q', '--q_value', dest='q', help='Minimum r^2 ' + \
            'distance for markers in a cloud cluster. [default=%default]',
            type=float, default=QVALUE_DEFAULT, metavar='[0, 1]')

    parser.add_option('-o', '--outdir', dest='outDir', default='.', type=str, 
            help='Direction to which the output files will be written. ' + \
                    '[default=%default]')

    (options, args) = parser.parse_args()

    if len(args) < 3:
        parser.print_help()
        exit(1)

    mauveFile = args[0]
    seqFiles = args[1:]

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    
    cf = logging.FileHandler('%s.log' %(basename(args[0]).rsplit('.', 1)[0]), mode='w')
    cf.setLevel(logging.INFO)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t++ %(message)s'))

    LOG.addHandler(cf)
    LOG.addHandler(ch)

    #
    # main 
    #

    LOG.info(('start partitioning fasta sequences according to mauve ' \
            + 'backbone data using minimum segment length %s') %options.minLength)
    genomeNames = basename(mauveFile).rsplit('.',1)[0].split('_')
    chr_locations = readChrLengths(seqFiles, genomeNames)
    LOG.info('processing backbone file')


    genomes, orthologies = parseMauveBackbone(chr_locations, mauveFile,
            options.minLength)

    LOG.info('#bp covered by mauve markers of min length %s:' %(options.minLength))
    for i in xrange(len(genomes)):
        LOG.info('\t%s\t%s' %(genomeNames[i], sum(map(lambda x: x[3]-x[2]+1,
            genomes[i]))))

    LOG.info('constructing backbone file')
    writeIadhoreFiles(genomes, genomeNames, orthologies, options.clusterType,
            options.q, options.outDir)
    LOG.info('finished')

