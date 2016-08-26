#!/usr/bin/env python

from sys import stdout, stderr, exit, maxint
from optparse import OptionParser
from itertools import product, combinations, izip
from os.path import basename, dirname, join
from random import shuffle
import logging
import csv
import re

MARKER_PAT = re.compile('^([^:]+):(\d+):(\d+)(\+$|-$|$)')

FONTSIZE_VIS = 20

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def readIadhoreConfig(data):

    res = list()
    for line in data:
        if not line.strip():
            continue
        if line.find('=') >= 0:
            k, v = line.split('=', 1)
            res.append((k.strip(), v.strip()))
        elif line.find(' ') >= 0 and len(res):
            if len(res[-1]) < 3:
                res[-1] = res[-1] + (list(), )
            k, v = line.split(' ', 1)
            res[-1][2].append((k.strip(), v.strip()))
        else:
            print >> stderr, ('Unable to parse line \"%s\" in i-AdHoRe ' + \
                    'config file. Exiting') %line.strip()
            exit(1)
    return res


def readBlastTable(data, gene2genome):
    res = dict()
    for g1i, g2k in csv.reader(data, delimiter='\t'):
        G1 = gene2genome[g1i]
        G2 = gene2genome[g2k]
        if not res.has_key((G1, G2)):
            res[(G1, G2)] = dict()
            res[(G2, G1)] = dict()
        if not res[(G1, G2)].has_key(g1i):
            res[(G1, G2)][g1i] = set()
        res[(G1, G2)][g1i].add(g2k)
        if not res[(G2, G1)].has_key(g2k):
            res[(G2, G1)][g2k] = set()
        res[(G2, G1)][g2k].add(g1i)
    return res


def readGenomes(iadhoreConfig, configPath):
    genomes = dict()
    gene2genome = dict()
    for x in iadhoreConfig:
        if len(x) != 3 or x[0] != 'genome':
            continue
        genomes[x[1]] = list() 
        for _, fPath in x[2]:
            f = join(configPath, fPath)
            for line in open(f):
                gx = line.strip()
                if gx[-1] not in ('+', '-'):
                    LOG.fatal(('Gene %s in file %s has unknown ' + \
                            'orientation.  Exiting') %(gx, f))
                    exit(1)
                genomes[x[1]].append(gx[:-1])
                gene2genome[gx[:-1]] = x[1]
    return genomes, gene2genome


def readSegments(data):

    isHeader = True
    
    res = dict()
    c = 0
    for line in csv.reader(data, delimiter='\t'):
        c += 1
        if isHeader:
            isHeader = False
            continue
        
        mid = int(line[1])
        if not res.has_key(mid):
            res[mid] = list()
        res[mid].append((line[2], line[4], line[5]))

    return res


def weightedScore(segments, genomes, gpos, blastMap):
    res = list()
    for mid in sorted(segments.keys()):
        c = 0.0
        s = 0
        for i in xrange(len(segments[mid])):
            G1, g1i, g1j = segments[mid][i]
            start1, end1 = sorted((gpos[G1][g1i], gpos[G1][g1j]))
            s += end1-start1+1
            for x in xrange(start1, end1+1):
                g1x = genomes[G1][x]
                hasHit = True
                j = 0
                while hasHit and j < len(segments[mid]):
                    j += 1
                    if i == j-1:
                        continue
                    G2, g2k, g2l = segments[mid][j-1]
                    hasHit = blastMap[(G1, G2)].has_key(g1x)
                    if not hasHit:
                        break
                    start2, end2 = sorted((gpos[G2][g2k], gpos[G2][g2l]))
                    hasHit_i = False
                    for g2y in blastMap[(G1, G2)][g1x]:
                        y = gpos[G2][g2y]
                        hasHit_i = y >= start2 and y <= end2
                        if hasHit_i:
                            break
                    hasHit = hasHit_i
                if hasHit:
                    c += 1
        res.append((mid, c/((len(segments[mid])-1) * s)))
    return res


def relaxedScore(segments, genomes, gpos, blastMap):
    res = list()
    for mid in sorted(segments.keys()):
        c = 0.0
        s = 0
        for i in xrange(len(segments[mid])):
            G1, g1i, g1j = segments[mid][i]
            start1, end1 = sorted((gpos[G1][g1i], gpos[G1][g1j]))
            s += end1-start1+1
            for x in xrange(start1, end1+1):
                g1x = genomes[G1][x]
                hasHit = False
                j = 0
                while not hasHit and j < len(segments[mid]):
                    j += 1
                    if i == j-1:
                        continue
                    G2, g2k, g2l = segments[mid][j-1]
                    if not blastMap[(G1, G2)].has_key(g1x):
                        break
                    start2, end2 = sorted((gpos[G2][g2k], gpos[G2][g2l]))
                    for g2y in blastMap[(G1, G2)][g1x]:
                        y = gpos[G2][g2y]
                        hasHit = y >= start2 and y <= end2
                        if hasHit:
                            break
                if hasHit:
                    c += 1
        res.append((mid, c/s))
    return res

def showHistogram(scores, stype, fileName):

    try: 
        import matplotlib.pyplot as plt
    except:
        LOG.fatal('Unable to import matplotlib.pyplot. Not installed? Exiting.')
        exit(1)

    # the histogram of the data
    n, bins, patches = plt.hist(map(lambda x: x[1], scores), 50,
            facecolor='#b30000', edgecolor='none', alpha=0.75)
    plt.xlabel(stype, fontsize=FONTSIZE_VIS)
    plt.ylabel('count', fontsize=FONTSIZE_VIS)
    #plt.axis([40, 160, 0, 0.03])
    #plt.grid(True)
    plt.savefig(fileName, formtat='eps', transparent=True)
    plt.show()

if __name__ == '__main__':

    usage = 'usage: %prog [options] <I-ADHORE CONFIGURAION FILE> <SEGMENTS FILE>'
    parser = OptionParser(usage=usage)
    parser.add_option('-v', '--visualize', dest='visual', default=False,
            action='store_true', help='Show histogram of scores. IMPORTANT: ' + \
                    'requires matplotlib.pyplot library! [default: %default]')
    parser.add_option('-t', '--type', dest='type', default='relaxed', type=str,
            help='Scoring type. [default: %default]', metavar=
            '(relaxed|weighted)')
    parser.add_option('-f', '--figure_name', dest='figName', default='(relaxed' + \
            '|weighted)_scores.eps', type=str, help='Name of output file of ' + \
            'histogram. Only applicable in combination with -v. [default:' + \
            '%default]')
    figNameOpt = parser.option_list[-1]

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        exit(1)


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

    if figNameOpt.default == options.figName:
        options.figName = '%s_scores.eps' %options.type

    iadhoreConfig = readIadhoreConfig(open(args[0]))

    iadhoreCMap = dict(x for x in iadhoreConfig if len(x) == 2)
    if iadhoreCMap.has_key('table_type') and \
            iadhoreCMap['table_type'] == 'family':
        LOG.fatal(('Unable to parse file %s with family assignments: not ' + \
            'implemented. Exiting.') %iadhoreCMap['blast_table'])
        exit(1)

    genomes, gene2genome = readGenomes(iadhoreConfig, dirname(args[0]))
    gpos = dict(map(lambda x: (x[0], dict(izip(x[1], xrange(len(x[1]))))),
        genomes.items()))
    segments = readSegments(open(args[1]))
    blastMap = readBlastTable(open(join(dirname(args[0]),
        iadhoreCMap['blast_table'])), gene2genome)


    scores = None
    if options.type == 'relaxed':
        scores = relaxedScore(segments, genomes, gpos, blastMap)
    elif options.type == 'weighted':
        scores = weightedScore(segments, genomes, gpos, blastMap)
    else:
        LOG.fatal('Scoring type must be either relaxed or weighted. Exiting')
        exit(1)

    if options.visual:
        showHistogram(scores, options.type + ' synteny score', options.figName)
    else:
        for s in scores:
            print >> stdout, '%s\t%s' %s

    LOG.info('finished')
