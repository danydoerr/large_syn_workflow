#!/usr/bin/env python

from sys import stdout, stderr, exit, maxint
from optparse import OptionParser
from itertools import product, combinations, izip
from os.path import basename, dirname, join, isfile
from random import shuffle
import logging
import csv
import re

MARKER_PAT = re.compile('^([^:]+):(\d+):(\d+)(\+$|-$|$)')

FONTSIZE_VIS = 20

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

def readMarkerSequences(data):

    res = list()
    
    isHeader = True
    for _, segid, gene, _, _, in csv.reader(data, delimiter='\t'):
        if isHeader:
            isHeader = False
            continue
        if int(segid) > len(res):
            res.append(set())
        res[-1].add(gene)

    return res

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
        if isHeader:
            isHeader = False
            continue
        
        mid = int(line[1])
        if not res.has_key(mid):
            res[mid] = list()
        res[mid].append((c, line[2], line[4], line[5]))
        c += 1

    return res


def weightedScore(segments, marker_seqs, blastMap, onlyAll=False):
    res = list()

    if onlyAll:
        n = len(set(reduce(lambda x,y: x+y, blastMap.keys())))

    for mid in sorted(segments.keys()):
        if onlyAll and len(set(x for _, x, _, _ in segments[mid])) < n:
            continue
        c = 0.0
        s = 0
        for i in xrange(len(segments[mid])):
            sid1, G1, _, _ = segments[mid][i]
            s += len(marker_seqs[sid1])
            for g1x in marker_seqs[sid1]:
                hasHit = True
                j = 0
                while hasHit and j < len(segments[mid]):
                    j += 1
                    if i == j-1:
                        continue
                    sid2, G2, _, _ = segments[mid][j-1]
                    hasHit = blastMap.has_key((G1, G2)) and \
                            blastMap[(G1, G2)].has_key(g1x)
                    if not hasHit:
                        break
                    hasHit_i = False
                    for g2y in blastMap[(G1, G2)][g1x]:
                        hasHit_i = g2y in marker_seqs[sid2]
                        if hasHit_i:
                            break
                    hasHit = hasHit_i
                if hasHit:
                    c += 1
        res.append((mid, c/s))
    return res


def relaxedScore(segments, marker_seqs, blastMap, onlyAll=False):
    res = list()
    for mid in sorted(segments.keys()):
        if onlyAll and len(set(x for _, x, _, _ in segments[mid])) < n:
            continue
        c = 0.0
        s = 0
        for i in xrange(len(segments[mid])):
            sid1, G1, _, _ = segments[mid][i]
            s += len(marker_seqs[sid1])
            for g1x in marker_seqs[sid1]: 
                hasHit = False
                j = 0
                while not hasHit and j < len(segments[mid]):
                    j += 1
                    if i == j-1:
                        continue
                    sid2, G2, _, _ = segments[mid][j-1]
                    if not blastMap.has_key((G1, G2)) or \
                            not blastMap[(G1, G2)].has_key(g1x):
                            break
                    for g2y in blastMap[(G1, G2)][g1x]:
                        hasHit = g2y in marker_seqs[sid2]
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
    plt.xlim([0, 1.005])
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
    parser.add_option('-a', '--only_all', dest='onlyAll', default=False,
            action='store_true', help='Consider only those multiplicons ' + \
                    'that span all genomes, not just a subset [default: ' + \
                    '%default]')
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
    
    cf = logging.FileHandler('%s.synteny_scores.log' %(basename(args[0]).rsplit('.', 1)[0]), mode='w')
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

    listElemFile = join(dirname(args[1]), 'list_elements.txt')

    if not isfile(listElemFile):
        LOG.fatal(('File %s is required, but does not exist. ' + \
                'Exiting.') %listElemFile)
        exit(1)

    marker_seqs = readMarkerSequences(open(listElemFile))
    _ , gene2genome = readGenomes(iadhoreConfig, dirname(args[0]))
    segments = readSegments(open(args[1]))
    
    blastMap = readBlastTable(open(join(dirname(args[0]),
        iadhoreCMap['blast_table'])), gene2genome)

    scores = None
    if options.type == 'relaxed':
        scores = relaxedScore(segments, marker_seqs, blastMap, options.onlyAll)
    elif options.type == 'weighted':
        scores = weightedScore(segments, marker_seqs, blastMap, options.onlyAll)
    else:
        LOG.fatal('Scoring type must be either relaxed or weighted. Exiting')
        exit(1)

    if options.visual:
        showHistogram(scores, options.type + ' synteny score', options.figName)
    else:
        for s in scores:
            print >> stdout, '%s\t%s' %s

    # compute coverage
    gNames = sorted(set(reduce(lambda x,y: x+y, blastMap.keys())))
    covered_markers = dict((G1, set()) for G1 in gNames)
    for seg in segments.values():
        if options.onlyAll and len(set(x for _, x, _, _ in seg)) < len(gNames):
            continue
        for sid, G1, _, _ in seg:
            covered_markers[G1].update(marker_seqs[sid])

    res = dict()
    for G1, markers in covered_markers.items(): 
        res[G1] = 0
        for g1i in markers:
            m = MARKER_PAT.match(g1i)
            _, start, end, _ = m.groups()
            res[G1] += int(end)-int(start)+1

    LOG.info('Coverage: \n%s' %('\n'.join(map(lambda x: '\t'.join(map(str, x)),
        sorted(res.items())))))

    LOG.info('finished')
