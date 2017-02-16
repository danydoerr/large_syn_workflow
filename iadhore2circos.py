#!/usr/bin/env python

from sys import stdout, stderr, exit, maxint
from optparse import OptionParser
from itertools import product, combinations
from os.path import basename, dirname, join
from random import shuffle
import csv
import re

MARKER_PAT = re.compile('^([^:]+):(\d+):(\d+)(\+$|-$|$)')

BREWER_COL = ['blues-%s-seq-3', 'bugn-%s-seq-3', 'bupu-%s-seq-3', 'gnbu-%s-seq-3', \
        'greens-%s-seq-3', 'greys-%s-seq-3', 'oranges-%s-seq-3', 'orrd-%s-seq-3', \
        'pubu-%s-seq-3', 'pubugn-%s-seq-3', 'purd-%s-seq-3', 'purples-%s-seq-3', \
        'rdpu-%s-seq-3', 'reds-%s-seq-3', 'ylgn-%s-seq-3', 'ylgnbu-%s-seq-3',\
        'ylorbr-%s-seq-3', 'ylorrd-%s-seq-3']
BREWER_COL_RANGE = range(3,10)

TICKS_CONF = ''' 
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p

# the tick label is derived by multiplying the tick position
# by 'multiplier' and casting it in 'format':
#
# sprintf(format,position*multiplier)
#

multiplier       = 1e-6

# %d   - integer
# %f   - float
# %.1f - float with one decimal
# %.2f - float with two decimals
#
# for other formats, see http://perldoc.perl.org/functions/sprintf.html

format           = %d

<tick>
spacing        = 5u
size           = 10p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>

</ticks>
'''

IDEOGRAM_CONF = '''
<ideogram>

<spacing>
default = 0.005r
</spacing>

# Ideogram position, fill and outline
radius           = 0.92r
thickness        = 100p
# fill             = yes
# fill_color       = black
# stroke_color     = black
stroke_thickness = 0

show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 0
# band_stroke_color     = white
band_transparency     = 1

# Minimum definition for ideogram labels.

show_label       = yes
# see etc/fonts.conf for list of font names
label_font       = default 
label_radius     = 1.045r  # if ideogram radius is constant, and you'd like labels close to image edge, 
                           # use the dims() function to access the size of the image
                           # label_radius  = dims(image,radius) - 60p
label_size       = 60
label_parallel   = yes

</ideogram>
'''

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

def readGenomes(iadhoreConfig, configPath):

    gNames = list()
    genomes = list()
    for x in iadhoreConfig:
        if len(x) != 3 or x[0] != 'genome':
            continue
        gNames.append(x[1])
        genomes.append(dict())
        for chr1, fPath in x[2]:
            f = join(configPath, fPath)
            genomes[-1][chr1] = list()
            for line in open(f):
                m = MARKER_PAT.match(line)
                if not m:
                    print >> stderr, ('Unable to parse marker ' + \
                            'identifier \"%s\" in file %s. Exiting') %(line.strip(), f)
                    #exit(1)
                    gid = line[:-2]
                    strand = line[-2]
                    start, end = 0, 0
                else:
                    gid, start, end, strand = m.groups()
                genomes[-1][chr1].append((gid, int(start), int(end), strand))
    return gNames, genomes



def readSegments(data):

    isHeader = True
    
    res = dict()
    c = 0
    for line in csv.reader(data, delimiter='\t'):
        c += 1
        if isHeader:
            isHeader = False
            continue

        if not res.has_key(line[1]):
            res[line[1]] = dict()
        if not res[line[1]].has_key(line[2]):
            res[line[1]][line[2]] = list()

        m1 = MARKER_PAT.match(line[4])
        m2 = MARKER_PAT.match(line[5])
        if not m1 or not m2:
            print >> stderr, ('Unable to parse marker identifier ' + \
                    '\"%s\" and/or \"%s\" in line %s of file %s. ' + \
                    'Exiting') %(line[4], line[5], c, data.name)
            #exit(1)
            start1, end2 = line[4:6]
        else:
            gid1, start1, end1, _ = m1.groups()
            gid2, start2, end2, _ = m2.groups()
            start1, end1, start2, end2 = map(int, (start1, end1, start2,
                end2))

            if start1 > start2:
                gid1, start1, end1, gid2, start2, end2 = gid2, start2, \
                        end2, gid1, start1, end1

        res[line[1]][line[2]].append((line[3], start1, end2))

    return res


def readCloudAP(gene2genome, data):
    isHeader = True
    res = dict()

    for line in csv.reader(data, delimiter='\t'):
        if isHeader:
            isHeader = False
            continue

        if not res.has_key(line[0]):
            res[line[0]] = dict()

        g1i, g2k = line[1:3]
        (Gx,chrx), (Gy,chry) = map(gene2genome.get, (g1i, g2k))

        m1 = MARKER_PAT.match(line[1])
        m2 = MARKER_PAT.match(line[2])
        if not m1 or not m2:
            print >> stderr, ('Unable to parse marker identifier ' + \
                    '\"%s\" and/or \"%s\" in line %s of file %s. ' + \
                    'Exiting') %(line[4], line[5], c, data.name)
            #exit(1)
            start1, end2 = line[4:6]
        else:
            gid1, start1, end1, _ = m1.groups()
            gid2, start2, end2, _ = m2.groups()
            start1, end1, start2, end2 = map(int, (start1, end1, start2, end2))

            if not res[line[0]].has_key(Gx):
                res[line[0]][Gx] = list()

            l = res[line[0]][Gx]
            p = [x for x in xrange(len(l)) if l[x][0] == chrx]
            if not p:
                l.append((chrx, start1, end1))
            else:
                l[p[0]] = (chrx, min(start1, l[p[0]][1]), max(end1,
                    l[p[0]][2]))

            if not res[line[0]].has_key(Gy):
                res[line[0]][Gy] = list()

            l = res[line[0]][Gy]
            p = [x for x in xrange(len(l)) if l[x][0] == chry]
            if not p:
                l.append((chry, start2, end2))
            else:
                l[p[0]] = (chry, min(start2, l[p[0]][1]), max(end2,
                    l[p[0]][2]))
    return res 
                    

def processMultiplicons(data, segments, allLevels, n, onlyAll=False):

    res = dict()
    isHeader = True
    for line in csv.reader(data, delimiter='\t'):
        if isHeader:
            isHeader = False
            continue

        mid = line[0]
        if not allLevels and int(line[6]) > 2:
            continue

        if not segments.has_key(mid):
            print >> stderr, ('Multiplicon with ID %s not found in ' + \
                    'segments file. Exiting') %mid
            exit(1)

        if onlyAll and len(segments[mid]) < n:
            continue

        level=1
        if len(line) > 5:
            try:
                level = int(line[6])
            except:
                pass

        for G1, G2 in combinations(sorted(segments[mid].keys()), 2):
            if not res.has_key((G1, G2)):
                res[(G1, G2)] = list()
            for (chr1, start1, end1), (chr2, start2, end2) in \
                    product(segments[mid][G1], segments[mid][G2]):
                id1 = '%s.%s' %(G1.lower(), chr1.lower())
                id2 = '%s.%s' %(G2.lower(), chr2.lower())
                res[(G1, G2)].append((id1, start1, end1, id2, start2, end2,
                    mid, level))

    return res

def writeKaryotype(genome, gName, out):

    glow = gName.lower()

    for chr1 in sorted(genome.keys()):
        chrlow = chr1.lower()
        print >> out, 'chr - %s.%s %s.%s 0 %s %s.%s' %(glow, chrlow, gName,
                chr1, genome[chr1][-1][2]+1, glow, chrlow)

    for chr1, seq1 in genome.items():
        chrlow = chr1.lower()

        prev_end = -1
        c = 0
        for gid, start, end, _ in seq1:
            if prev_end + 1 < start:
                print >> out, ('band %s.%s %s.%s.no_seg%s %s.%s.no_seq%s ' + \
                        '%s %s %s.no_seg') %(glow, chrlow, glow, chrlow, c,
                                glow, chrlow, c, prev_end + 1, start, glow)
            print >> out, 'band %s.%s %s.%s %s.%s %s %s %s.seg' %(glow, chrlow,
                    glow, gid.lower(), gName, gid, start, end+1, glow)
            c += 1
            prev_end = end


def writeCircosConf(g1name, g1, g2name, g2, multiplicons, sbs, out):

    g1low = g1name.lower()
    g2low = g2name.lower()

    chr1s = map(lambda x: '.'.join(x), zip([g1low] * len(g1), sorted(map(lambda
        x: x.lower(), g1.keys()))))
    chr2s = map(lambda x: '.'.join(x), zip([g2low] * len(g2), sorted(map(lambda
        x: x.lower(), g2.keys()))))

    print >> out, ('karyotype = karyotype.%s.txt,karyotype.%s.txt') %(g1name,
            g2name)
    print >> out, 'chromosomes_units = 1000000'
    print >> out, 'chromosomes_display_default = no' 
    print >> out, 'chromosomes                 = %s' %(';'.join(chr1s + chr2s))

    print >> out, '<colors>'
    for c in chr1s:
        print >> out, '%s = blue' %c
    for c in chr2s:
        print >> out, '%s = green' %c

    print >> out, '%s.seg = black' %g1low
    print >> out, '%s.no_seg = blue' %g1low
    print >> out, '%s.seg = black' %g2low
    print >> out, '%s.no_seg = green' %g2low

    c_max = len(BREWER_COL_RANGE) * len(BREWER_COL)
    if options.colorLinks:
        shuffle(multiplicons)
        for x in xrange(len(multiplicons)):
            r = (x % c_max) / len(BREWER_COL)
            i = (x % c_max) % len(BREWER_COL)
            print >> out, 'color.%s = %s' %(multiplicons[x], BREWER_COL[i] %
                    BREWER_COL_RANGE[r])
    else:
        c_min = len(chr1s) < len(chr2s) and chr1s or chr2s
        if len(c_min) > len(BREWER_COL_RANGE) * len(BREWER_COL):
            print >> stderr, '!! Too many chromosomes. Not enough colors ' + \
                    'available to color all links! Exiting'
            exit(1)

        x = 0
        for c in c_min:
            r = x / len(BREWER_COL)
            i = x % len(BREWER_COL)
            print >> out, '%s.l = %s' %(c, BREWER_COL[i] %BREWER_COL_RANGE[r])
            x += 1

    print >> out, '</colors>'
    print >> out, '<links>\n<link>'
    print >> out, 'file          = %s_%s.links' %(g1name, g2name)
    print >> out, 'radius        = 0.99r'
    print >> out, 'bezier_radius = 0r'
    print >> out, 'ribbon = yes'
    print >> out, 'flat   = yes'
    print >> out, '<rules>\n<rule>'
    print >> out, 'condition     = 1'
    if options.colorLinks:
        print >> out, 'color         = eval(\'color.\' . var(mid) . ' + \
                '\'_a\' .  var(level))' 
    else:
        print >> out, 'color         = eval(var(chr%s) . \'.l_a3\')' %(c_min ==
                chr1s and 1 or 2)
    print >> out, 'flow          = continue'
    print >> out, '</rule>\n</rules>'
    print >> out, '</link>\n</links>'
    print >> out, '<<include ideogram.conf>>'
    print >> out, '<<include ticks.conf>>'
    print >> out, '<image>\nfile  = %s.png' %(basename(out.name).split('.', 1)[0])
    print >> out, 'dir   = .\npng   = yes\nradius         = 1500p' 
    print >> out, 'angle_offset      = -90'
    print >> out, 'auto_alpha_colors = yes\nauto_alpha_steps  = 5'
    print >> out, 'background = white\n</image>'
    print >> out, '<<include etc/colors_fonts_patterns.conf>>'
    print >> out, '<<include etc/colors.brewer.conf>>'
    print >> out, '<<include etc/housekeeping.conf>>'
    print >> out, 'data_out_of_range* = warning'


if __name__ == '__main__':
    usage = 'usage: %prog [options] <I-ADHORE CONFIGURAION FILE> ' + \
            '[<MULTIPLICONS FILE> <SEGMENTS FILE>] OR [<CLOUDS FILE> ' + \
            '<CLOUD_AP FILE>]'
    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--color_links', dest='colorLinks', default=False,
            action='store_true', help='Color links individually')
    parser.add_option('-o', '--outdir', dest='outDir', default='.', type=str, 
            help='Direction to which the output files will be written. ' + \
                    '[default=%default]')
    parser.add_option('-a', '--only_all', dest='onlyAll', default=False,
            action='store_true', help='Consider only those multiplicons ' + \
                    'that span all genomes, not just a subset [default: ' + \
                    '%default]')
    parser.add_option('-l', '--plot_all_levels', dest='allLevels', default=False,
            action='store_true', help='Not only draw level 2 intervals, ' + \
                    'but also intervals from higher levels [default: %default]')

    (options, args) = parser.parse_args()
    
    if len(args) != 3:
        parser.print_help()
        exit(1)

    iadhoreConfig = readIadhoreConfig(open(args[0]))
    gNames, genomes = readGenomes(iadhoreConfig, dirname(args[0]))
    #
    # write ticks & ideogram conf 
    #

    ticks_out = open(join(options.outDir, 'ticks.conf'), 'w')
    ticks_out.write(TICKS_CONF)
    ticks_out.close()

    ideo_out = open(join(options.outDir, 'ideogram.conf'), 'w')
    ideo_out.write(IDEOGRAM_CONF)
    ideo_out.close()
    #
    # write karyotypes
    #

    for i in xrange(len(genomes)):
        out = open(join(options.outDir, 'karyotype.%s.txt' %gNames[i]), 'w')
        writeKaryotype(genomes[i], gNames[i], out)
        out.close() 

    #
    # write links and circos conf
    #

    if basename(args[2]).lower() == 'cloudap.txt':

        gene2genome = dict()
        for x in xrange(len(genomes)):
            Gx = gNames[x]
            for chrx, genes in genomes[x].items():
                for g in genes:
                    gid = '%s:%s:%s' %g[:3]
                    gene2genome[gid] = (Gx, chrx)

        segments = readCloudAP(gene2genome, open(args[2]))
        syn_blocks = processMultiplicons(open(args[1]), segments, True,
                len(gNames), options.onlyAll)

    else:
        segments = readSegments(open(args[2]))
        syn_blocks = processMultiplicons(open(args[1]), segments,
                options.allLevels, len(gNames), options.onlyAll)

    gname2id = dict(zip(gNames, range(len(gNames))))
    for (G1, G2), sbs in syn_blocks.items():
        out = open(join(options.outDir, '%s_%s.links' %(G1, G2)), 'w')
        sbs.sort()
        for x in sbs:
            print >> out, ' '.join(map(str, x[:-2])), ('mid=%s,' + \
                    'level=%s') %(x[-2], min(x[-1], 5))
        out.flush()
        out.close()
   
        out = open(join(options.outDir, '%s_%s.circos.conf' %(G1, G2)), 'w')
        writeCircosConf(G1, genomes[gname2id[G1]], G2, genomes[gname2id[G2]],
                segments.keys(), sbs, out)
