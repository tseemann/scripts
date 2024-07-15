#!/usr/bin/env python

"""

Changes:

0.2  - Added --titleX options.

"""

VERSION = '0.2'

import sys, re, numpy, itertools, os

USAGE = """\
Usage:

  fisher_diff.py [options] \\
      1st-evidence.txt \\
      2nd-evidence.txt \\
      significance_cutoff >output_file.csv

Options:

  --title1 "xyz"   - Column title for first evidence file in output
  --title2 "xyz"   - Column title for second evidence file in output

1st-evidence.txt and 2nd-evidence.txt should be *-evidence.txt files 
produced by shrimp_consensus.py

significance_cutoff can be in scientific notation, eg 1e-5
"""

def get_option_value(argv, option, conversion_function, default):
    argv = argv[:]
    value = default
    while True:
        try:
            location = argv.index(option)
        except ValueError: #Not found
            break
            
        if location == len(argv)-1 :
            raise Error('Option %s requires a paramter' % option)
        
        try:
            value = conversion_function(argv[location+1])
        except Exception:
            raise Error('Option for %s not in expected format' % option)
        
        del argv[location:location+2]

    return value, argv

def status(string):
    """ Display a status string. """
    
    sys.stderr.write('\r\x1b[K' + string)
    sys.stderr.flush()

def decode_evidence(desc):
    result = [ ]
    for item in desc.split():
        match = re.match('"([^"]*)"x([0-9]+)$', item)
        result.append( (match.group(1), int(match.group(2))) )
    return result

def read_file(filename):
    f = open(filename,'rb')
    f.readline()
    for line in f:
        parts = line.rstrip('\n').split('\t')
        yield int(parts[0]), parts[1], parts[2], parts[3]



LOG_FAC_CACHE = [ 0.0 ]
def log_fac(n):
    try:
        return LOG_FAC_CACHE[n]
    except IndexError:
        while len(LOG_FAC_CACHE) <= n:
            LOG_FAC_CACHE.append(LOG_FAC_CACHE[-1]+numpy.log(len(LOG_FAC_CACHE)))
        return LOG_FAC_CACHE[n]  

SUM_ERROR_MARGIN = 0.999999

class Cutoff_exceeded(Exception): pass

def fexact(matrix, significance_cutoff):
    matrix = numpy.asarray(matrix)
    n_row, n_col = matrix.shape
    row_sum = numpy.sum(matrix,1)
    col_sum = numpy.sum(matrix,0)
    n = numpy.sum(row_sum)

    cutoff = sum([ log_fac(item) for item in matrix.ravel() ]) * SUM_ERROR_MARGIN

    const_part = sum([ log_fac(item) for item in row_sum ]) + sum([ log_fac(item) for item in col_sum ]) - log_fac(n)
    
    significance = [ 0.0 ]
    
    row_remainders = row_sum.copy()
    def generate(row, col, col_remainder, total):
        cell_min = max(0, col_remainder - numpy.sum(row_remainders[row+1:]))
        cell_max = min(col_remainder, row_remainders[row])
        
        #next_row_remainders = row_remainders.copy()
        old_row_remainder = row_remainders[row]
        for i in xrange(cell_min,cell_max+1):
            next_total = total + log_fac(i)            
            row_remainders[row] = old_row_remainder - i
            if row+1 >= n_row:
                if col+1 >= n_col:
                    if next_total >= cutoff:
                        significance[0] += numpy.exp( const_part-next_total )
                        if significance[0] > significance_cutoff:
                            raise Cutoff_exceeded
                else:
                    generate(0, col+1, col_sum[col+1], next_total)
            else:
                generate(row+1, col, col_remainder-i, next_total)
        row_remainders[row] = old_row_remainder

    try:
        generate(0,0,col_sum[0],0.0)
    except Cutoff_exceeded:
        return None
    
    return significance[0]

SIG_CACHE = { }
def significance(ev1, ev2, cutoff):
    # Avoid duplicates in SIG_CACHE
    if ev2 < ev1:
        ev1,ev2 = ev2,ev1

    options = { }
    for item in ev1:
        options[item[0]] = len(options)
    for item in ev2:
        if item[0] not in options:
            options[item[0]] = len(options)
    
    n = len(options)
    
    matrix = numpy.zeros((n,2),'int')
    for item in ev1:
        matrix[ options[item[0]], 0 ] = item[1]
    for item in ev2:
        matrix[ options[item[0]], 1 ] = item[1]
    
    #s = fexact(matrix)
    #if s == 0:
    #    print matrix, s
    #return fexact(matrix)
    
    key = (tuple(matrix.ravel()),cutoff)
    if key not in SIG_CACHE:
        SIG_CACHE[key] = fexact(matrix,cutoff)    
    return SIG_CACHE[key]


def title_from_filename(filename):
    filename = os.path.abspath(filename)
    a,b = os.path.split(filename)
    a = os.path.split(a)[1]
    result = a + '/' + b
    if result.endswith('-evidence.txt'):
        result = result[:-13]
    return result


args = sys.argv[1:]
title1, args = get_option_value(args, '--title1', str, None)
title2, args = get_option_value(args, '--title1', str, None)

if len(args) != 3:
    print >> sys.stderr, USAGE
    sys.exit(1)

filename1 = args[0]
filename2 = args[1]
cutoff = float(args[2])

if title1 is None: title1 = title_from_filename(filename1)
if title2 is None: title2 = title_from_filename(filename2)

n = 1
while significance([('A',n)],[('T',n)],1.0) > cutoff:
    n += 1

print '%g,significance cutoff' % cutoff
print '%d,depth required to call substitution (greater if there are errors in the reads)' % n


print 'Position,Type,Reference,%s,%s,p-value (no correction for multiple testing)' % (title1, title2)

for (pos1, ins1, sub1, ref1), (pos2, ins2, sub2, ref2) in itertools.izip(read_file(filename1), read_file(filename2)):
    assert pos1 == pos2 and ref1 == ref2

    if pos1 % 1000 == 0:
        status('Testing %d' % pos1)

    #if pos1 not in (692366, 917937, 1349539, 1110346, 1411407): continue

    dec_ins1 = decode_evidence(ins1)
    dec_ins2 = decode_evidence(ins2)
    if dec_ins1 and dec_ins2:
        sig = significance(decode_evidence(ins1), decode_evidence(ins2), cutoff)    
        if sig is not None and sig <= cutoff:
            status('')
            print '%d,%s,"","%s","%s",%g' % (pos1, 'insertion-before', ins1.replace('"',"'"), ins2.replace('"',"'"), sig)

    dec_sub1 = decode_evidence(sub1)
    dec_sub2 = decode_evidence(sub2)
    if dec_sub1 and dec_sub2:
        sig = significance(dec_sub1, dec_sub2, cutoff)        
        if sig is not None and sig <= cutoff:
            if dec_sub1[0][0] == '-' or dec_sub2[0][0] == '-':
                what = 'deletion'
            elif dec_sub1[0][0] != dec_sub2[0][0]:
                what = 'substitution'
            else:
                what = 'different mix'
            status('')
            print '%d,%s,%s,"%s","%s",%g' % (pos1, what, ref1, sub1.replace('"',"'"), sub2.replace('"',"'"), sig)

status('')
