#!/usr/bin/env python

# eg blastall -p blastn -m 8 -d 295-consensus.fna -i 256-LIKE_seqA.dna |python ../blast2tag.py

import sys, optparse

option_parser = optparse.OptionParser("""usage: %prog [options] -t TYPE <blastn-output >tag-file

blastn-output:
    The consensus sequences should be the database (-d) and the thing you
    are searching for the query (-i). Output should be in tabular (-m 8)
    format. 
	
When importing into Gap4, you should set "Unpadded tag positions" to yes.""")

option_parser.add_option("-e", dest="e_value", metavar=" NNN", type="float",
                         help="e-value cutoff for blast tables, default no cutoff")

option_parser.add_option("-b", dest="bit_score", metavar=" NNN", type="float",
                         help="bit-score cutoff for blast tables, default no cutoff")

option_parser.add_option("-t", dest="type", metavar=" TYPE",
                         help="Tag type, eg COMM OLIG ENZ0 ENZ1")

options, args = option_parser.parse_args()

if len(args) != 0:
    option_parser.print_help(sys.stderr)
    sys.exit(0)

if not options.type:
    option_parser.error('A tag type is required, eg -t COMM') 

if len(options.type) != 4:
    option_parser.error('Tag type must be 4 letters')

for line in sys.stdin:
    fields = line.strip().split('\t')
    
    if (options.e_value and float(fields[10]) > options.e_value) or \
       (options.bit_score and float(fields[11]) < options.bit_score):
        continue
    
    query_id = fields[0]
    id = fields[1]
    q_start = int(fields[6])
    q_end = int(fields[7])
    start = int(fields[8])
    end = int(fields[9])
    if start < end:
        strand = '+'
    else:
        strand = '-'
	end, start = start, end
	
    print 'ID   %s' % id
    print 'TC   %s %s %d..%d' % (options.type, strand, start, end)
    print 'TC        Blast hit: %s %d..%d ' % (query_id, q_start, q_end)
    print 'TC        E-value: %s' % fields[10]
    print 'TC        Bit score: %s' % fields[11]


