#!/usr/bin/env python2.6

"""

SHRiMP read consensus caller

To call a consensus we require:
- greater than a certain depth (to be counted in depth, must match the consensus)
- winner to comprise a certain percentage of all (purity)


Changes:
0.1 - initial version

0.2 - use paired end information
      output a full "evidence" dump
      added --trim option

0.3 - DNA ambiguity codes in consensus
      don't enforce purity on insertions

0.4 - added options to allow SOLiD paired end reads
      include single reads when processing a paired end run

0.5 - replaced "has_consensus" file with "reference_having_consensus" file
      added alignments.maf output

0.6 - added option --infidelty
      fixed bug that caused it to declare more read pairs ambiguous than it should have

0.7 - use Cython if present

0.8 - output consensus_masked.fa
      warn if reads look like paired and but no --max-pair-sep given

0.9 - more efficient use of memory (merely scan hits and record location in file)

0.10 - Added option --ambiguity-codes
       Reads not following pair naming convention will no longer cause an error

Note: Python's garbage collection means allocation is O(n^2) (still true as at Python 2.6.1)
      and is disabled in this script
      (this still leaves reference counting, so there should not be leaks)

"""

VERSION = '0.10'


import sys, os, numpy, heapq, string, math
from Bio import Seq, SeqRecord, SeqIO, Alphabet


USAGE = """\

You can also specify:

  --solid

to use SOLiD read name suffixes and directions (the defaults 
work with Illumina paired end). You will still need to specify 
"--max-pair-sep NNN" to enable paired-end processing. 

Usage: 

  shrimp_consensus.py [options] working_dir

where working_dir was created by shrimp_run.py.

"""

class Error(Exception): 
    pass

def status(string):
    """ Display a status string. """
    
    sys.stderr.write('\r\x1b[K' + string)
    sys.stderr.flush()

COMPLEMENTER = string.maketrans('ACGTacgt','TGCAtgca')

def reverse_complement(seq):
    """ Reverse complement of a string, works with ACGTN- """
    
    return seq.translate(COMPLEMENTER)[::-1]
    #return Seq.Seq(seq, Alphabet.IUPAC.ambiguous_dna).reverse_complement().tostring()


def filesystem_friendly_name(name):
    """ Remove special characters from a name """

    for char in '\'"<>&|/\\_ .':
        name = name.replace(char,'_')
    return name

def open_possibly_compressed_file(filename):
    if filename.endswith('.gz'):
        #return gzip.open(filename,'rb')
        #TODO: error checking
        quoted = "'" + filename.replace("'","'\\''") + "'"
        return os.popen('gunzip <' + quoted, 'rb')
    else:
        return open(filename,'rb')
    
def edit_string_to_alignment(edit_string, corresp_seq):
    """ Convert an edit string as output by SHRiMP and ELAND to
        a conventional alignment. """

    len_edit_string = len(edit_string)
    read_ali = ''
    contig_ali = ''
    i = 0    
    while i < len_edit_string:
        if edit_string[i] == '(':
            j = i + 1
            while edit_string[j] != ')': j += 1
            gap = edit_string[i+1:j].upper()
            
            #print 'insertion', gap
            read_ali += gap
            contig_ali += '-' * len(gap)
            
            i = j + 1
        elif edit_string[i].isdigit():
            j = i
            while j < len_edit_string and edit_string[j].isdigit(): j += 1
            n_matches = int(edit_string[i:j].upper())

            #print 'match', n_matches            
            read_ali += corresp_seq[:n_matches]
            contig_ali += corresp_seq[:n_matches]
            corresp_seq = corresp_seq[n_matches:]

            i = j            
        elif edit_string[i] == '-':
            #print 'deletion'
            j = i
            while j < len_edit_string and edit_string[j] == '-': j += 1
            
            n_deleted = j-i
            contig_ali += corresp_seq[:n_deleted]
            read_ali += '-' * n_deleted
            corresp_seq = corresp_seq[n_deleted:]
        
            i = j
        else:
            #print 'substitution', edit_string[i]
            contig_ali += corresp_seq[:1]
            read_ali += edit_string[i].upper()
            corresp_seq = corresp_seq[1:]
            
            i += 1

    return contig_ali, read_ali

def roll_alignment(ali, ali_other):
    """ Normalize one half of an alignment by shifting "-"s as far right as possible.
    
        (SHRiMP alignments can differ depending on whether the read is forward 
         or reverse, when the location of an insertion or deletion is ambiguous.)
    """

    if '-' not in ali: return ali

    len_ali = len(ali)

    i = 0
    while i < len_ali:
        if ali[i] != '-': 
            i += 1
            continue
            
        j = i+1
        while j < len_ali and ali[j] == '-': j += 1
        if j >= len_ali: 
            break
            
        if ali_other[i] == ali_other[j] or ali[j] != ali_other[j]:
            ali = ali[:i] + ali[j] + ali[i+1:j] + '-' + ali[j+1:]            
            i = i + 1
        else:
            i = j
            
    return ali


def _get1(x): return x[1]

def consensus(counts, min_depth,min_purity):
    """ Call a consensus, if it meets minimum depth and purity 
        (proportion of total) requirements. """

    if not counts: return None

    total = sum(counts.values())
    counts = counts.items()
    counts.sort(key=_get1)
    
    if counts[-1][1] < min_depth or \
       counts[-1][1] < min_purity * total: 
        return None
    
    return counts[-1][0]


AMBIGUITY_CODES = {
    '-' : '-',
    'A' : 'A',
    'T' : 'T',
    'G' : 'G',
    'C' : 'C',
    'GT' : 'K',
    'AC' : 'M',
    'AG' : 'R',
    'CT' : 'Y',
    'CG' : 'S',
    'AT' : 'W',
    'CGT' : 'B',
    'ACG' : 'V',
    'ACT' : 'H',
    'AGT' : 'D',
}

AMBIGUITY_DECODE = { }
for decode in AMBIGUITY_CODES:
    AMBIGUITY_DECODE[ AMBIGUITY_CODES[decode] ] = decode

def ambiguity_code_consensus(counts, min_depth,min_purity):
    """ Call a consensus, if it meets minimum depth and purity 
        (proportion of total) requirements. 
        
        Allow DNA ambiguity codes. 
        
        """

    if not counts: return None

    total = sum(counts.values())
    counts = counts.items()
    counts.sort(key=_get1, reverse=True)
    
    total_count = 0
    bases = [ ]
    cutoff = None #<- if multiple items have the same count as that required to reach purity/depth
                  #   be sure to include them all in the ambiguity code
    for base, count in counts:
        if cutoff is not None and count < cutoff:
            break
    
        total_count += count
        bases.append(base)
        
        if total_count >= min_depth and \
           total_count >= min_purity * total: 
            cutoff = count
    
    if cutoff is None:
        return None

    bases.sort()
    bases = ''.join(bases)
    
    return AMBIGUITY_CODES.get(bases, None)

def entropy(counts):
    total = sum(counts.values())
    result = 0.0
    for value in counts.values():
        if value == 0: continue
        p = float(value) / total
        result -= p * math.log(p)
    return result / math.log(2.0)


#def iter_increasing_pairs(list1, list2, scorer=lambda a,b: a+b):
#    """ Yield pairs from list1 and list2, 
#        in increasing order of combined score.
#    
#        lists must be sorted and score combiner
#        must be monotonic increasing in both parameters. """
#
#    if not list1 or not list2: return
#    
#    heap = [ (scorer(list1[0],list2[0]),0,0) ]
#    while heap:
#        score, pos1, pos2 = heapq.heappop(heap)        
#        yield score, pos1, pos2
#        
#        if not pos2 and pos1+1 < len(list1):
#            heapq.heappush(heap,
#                (scorer(list1[pos1+1],list2[pos2]), pos1+1, pos2) )
#        if pos2+1 < len(list2):
#            heapq.heappush(heap,
#                (scorer(list1[pos1],list2[pos2+1]), pos1, pos2+1) )

class Hit_pair_iter(object):
    """ Yield pairs of hits from list1 and list2, 
        in decreasing order of combined score.
    
        lists must be sorted. """
        
    def __init__(self, list1, list2):
        if not list1 or not list2:
            self.heap = [ ]
        else:
            self.heap = [ (-list1[0].score-list2[0].score,0,0) ]
        self.list1 = list1
        self.list2 = list2
    
    def __iter__(self):
        return self
    
    def next(self):
        if not self.heap:
            raise StopIteration()
    
        neg_score, pos1, pos2 = heapq.heappop(self.heap)        
        
        if not pos2 and pos1+1 < len(self.list1):
            heapq.heappush(self.heap,
                (-self.list1[pos1+1].score-self.list2[pos2].score, pos1+1, pos2) )
        if pos2+1 < len(self.list2):
            heapq.heappush(self.heap,
                (-self.list1[pos1].score-self.list2[pos2+1].score, pos1, pos2+1) )

        return -neg_score, pos1, pos2

def pretty_number(number, width=0):
    """ Adds commas for readability. """
    assert isinstance(number,int)    
    s = str(abs(number))[::-1]
    groups = [ s[i:i+3] for i in xrange(0,len(s),3) ]
    result = ','.join(groups)
    if number < 0: 
        result += '-'
    if len(result) < width:
        result += ' '*(width-len(result))
    return result[::-1] 

def pretty_evidence(counts):
    counts = counts.items()
    counts.sort(key=_get1, reverse=True)
    return ' '.join([ '"%s"x%d' % item for item in counts ])

def _hit_score(hit): return hit.score

class Refseqset(object):
    """ A collection of reference seqeunces.
    
        Properties:
    
            seqs - { name : Refseq }
            hits - { read_name : [ Hit ] }

    """

    def __init__(self):
        self.seqs = { }
        self.hits = { }
        self.shrimp_file = None
    
    def add_sequence(self, name, sequence):
        assert name not in self.seqs, 'Duplicate reference sequence name'
        
        self.seqs[name] = Refseq(name, sequence)

    def read_shrimp(self, filename, warn_about_paired_end=False):
        assert self.shrimp_file is None
        self.shrimp_file = open(filename,'rb')
        
        n = 0        
        while True:
            position = self.shrimp_file.tell()
            line = self.shrimp_file.readline()
            if not line: break
            
            if line.startswith('#'): continue
            
            n += 1
            if n % 10000 == 0:
                status('Scanning hit %s of %s' % (pretty_number(n), filename))
            
            read_name = line[:line.index('\t')].split()[0]
            if read_name not in self.hits: 
                self.hits[read_name] = [ ]
            self.hits[read_name].append(position)
            
            if warn_about_paired_end and \
               (read_name[-3:] in ('_F3','_R3') or
                read_name[-2:] in ('/1','/2')):
                status('')
                
                sys.stderr.write('\n*** WARNING: read names look like paired end, but no --max-pair-sep given ***\n\n')
                warn_about_paired_end = False 
                
        status('')

    def get_hits(self, read_name):
        result = [ ]
        for position in self.hits[read_name]:
            self.shrimp_file.seek(position)
            line = self.shrimp_file.readline()
            
            (read_name, contig_name, strand, 
             contig_start, contig_end, 
             read_start, read_end, read_length, 
             score, edit_string) = line.rstrip().split('\t')

            hit = Hit()
        
            hit.read_name = read_name.split()[0]
            hit.ref_name = contig_name.split()[0]
            hit.ref_start = int(contig_start)-1
            hit.ref_end = int(contig_end)
            hit.read_start = int(read_start)-1
            hit.read_end = int(read_end)
            hit.read_length = int(read_length)
            hit.score = int(score)
            hit.forward = (strand == '+')
            
            assert hit.read_name == read_name
            
            #if not read_name.endswith('/1') and not read_name.endswith('/2'):
            #    print repr(line)
            
            corresp_seq = self.seqs[hit.ref_name].reference[hit.ref_start:hit.ref_end]
    
            if not hit.forward:
                corresp_seq = reverse_complement(corresp_seq)

            hit.ref_ali, hit.read_ali = edit_string_to_alignment(edit_string, corresp_seq)	    
            
            if not hit.forward:
                hit.ref_ali = reverse_complement(hit.ref_ali)
                hit.read_ali = reverse_complement(hit.read_ali)

            #Normalization -- move "-"s as far right as possible
            hit.read_ali = roll_alignment(hit.read_ali, hit.ref_ali)
            hit.ref_ali = roll_alignment(hit.ref_ali, hit.read_ali)
            
            result.append(hit)
                
        return result

    def process_hits(self, trim, infidelity, read_names=None):
        """ Update counts in Refseq objects based on hits read in. """
        
        if read_names is None: 
            read_names = self.hits
        
        n_total = 0
        n_unambiguous = 0
        
        for i, read_name in enumerate(read_names):
            if (i % 10000) == 0:
                status('Processing read %s of %s' % (pretty_number(i), pretty_number(len(read_names))))
        
            #hits = self.hits[read_name]
            hits = self.get_hits(read_name)
            hits.sort(key=_hit_score, reverse=True)
            
            n = 1
            while n < len(hits) and hits[n].score >= hits[0].score*infidelity:
                n += 1
                
            if n == 1:
                self.seqs[ hits[0].ref_name ].process_unambiguous_hit(hits[0], trim)
            
            for hit in hits[:n]:
                self.seqs[ hit.ref_name ].process_hit(hit, n, trim)

        status('')
    
    def process_paired_hits(self, max_pair_sep, same_direction, suffix1, suffix2, trim, infidelity):
        """ Update counts in Refseq objects based on hits read in. """
        
        n_total = 0
        n_valid_total = 0        
        n_valid_unambiguous = 0
        unambiguous_seps = [ ]
        
        weird_pair_report = [ ]
        
        orphans = [ ]
        unpaired = [ ]
        
        for i, read_name_1 in enumerate(self.hits):
            if (i % 10000) == 0:
                status('Processing reads as pairs %s of %s' % (pretty_number(i), pretty_number(len(self.hits))))
                
            if read_name_1.endswith(suffix2): 
                read_name_2 = read_name_1[:len(read_name_1)-len(suffix2)] + suffix1
                if read_name_2 not in self.hits:
                    orphans.append(read_name_1)
                continue
            
            #assert read_name_1.endswith(suffix1), 'Paired read name does not follow given naming convention: %s' % read_name_1
            if not read_name_1.endswith(suffix1):
                unpaired.append(read_name_1)
                continue
            
            read_name_2 = read_name_1[:len(read_name_1)-len(suffix1)] + suffix2
            if read_name_2 not in self.hits:
                #TODO: report / treat as singleton
                # (also for /2 but not /1)
                orphans.append(read_name_1)
                continue
            
            n_total += 1
            
            #hits_1 = self.hits[read_name_1]
            hits_1 = self.get_hits(read_name_1)
            hits_1.sort(key=_hit_score, reverse=True)
            #hits_2 = self.hits[read_name_2]
            hits_2 = self.get_hits(read_name_2)
            hits_2.sort(key=_hit_score, reverse=True)
            
            pair_iter = Hit_pair_iter(hits_1,hits_2)
            
            top_hits = [ ]
            for score, pos_1, pos_2 in pair_iter:
                if top_hits and score < top_hits[0][0]*infidelity: break
                
                hit_1 = hits_1[pos_1]
                hit_2 = hits_2[pos_2]
                
                if hit_1.ref_name != hit_2.ref_name:
                    continue
                
                if same_direction:
                    if hit_1.forward != hit_2.forward: continue
                else:
                    if hit_1.forward == hit_2.forward: continue
                    
                if hit_1.forward:
                    sep = hit_2.ref_start - hit_1.ref_end
                else:
                    sep = hit_1.ref_start - hit_2.ref_end
                
                if sep > max_pair_sep or sep < 0:
                    continue
                
                top_hits.append((score, hit_1, hit_2))
            
            if not top_hits:
                if (len(hits_1) == 1 or hits_1[0].score*infidelity > hits_1[1].score) and \
                   (len(hits_2) == 1 or hits_2[0].score*infidelity > hits_2[1].score):
                    hit_1 = hits_1[0]
                    hit_2 = hits_2[0]
                    if hit_1.forward:
                        hit_1_end = hit_1.ref_end
                    else:
                        hit_1_end = hit_1.ref_start
                    if hit_2.forward:
                        hit_2_end = hit_2.ref_end
                    else:
                        hit_2_end = hit_2.ref_start
                    #out by one error above?
                    
                    weird_pair_report.append((hit_1.ref_name, hit_1.forward, hit_1_end, hit_2.ref_name, hit_2.forward, hit_2_end))
            
                continue
            
            n_valid_total += 1
            
            if len(top_hits) == 1:
                n_valid_unambiguous += 1

                hit_1 = top_hits[0][1]
                hit_2 = top_hits[0][2]
                if hit_1.forward:
                    sep = hit_2.ref_start - hit_1.ref_end
                else:
                    sep = hit_1.ref_start - hit_2.ref_end
                unambiguous_seps.append(sep)
                
                self.seqs[ hit_1.ref_name ].process_unambiguous_hit( hit_1, trim )
                self.seqs[ hit_2.ref_name ].process_unambiguous_hit( hit_2, trim )
        
            for neg_score, hit_1, hit_2 in top_hits:
                self.seqs[ hit_1.ref_name ].process_hit(hit_1, len(top_hits), trim)
                self.seqs[ hit_2.ref_name ].process_hit(hit_2, len(top_hits), trim)

        status('')
        
        self.process_hits(trim, infidelity, orphans + unpaired)
        
        stats_text = (        
           pretty_number(n_total,20) + ' read pairs where both reads hit something\n' +
           pretty_number(n_valid_total,20) + ' read pairs validly oriented and spaced\n' +
           pretty_number(n_valid_unambiguous,20) + ' unambiguously\n'
        )
        
        if unambiguous_seps:
            stats_text += pretty_number(int(numpy.median(unambiguous_seps)),20) + \
                          ' median separation of paired reads (limit %d)\n' % max_pair_sep

        stats_text += '\n' + pretty_number(len(orphans),20) + ' reads with no hits to their pair\n'
        
        if unpaired:
            stats_text += '\n' + pretty_number(len(unpaired),20) + ' without a read-pair suffix, treated as single reads\n'
        
        return stats_text, weird_pair_report

class Hit(object):
    """ Structure to store hits """

    __slots__ = (
        'ref_name',
        'read_name',
        'score',
        'forward',
        'ref_start',
        'ref_end',
        'ref_ali',
        'read_start',
        'read_end',
        'read_ali',
        'read_length',
    )

    
class Refseq(object):
    """ Reference sequence and statistics on alignments thereto.
    
        Properties:

            name
            reference - reference sequence
        
            depth - unambiguous hit depth
            depth_ambiguous - estimated true depth
        
            base_counts - { letter : count } (includes deletions as '-')
            insertions - { position : [ sequence ] }
                
    """

    def __init__(self, name, reference):
        self.name = name
        self.reference = reference

        self.depth = numpy.zeros(len(self.reference), 'int')
        self.depth_ambiguous = numpy.zeros(len(self.reference), 'float')
        
        #self.incomplete_ends_count = numpy.zeros(len(self.reference), 'int')
        #self.bias = numpy.zeros(len(self.reference), 'float')

        #self.base_counts = { }        
        #for base in 'ACGT-':
        #    self.base_counts[base] = numpy.zeros(len(self.reference), 'int')
            
        #self.insertions = { } # index *after* insertion -> { inserted sequence : count }
        
        self.base_counts = [ {} for i in xrange(len(reference)) ]
        self.insertions = [ {} for i in xrange(len(reference)) ]

    def process_unambiguous_hit(self, hit, trim):
        #self.depth[hit.ref_start:hit.ref_end] += 1
        
        #if hit.read_start > trim:
        #    if hit.forward:
        #        self.incomplete_ends_count[hit.ref_start] += 1
        #    else:
        #        self.incomplete_ends_count[hit.ref_end-1] += 1
        #if hit.read_end < hit.read_length-trim:
        #    if hit.forward:
        #        self.incomplete_ends_count[hit.ref_end-1] += 1
        #    else:
        #        self.incomplete_ends_count[hit.ref_start-1] += 1
        
        i = 0
        position = hit.ref_start
        
        hit_ref_ali = hit.ref_ali
        hit_read_ali = hit.read_ali
        len_hit_ref_ali = len(hit_ref_ali)
        
        while i < trim:
            if hit_ref_ali[i] != '-':
                position += 1
            i += 1

        #scaler = 1.0/(len(hit.ref_ali)-1.0)
        while i < len_hit_ref_ali-trim:
            if hit_ref_ali[i] == '-':
                j = i + 1
                while j < len_hit_ref_ali and hit_ref_ali[j] == '-': j += 1
                what = hit_read_ali[i:j]
                
                #if position not in self.insertions: self.insertions[position] = { }
                counter = self.insertions[position]
                counter[what] = counter.get(what,0) + 1
                
                i = j	        
            else:
                if hit_read_ali[i] not in 'NX': #X is a shrimp colorspace confuzzle
                    #self.base_counts[ hit_read_ali[i] ][position] += 1
                    counter = self.base_counts[position]
                    base = hit_read_ali[i]
                    counter[base] = counter.get(base,0)+1
                self.depth[position] += 1
                
                #self.bias[position] += i*scaler - 0.5 
                
                position += 1        
                i += 1
    
    def process_hit(self, hit, one_of_n, trim):
        #TODO: don't ignore trim!
        #self.depth_ambiguous[hit.ref_start:hit.ref_end] += 1.0 / one_of_n
        self.depth_ambiguous[hit.ref_start:hit.ref_end] = self.depth_ambiguous[hit.ref_start:hit.ref_end] + 1.0 / one_of_n

    def consensus(self, min_depth,min_purity,use_ambiguity_codes):
        """ 
        Returns:        
          consensus sequence        
          snp/indel report        
          For each position in the reference, whether there was a consensus        
        """
    
        result = [ ]
        result_masked_only = [ ]
        
        alignment_result    = [ ]
        alignment_reference = [ ]
        
        insertion_evidence = [ ]
        substitution_evidence = [ ]
        
        has_consensus = [ ]
        
        status('Consensus %s' % self.name)
        
        report = [ ]
        
        for i in xrange(len(self.reference)):
            if i % 10000 == 0:
                status('Consensus %s %s' % (self.name, pretty_number(i)))
            #if i in self.insertions:
            # Need to count absense of insertions in consensus
            #if i in self.insertions:
            #insertions = self.insertions.get(i,{}).copy()
            insertions = self.insertions[i]
            total = sum(insertions.values())
            #else:
            #    insertions = { }
            #    total = 0
            depth = self.depth[i]
            if i: depth = min(depth, self.depth[i-1])
            if depth > total: insertions['-'] = depth-total
                
            #c = consensus(insertions, min_depth,min_purity)
            c = consensus(insertions, min_depth,0.0)
            if c is not None and c != '-': 
                result.append(c)
                result_masked_only.append(c)
                report.append(('insertion-before', i, '-', c, insertions))
                
                alignment_result.append(c)
                alignment_reference.append('-' * len(c))
            
            insertion_evidence.append(insertions)
            
            
            
            #counts = { }
            #for base in self.base_counts:
            #    if self.base_counts[base][i]: 
            #        counts[base] = int(self.base_counts[base][i]) #Typecast to strip weird slow printing numpy type
            counts = self.base_counts[i]
            
            if use_ambiguity_codes:
                c = ambiguity_code_consensus(counts, min_depth,min_purity)
            else:
                c = consensus(counts, min_depth,min_purity)
            
            if c is None:
                result.append('N')
                result_masked_only.append( self.reference[i].lower() )
                has_consensus.append(False)
            elif c == '-':
                report.append(('deletion',i, self.reference[i], '-', counts))
                has_consensus.append(True)
            else:
                result.append(c)
                result_masked_only.append(c)
                has_consensus.append(c in 'ACGT') #Exclude ambiguity codes
                
                #if c != self.reference[i]:
                if self.reference[i] not in AMBIGUITY_DECODE[c]:
                    report.append(('substitution', i, self.reference[i], c, counts))
                    
            substitution_evidence.append(counts)
            alignment_result.append(c or 'N')
            alignment_reference.append(self.reference[i])
            
        
        status('')
        
        return (
            ''.join(result),
            ''.join(result_masked_only), 
            report, 
            has_consensus, 
            insertion_evidence, 
            substitution_evidence,
            ''.join(alignment_reference), 
            ''.join(alignment_result)
        )



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


def get_flag(argv, flag):
    argv = argv[:]
    any = False
    while True:
        try:
            location = argv.index(flag)
            any = True
        except ValueError: #Not found
            break
        
        del argv[location]
    return any, argv


def display_table(table, left_pad, alignment, output_file):
    if not table: return
    
    widths = [ 0 ] * len(table[0])
    for line in table:
        for i, item in enumerate(line):
            widths[i] = max(widths[i], len(item))

    for line in table:
        out_line = left_pad 
        for i, item in enumerate(line):
            pad = ' '*(widths[i]-len(item))
            if alignment[i] == 'L':
                item = item + pad
            else:
                assert alignment[i] == 'R'
                item = pad + item
            out_line += item + '  '
        output_file.write( out_line.rstrip(' ') + '\n' )
    

def main(args):
    suffix1 = '/1'
    suffix2 = '/2'
    same_direction = 0
    
    flag_solid, args = get_flag(args,'--solid')
    if flag_solid:
        suffix1 = '_R3'
        suffix2 = '_F3'
        same_direction = 1

    min_depth, args = get_option_value(args, '--depth', int, 2)
    min_purity, args = get_option_value(args, '--purity', float, 0.5)
    trim, args = get_option_value(args,'--trim', int, 0)
    infidelity, args = get_option_value(args,'--infidelity', float, 1.0)
    use_ambiguity_codes, args = get_option_value(args,'--ambiguity-codes', int, 1)
    
    max_pair_sep, args = get_option_value(args, '--max-pair-sep', int, None)
    suffix1, args = get_option_value(args, '--suffix1', str, suffix1)
    suffix2, args = get_option_value(args, '--suffix2', str, suffix2)
    same_direction, args = get_option_value(args, '--same-dir', int, same_direction)

    sys.stderr.write( 'Options:\n\n' )
    
    if max_pair_sep is None:
         max_pair_sep_text = 'not given'
    else:
         max_pair_sep_text = '%d' % max_pair_sep
    table = [
        ['--depth', '%d' % min_depth, '(minimum depth required for consensus)'],
        ['--purity', '%.2f' % min_purity, '(minimum purity required for consensus'],
        ['','',' Note: does not apply to insertions)'],
        ['--trim', '%d' % trim, '(amount to trim from start/end of alignments)'],
        ['--infidelity', '%.2f' % infidelity, '(for a hit to a read/read-pair to be declared '],
        ['','',                               ' unambiguous and counted toward consensus,'],
        ['','',                               ' any runner-up hits must score less than this'],
        ['','',                               ' proportion of the best hit\'s score)'],
        ['--ambiguity-codes', '%d' % use_ambiguity_codes, '(use IUPAC ambiguity codes.'],
        ['','',                                     ' 0-no just use Ns, 1-yes)'],
        ['','',''],
        ['--max-pair-sep', max_pair_sep_text, '(maximum distance between paired ends,'],
        ['', '', ' will treat reads as unpaired if not given)'],
        ['--suffix1', suffix1, '(suffix for first read in read pairs)'],
        ['--suffix2', suffix2, '(suffix for second read in read pairs)'],
        ['--same-dir',  '%d' % same_direction , '(read pairs have the same orientation'],
        ['','',                                 ' 0-no 1-yes)'],
    ]
    display_table(table, '  ', 'LLL', sys.stderr)

    sys.stderr.write( '\n' )
    
    if len(args) != 2:
        sys.stderr.write( USAGE )
        return 1
    
    output_dir = args[1]
    #reference_filename = args[2]
    #shrimp_filenames = args[3:]
    
    #i = 2
    #reference_filenames =  [ ]
    #while i < len(args) and args[i] != '--hits':
    #    reference_filenames.append(args[i])
    #    i += 1
    
    #shrimp_filenames = [ ]
    #if i < len(args):
    #    shrimp_filenames = args[i+1:]
    
    #if not reference_filenames:
    reference_filename = os.path.join(output_dir,'reference.fa')
        
    #if not shrimp_filenames:
    #    #shrimp_filenames = [ os.path.join(output_dir,'shrimp_hits.txt.gz') ]
    #    shrimp_filenames = [ os.path.join(output_dir,'shrimp_hits.txt') ]
    
    shrimp_filename = os.path.join(output_dir,'shrimp_hits.txt')
    
    if not os.path.exists(shrimp_filename) and \
       os.path.exists(os.path.join(output_dir,'shrimp_hits.txt.gz')):
        sys.stderr.write('shrimp_hits.txt.gz is no longer compressed,\n')
        sys.stderr.write('  you need to gunzip it or re-run shrimp_run.py\n')
        return 1
        
    for filename in [reference_filename, shrimp_filename]:
        assert os.path.exists(filename), filename + ' does not exist'
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    seqset = Refseqset()
    
    seq_order = [ ]
    for record in SeqIO.parse(open(reference_filename,'rU'), 'fasta'):
        seqset.add_sequence(record.id, record.seq.tostring().upper())
        seq_order.append(record.id)

    seqset.read_shrimp(shrimp_filename, max_pair_sep is None)

    if max_pair_sep is None:
        seqset.process_hits(trim, infidelity)
    else:
        stats_text, weird_pair_report = \
            seqset.process_paired_hits(
                max_pair_sep, same_direction, suffix1, suffix2, trim, infidelity)

        stats_file = open(os.path.join(output_dir, 'pair_stats.txt'), 'wb')
        stats_file.write( stats_text )
        stats_file.close()
        print stats_text

        weird_pair_file = open(os.path.join(output_dir, 'weird_pairs.txt'), 'wb')
        weird_pair_file.write('Left sequence\tLeft strand\tLeft position\tRight sequence\tRight strand\tRight position\n')
        for left_name,left_fwd,left_pos,right_name,right_fwd,right_pos in weird_pair_report:
            weird_pair_file.write('%s\t%s\t%d\t%s\t%s\t%d\n' % (
                left_name, left_fwd and '+' or '-', left_pos,
                right_name, right_fwd and '+' or '-', right_pos,
            ))
        weird_pair_file.close()
        
    consensus_file = open(os.path.join(output_dir, 'consensus.fa'), 'wb')
    consensus_masked_file = open(os.path.join(output_dir, 'consensus_masked.fa'), 'wb')
    reference_having_consensus_file = open(os.path.join(output_dir, 'reference_having_consensus.fa'), 'wb')
    report_file = open(os.path.join(output_dir, 'report.txt'), 'wb')
    report_gff_file = open(os.path.join(output_dir, 'report.gff'), 'wb')
    alignment_file = open(os.path.join(output_dir, 'alignment.maf'), 'wb')
    
    report_file.write( 'Sequence\tPosition in reference\tChange type\tOld\tNew\tEvidence\n' )
    
    alignment_file.write( '##maf version=1\n' )
    alignment_file.write( '#shrimp_consensus.py %s\n' % VERSION )
    
    for name in seq_order:
        consensus, consensus_masked_only, report, has_consensus, insertion_evidence, \
        substitution_evidence, alignment_reference, alignment_result = \
            seqset.seqs[name].consensus(min_depth,min_purity,use_ambiguity_codes)
        
        status('Write results for ' + name)
            
        seq = Seq.Seq(consensus, Alphabet.IUPAC.ambiguous_dna)
        record = SeqRecord.SeqRecord( seq, id=name, description='' )
        SeqIO.write([record], consensus_file, 'fasta')
        
        seq = Seq.Seq(consensus_masked_only, Alphabet.IUPAC.ambiguous_dna)
        record = SeqRecord.SeqRecord( seq, id=name, description='' )
        SeqIO.write([record], consensus_masked_file, 'fasta')
        
        #seq = Seq.Seq(''.join([ item and '1' or '0' for item in has_consensus ]))
        #record = SeqRecord.SeqRecord( seq, id=name, description='' )
        #SeqIO.write([record], has_consensus_file, 'fasta')

        ref = seqset.seqs[name].reference
        seq = Seq.Seq(''.join([ item and ref[i] or 'n' for i, item in enumerate(has_consensus) ]))
        record = SeqRecord.SeqRecord( seq, id=name, description='' )
        SeqIO.write([record], reference_having_consensus_file, 'fasta')
        
        alignment_file.write( '\n' )
        alignment_file.write( 'a\n')
        alignment_file.write( 's %s           0 %d + %d %s\n' % (name, len(ref), len(ref), alignment_reference))
        alignment_file.write( 's %s-consensus 0 %d + %d %s\n' % (name, len(ref), len(ref), alignment_result))
                
        status('Write report for ' + name)
            
        for change_type, position, old, new, counts in report:
            report_file.write( '%s\t%d\t%s\t%s\t%s\t%s\n' % (
                name,
                position+1,
                change_type,
                old,
                new,
                pretty_evidence(counts)
            ))
            
            if change_type == 'deletion':
                start = position+1
                end = position+1
                product = 'Base deleted: %s' % (old)
            elif change_type == 'insertion-before':
                #Bracket insertion
                start = position
                end = position+1
                product = 'Insertion: .' + new + '.'
            else:
                start = position+1
                end = position+1
                product = 'Substitution: %s became %s' % (old,new)
            
            product += ' ('+pretty_evidence(counts)+')'
            
            report_gff_file.write( '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n' % (
                name,
                'shrimp_consensus',
                'variation',
                start,
                end,
                '.', #score
                '+', #strand
                '.', #frame
                'product='+product
            ))
            
        status('Write evidence for ' + name)
            
        f = open(os.path.join(output_dir, filesystem_friendly_name(name) + '-evidence.txt'),'wb')
        f.write( 'Position\tInsertion-before evidence\tSubstitution evidence\tReference\n' )
        for i in xrange(len(insertion_evidence)):
            f.write( '%d\t%s\t%s\t%s\n' % (
                i+1,
                pretty_evidence(insertion_evidence[i]),
                pretty_evidence(substitution_evidence[i]),
                seqset.seqs[name].reference[i]
            ))
        
        status('Write unambiguous depth for ' + name)
            
        f = open(os.path.join(output_dir, filesystem_friendly_name(name) + '-unambiguous-depth.userplot'),'wb')
        for x in seqset.seqs[name].depth:
            f.write( '%d\n' % int(x) )
            #Note: typecast necessary for speed, numpy array printing is slow
        f.close()
        
        status('Write ambiguous depth for ' + name)
            
        f = open(os.path.join(output_dir, filesystem_friendly_name(name) + '-ambiguous-depth.userplot'),'wb')
        for x in seqset.seqs[name].depth_ambiguous:
            f.write( '%.1f\n' % float(x) )
            #Note: typecast necessary for speed, numpy array printing is slow
        f.close()

        #f = open(os.path.join(output_dir, filesystem_friendly_name(name) + '-partial-alignment.userplot'),'wb')
        #for x in seqset.seqs[name].incomplete_ends_count:
        #    f.write( '%d\n' % x ) 
        #f.close()
        #
        #f = open(os.path.join(output_dir, filesystem_friendly_name(name) + '-bias.userplot'),'wb')
        #for x in seqset.seqs[name].bias / numpy.maximum(1.0,numpy.sqrt(seqset.seqs[name].depth)):
        #    f.write( '%.1f\n' % x )
        #f.close()
        #
        #f = open(os.path.join(output_dir, filesystem_friendly_name(name) + '-entropy.userplot'),'wb')
        #for i in xrange(len(substitution_evidence)):
        #    f.write( '%.1f\n' % (
        #        entropy(substitution_evidence[i]) +
        #        entropy(insertion_evidence[i])
        #    ))
        #f.close()
        
        status('')

    consensus_file.close()
    consensus_masked_file.close()
    reference_having_consensus_file.close()
    report_file.close()
    report_gff_file.close()
    alignment_file.close()

    return 0
    
if __name__ == '__main__':
    import gc
    gc.disable()

    # Attempt Cython build    
    try:
        failure = 'not present'
        import pyximport
        failure = 'couldn\'t build'
        module = pyximport.PyxLoader('shrimp_consensus', __file__).load_module('shrimp_consensus')
        main = module.main
    except:
        sys.stderr.write('\nCython %s, falling back to Python\n\n' % failure)
        #import traceback
        #traceback.print_exc()
    
    sys.exit( main(sys.argv) )


