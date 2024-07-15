#!/usr/bin/env python

"""

Determine consequences of SNPs for amino acid coding

0.2 - Use BioPython CodonTable start codons
      Open non-gzipped genbank files

0.3 - Allow a default transl_table to be specified
      Don't die if /translation absent
      /codon_start respected (badly, but should fail on error assuming /translation present)
      Look for new end of protein if end codon changed
      Option to be --verbose or not      

0.4 - Issue a warning rather than dying on invalid start codon, lack of end codon,
      length not a multiple of 3.

0.5 - Spreadsheet friendly output: use tabs and --tabular option

0.6 - --use-coverage option

0.7 - Unicode graphlets were way too cute, removed.

0.8 - Warn only about translation mismatch
      Warnings sent to stderr

"""

VERSION = '0.8'

from Bio import Seq, SeqIO
from Bio.Data import CodonTable

import gzip, string, bisect, sys, os, numpy

class Error(Exception): pass

def warn(message):
    sys.stderr.write('* Warning: %s\n' % message)

def filesystem_friendly_name(name):
    """ Remove special characters from a name """

    for char in '\'"<>&|/\\_ .':
        name = name.replace(char,'_')
    return name

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



COMPLEMENTER = string.maketrans('ABCDGHKMRSTVWY','TVGHCDMKYSABWR')
def reverse_complement(seq):
    return seq.translate(COMPLEMENTER)[::-1]

def sequence_slice(seq, start,end):
    if start > end:
        return reverse_complement(sequence_slice(seq,end,start))
    
    start = min(len(seq),max(0,start))
    end = min(len(seq),max(0,end))
    return seq[start:end]


class Half_alignment:
    def __init__(self, text):
        self.data = [0]
        
        for char in text:
            if char == '-':
                self.data.append(self.data[-1])
            else:
                self.data.append(self.data[-1] + 1) 
                
    def __getitem__(self,key):
        if key < 0:
            return key
        if key >= len(self.data):
            return self.data[-1]+key-(len(self.data)-1)
        return self.data[key]

    def find(self, value, left=True):
        if value <= 0:
            return value
        last_value = self.data[-1]
        if value >= last_value:
            return value-last_value+(len(self.data)-1)
        
        if left:
            return bisect.bisect_left(self.data, value)
        else:
            return bisect.bisect_right(self.data, value)-1

class Half_alignment_simple:
    def __getitem__(self, key): return key
    def find(self, value, left=True): return value
half_alignment_simple = Half_alignment_simple()


class Alignment:
    def __init__(self, start1, end1, forward1, text1, start2, end2, forward2, text2):
        self.start1 = start1
        self.end1 = end1
        self.forward1 = forward1
        if text1 is None:
            self.ali1 = half_alignment_simple
        else:
            self.ali1 = Half_alignment(text1)
        
        self.start2 = start2
        self.end2 = end2
        self.forward2 = forward2
        if text2 is None:
            self.ali2 = half_alignment_simple
        else:
            self.ali2 = Half_alignment(text2)

    def project(self, position, left=True):
        """ Positions lie between bases.
            start and end are positions. """
            
        if self.forward1:
            position = position - self.start1
        else:
            position = self.end1 - position
        
        position = self.ali2[ self.ali1.find(position, left) ]
        
        if self.forward2:
            position = position + self.start2
        else:
            position = self.end2 - position
        
        return position

    def back_project(self, position, left=True):
        """ Positions lie between bases.
            start and end are positions. """
            
        if self.forward2:
            position = position - self.start2
        else:
            position = self.end2 - position
        
        position = self.ali1[ self.ali2.find(position,left) ]
        
        if self.forward1:
            position = position + self.start1
        else:
            position = self.end1 - position
        
        return position


def open_possibly_compressed_file(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename,'rb')
    else:
        return open(filename,'rb')




def get_dna(sequence, feature):
    if not feature.sub_features:
        return record.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]
    
    parts = [ get_dna(sequence, sub_feature) for sub_feature in feature.sub_features ]
    result = parts[0]
    for item in parts[1:]: result = result + item
    return result 


def half_alignment_from_feature(feature):
    if not feature.sub_features:
        xs = 'X' * (feature.location.nofuzzy_end-feature.location.nofuzzy_start)
        return xs
    
    result = ''
    position = feature.location.nofuzzy_start
    
    for sub_feature in feature.sub_features:
        assert sub_feature.strand == feature.strand 

        step = sub_feature.location.nofuzzy_start - position
        assert step >= 0
        result += '-' * step
        
        result += half_alignment_from_feature(sub_feature)
        
        position = sub_feature.location.nofuzzy_end
    
    return result
        
def alignment_from_feature(seq, feature):
    hal2 = half_alignment_from_feature(feature)
    hal1 = 'X' * len(hal2)
    assert len(hal2) == feature.location.nofuzzy_end-feature.location.nofuzzy_start
    
    forward = feature.strand != -1
    if not forward: hal2 = hal2[::-1]
    return Alignment(feature.location.nofuzzy_start,feature.location.nofuzzy_end,forward,hal1,
                     0, hal2.count('X'),True, hal2)

class Half: pass

def load_alignments(filename):
    f = open(filename, 'rU')
    line = f.readline()
    
    result = [ ]
    while True:        
        while line and (line.startswith('#') or not line.strip()):
            line = f.readline()
        if not line: break
        
        assert line.startswith('a')
        
        line = f.readline()
        lines = [ ]
        while line.startswith('s'):
            lines.append(line)
            line = f.readline()
        
        assert len(lines) == 2, 'Only expecting two sequences in alignment'

        halves = [ ]
        for item in lines:
            half = Half()
            halves.append(half)
                    
            parts = item.rstrip().split()
            half.name = parts[1]
            half.ali = parts[6]
            half.seq = half.ali.replace('-','')
            assert int(parts[2]) == 0
            assert parts[4] == '+'
            assert int(parts[3]) == int(parts[5])
        
        result.append((
            halves[0].name,
            halves[0].seq, 
            halves[1].seq,
            Alignment(
                0,len(halves[0].seq),True, halves[0].ali,
                0,len(halves[1].seq),True, halves[1].ali
            )
        )) 
    
    f.close()
    
    return result    
        


def get_graph(path, name, suffix):
    filename = os.path.join(path, filesystem_friendly_name(name) + '-' + suffix + '.userplot')
    
    result = [ ]
    for item in open(filename,'rb'):
        result.append( float(item.strip()) )
    return numpy.array(result)


def graphlet(data):
    n = min(len(data),10)
    result = ''
    for i in xrange(n):
        value = numpy.average(data[i*len(data)//n:(i+1)*len(data)//n])
        
        if value <= 0.0:
            result += '_'
        else:
            value = int(value+0.5)
            if value > 9:
                result += '!'
            else:
                result += str(value) 
        
        #x = int(value*3.5)
        #
        #if value <= 0.0:
        #    result += unichr(0x2594)
        #elif value >= 2.0:
        #    result += unichr(0x2592)
        #else:
        #    result += unichr(0x2581+max(0,min(7,x)))
    return result #( unichr(0x2595) + result + unichr(0x258f) ).encode('utf-8') 
        
        

USAGE = """
USAGE:

    consequences.py [options] genbank_file.gbk[.gz] working_dir/alignments.maf
    
Options:

    --use-coverage       - Use depth of coverage information.
                           Produces small graphs of coverage depth for each CDS,
                           with areas greater than 2 x the median highlighted.
                           Note: this will produce output for every CDS

    --transl_table NN    - Translation table to use, 
                           will be overridden by /transl_table qualifier if 
                           present. If not specified, /transl_table qualifier 
                           in each CDS is mandatory.
                           "--transl_table 1" is safe for NCBI GenBank files.
    
    --tabular            - Spreadsheet friendly output.
    
    --noheader           - No title row.
    
    --verbose            - List amino acid changes, rather than just counting 
                           them.
        
"""

def main(argv):
    default_transl_table, argv = get_option_value(argv, '--transl_table', int, None)
    use_coverage, argv = get_flag(argv, '--use-coverage')
    tabular, argv = get_flag(argv, '--tabular')
    noheader, argv = get_flag(argv, '--noheader')
    verbose, argv = get_flag(argv, '--verbose')

    if len(argv) != 3:
        print USAGE
        return 1
    
    genbank_filename = argv[1]
    alignment_filename = argv[2]
    
    working_dir = os.path.split(alignment_filename)[0]
    
    alignments = load_alignments(alignment_filename)
    
    summaries = [ ]
    details = [ ]
    
    if not noheader:
        fields = 'Sequence\tLocus tag\tOld length (aa)\tNew length (aa)\tAmino acid changes\t'
        if use_coverage: fields += 'Unambiguous coverage vs median\t\tAmbiguous coverage vs median\t\tAmbiguous percent with any hits\t'
        fields += 'Gene\tProduct'
        if tabular: fields += '\tChanges of note'
        print fields
    
    for record in SeqIO.parse(open_possibly_compressed_file(genbank_filename),'genbank'):
        #print dir(record)
        #print record.id
        
        sequence = record.seq.tostring()
    
        for name, seq1, seq2, alignment in alignments:
            if seq1 != sequence: continue
            #print
            #print '=====', name, '=====' 
 
            if use_coverage:       
                depth = get_graph(working_dir, name, 'unambiguous-depth')
                ambiguous_depth = get_graph(working_dir, name, 'ambiguous-depth')
                median_depth = numpy.median(depth)
                median_ambiguous_depth = numpy.median(ambiguous_depth)
    
            for feature in record.features:
                if feature.type != 'CDS': continue
                
                locus_tag = feature.qualifiers['locus_tag'][0]
                
                if 'transl_table' in feature.qualifiers:
                    transl_table_no = int(feature.qualifiers['transl_table'][0])
                else:
                    assert default_transl_table is not None, 'No /transl_table for CDS, and default transl_table not given'
                    transl_table_no = default_transl_table
                
                transl_table = CodonTable.ambiguous_dna_by_id[transl_table_no]
                start_codons = transl_table.start_codons
                
                feature_alignment = alignment_from_feature(sequence, feature)
                
                dna = [ ]
                new_dna = [ ]
                shifts = [ ]
                for i in xrange(feature_alignment.end2):
                    p1 = feature_alignment.back_project(i, left=False)
                    p2 = feature_alignment.back_project(i+1, left=True)
                    assert abs(p2-p1) < 2
                    dna.append( sequence_slice(sequence,p1,p2) )
                    
                    p1a = alignment.project(p1, left=False)
                    p2a = alignment.project(p2, left=False) #Hmm
                    
                    diff = (p2-p1)-(p2a-p1a)
                    #if diff:
                    #    if diff%3:
                    #        frame_shift = True
                    #    else:
                    #        frame_preserving_shift = True
                    new_dna.append( sequence_slice(seq2,p1a,p2a) )
                    
                    if diff:
                        shifts.append((i,dna[-1],new_dna[-1]))
                    
                dna = ''.join(dna)
                new_dna = ''.join(new_dna)
                
                # This usually indicated a CDS truncated at the start?
                # in which case, will probably fail some way or other down the line.
                if 'codon_start' in feature.qualifiers:
                    codon_start = int(feature.qualifiers['codon_start'][0]) - 1
                else:
                    codon_start = 0
                dna = dna[codon_start:]
                new_dna = new_dna[codon_start:]
                
                if len(dna) % 3 != 0:
                    warn(locus_tag + ' length not a multiple of 3')
                #assert len(new_dna) % 3 == 0
                
                protein = Seq.Seq(dna).translate(table=transl_table_no).tostring()            
                # http://en.wikipedia.org/wiki/Start_codon is always translated to M
                protein = 'M' + protein[1:]
                
                if dna[:3] not in start_codons:
                    warn(locus_tag + ' has unknown start codon: ' + dna[:3])
                    
                if not protein.endswith('*'):
                    warn(locus_tag + ' lacks end codon')
                                
                if 'translation' in feature.qualifiers:
                    expect = feature.qualifiers['translation'][0]
                    if protein[:-1] != expect:
                        warn(locus_tag + ' translation given in feature does not match translation from DNA')                
    
                new_protein = Seq.Seq(new_dna).translate(table=transl_table_no).tostring()            
                new_protein = 'M' + new_protein[1:]

                # If end codon changed, find new end                
                # Don't bother if there are unknown amino acids
                if 'X' not in new_protein and '*' not in new_protein:
                    #This is very inefficient
                    i = feature_alignment.end2
                    while True:
                        p1 = feature_alignment.back_project(i, left=False)
                        p2 = feature_alignment.back_project(i+1, left=True)
                        p1a = alignment.project(p1, left=False)
                        p2a = alignment.project(p2, left=False) #Hmm
                        if p1a < 0 or p2a < 0 or p1a > len(seq2) or p2a > len(seq2):
                            break
                            
                        new_dna += sequence_slice(seq2,p1a,p2a)                        
                        new_protein = Seq.Seq(new_dna).translate(table=transl_table_no).tostring()            
                        new_protein = 'M' + new_protein[1:]
                        if 'X' in new_protein or '*' in new_protein: break
                        
                        i += 1
                
                if '*' in new_protein:
                    new_protein = new_protein[:new_protein.index('*')+1] 
    

                diffs = [ ]
                for i in xrange(min(len(new_protein),len(protein))):
                    if protein[i] != new_protein[i] and new_protein[i] != 'X':
                        diffs.append(i)
    
                diff_start = new_dna[0] not in ('N',dna[0]) or \
                             new_dna[1] not in ('N',dna[1]) or \
                             new_dna[2] not in ('N',dna[2]) 
    
                if use_coverage or diffs or diff_start or shifts or len(new_protein) != len(protein):
                    line = name + '\t' + locus_tag + '\t' + \
                          '%d\t' % (len(protein)-1) + \
                          '%d\t' % (len(new_protein)-1) + \
                          '%d\t' % len(diffs)
                    
                    if use_coverage:
                        cds_depth = depth[feature_alignment.start1:feature_alignment.end1] / median_depth
                        if not feature_alignment.forward1: cds_depth = cds_depth[::-1]
                        cds_ambiguous_depth = ambiguous_depth[feature_alignment.start1:feature_alignment.end1] / median_ambiguous_depth
                        if not feature_alignment.forward1: cds_ambiguous_depth = cds_ambiguous_depth[::-1]
                        #cds_average_depth_ratio = numpy.average(depth[feature_alignment.start1:feature_alignment.end1]) / median_depth 
                        #cds_average_ambiguous_depth_ratio = numpy.average(ambiguous_depth[feature_alignment.start1:feature_alignment.end1]) / median_ambiguous_depth                        
                        #line += '%.1f\t' % cds_average_depth_ratio 
                        #line += '%.1f\t' % cds_average_ambiguous_depth_ratio
                        
                        #line += '%.1f..%.1f\t' % (numpy.minimum.reduce(cds_depth)/median_depth, numpy.maximum.reduce(cds_depth)/median_depth) 
                        #line += '%.1f+/-%.1f\t' % (numpy.average(cds_depth)/median_depth, numpy.var(cds_depth)**0.5/median_depth) 
                        #line += '%.1f..%.1f\t' % (numpy.minimum.reduce(cds_ambiguous_depth)/median_ambiguous_depth, numpy.maximum.reduce(cds_ambiguous_depth)/median_ambiguous_depth)
                        
                        line += '%.1f\t' % numpy.average(cds_depth) + graphlet(cds_depth)+'\t' 
                        line += '%.1f\t' % numpy.average(cds_depth) + graphlet(cds_ambiguous_depth)+'\t'
                        line += '%.1f%%\t' % (numpy.average(cds_ambiguous_depth > 0.0)*100.0)
                    
                    line += '%s\t' % feature.qualifiers.get('gene',[''])[0] + \
                            '%s' % feature.qualifiers.get('product',[''])[0]
                    
                    notes = [ ]
                    
                    if use_coverage and 'X' in new_protein:
                        xs = new_protein.count('X')
                        if xs == len(new_protein)-1: #First is M, so len-1
                            notes.append('\ No consensus')
                        else:
                            notes.append('\ No consensus for %d aa' % (new_protein.count('X')))
                                       
                    if len(new_protein) < len(protein):
                        notes.append('\ Shorter by %d aa' % (len(protein)-len(new_protein)))

                    if len(new_protein) > len(protein):
                        notes.append('\ Longer by %d aa' % (len(new_protein)-len(protein)))
                    
                    if diff_start:
                        notes.append('\ Start changed: %s -> %s' % (dna[:3], new_dna[:3]))
                        if new_dna[:3] not in start_codons:
                            notes.append('  No longer a start codon!')
                            
                    if shifts:
                        notes.append('\ Indels:')
                    
                        for pos, old, new in shifts:
                            notes.append('    base %5d / codon %5d   %s -> %s' % (pos+1,(pos//3)+1,old,new or '-'))
                        
                    if diffs:
                        if verbose:
                            notes.append('\ Amino acid changes:')
                            for i in diffs:
                                notes.append('    codon %5d   %s->%s   (%s->%s)' % (i+1, protein[i], new_protein[i], dna[i*3:i*3+3], new_dna[i*3:i*3+3]))
                    
                    #if len(new_protein) > len(protein):
                    #    print 'New protein is longer:', new_protein[len(protein):]
                    #if len(new_protein) < len(protein):
                    #    print 'New protein is shorter:', protein[len(new_protein):]
                    #print protein
                    #print new_protein
                    
                    if tabular:
                        print line + '\t' + ' '.join([ ' '.join(note.strip().split()) for note in notes ])
                    else:
                        print line
                        for note in notes:
                            print '\t' + note
    return 0
    
if __name__ == '__main__':
    sys.exit( main( sys.argv ) )


