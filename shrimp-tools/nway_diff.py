#!/usr/bin/env python2.6

import sys, os

from Bio import SeqIO

USAGE = """\

Usage:

    nway_diff.py [options] workingdir1 workingdir2 [...]

Options:

    --indels 0/1       count insertions and deletions, 0=no, 1=yes (default)

    --reference 0/1    include reference, 0=no, 1=yes (default) 
    
    --list 0/1         list differences rather than just counting them
                       0=no (default), 1=yes

Lack of consensus or ambiguity codes in any data set will cause differences 
at that position to be ignored.

Three or more different bases/insertions at a given position will cause it to be
ignored.

"""


AMBIGUOUS = list('KMRYSWBVHD')


class Error(Exception): 
    pass

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





def main(argv):
    use_indels, argv = get_option_value(argv,'--indels',int,1)
    use_reference, argv = get_option_value(argv,'--reference',int,1)
    make_list, argv = get_option_value(argv,'--list',int,0)

    if len(argv) < 2:
        sys.stderr.write(USAGE)
        return 1
    
    working_dirs = argv[1:]
    
    reference_data = { } # (ref_name, position, change_type) -> string
    strain_data = { } # working_dir -> (ref_name, position, change_type) -> string
    
    strain_reference_having_consensus = { } # working_dir -> ref_name -> string
    
    for working_dir in working_dirs:
        assert working_dir not in strain_data, 'Working directory given twice'
        strain_data[working_dir] = { }
        
        report_file = open(os.path.join(working_dir, 'report.txt'), 'rU')
        report_file.readline()
        for line in report_file:
            ref_name, position, change_type, old, new, evidence = \
                line.rstrip('\n').split('\t')
            
            if change_type == 'deletion':
                change_type = 'substitution'
            
            if not use_indels and \
               (change_type == 'insertion-before' or new == '-'):
                continue
            
            key = (ref_name, int(position), change_type)
            if key in reference_data:
                assert reference_data[key] == old
            else:
                reference_data[key] = old
            
            strain_data[working_dir][key] = new
        report_file.close()
        
        strain_reference_having_consensus[working_dir] = { }
        ref_have_con_file = open(os.path.join(working_dir, 'reference_having_consensus.fa'), 'rU')
        for record in SeqIO.parse(ref_have_con_file, 'fasta'):
            strain_reference_having_consensus[working_dir][record.id] = \
                record.seq.tostring().upper()
        ref_have_con_file.close()

    keys = sorted(reference_data)

    #Fill in any blanks
    for working_dir in working_dirs:
        for key in keys:
            if key in strain_data[working_dir]: continue
        
            # - Positions in report files start from 1 not 0
            # - Insertions must be bracketed
            lacks_consensus = (
                strain_reference_having_consensus[working_dir][key[0]][key[1]-1] == 'N' or
                (key[2] == 'insertion-before' and key[1] > 1 and
                 strain_reference_having_consensus[working_dir][key[0]][key[1]-2] == 'N')
            )
            
            #If there's no consensus, record it as ambiguous
            if lacks_consensus:
                strain_data[working_dir][key] = 'N'                
            else:
                strain_data[working_dir][key] = reference_data[key]

 
    all_data_names = ([ 'reference' ] if use_reference else []) + working_dirs
    all_data = ([ reference_data ] if use_reference else []) + \
               [ strain_data[working_dir] for working_dir in working_dirs ] 
    
    
    ones = ( 1 << len(all_data) )-1
    
    total_differences = 0
    
    if make_list:
        print '\t'.join(['Partition','Sequence','Position in reference','Change type'] + all_data_names) 
    
    for i in xrange(1,(1<<len(all_data))-1,2):
        set1 = [ ]
        set2 = [ ]
        for j in xrange(len(all_data)):
            if i & (1<<j):
                set1.append(j)
            else:
                set2.append(j)

        if make_list:
            print
            print ', '.join( all_data_names[i] for i in set1 ) + '   vs   ' + \
                  ', '.join( all_data_names[i] for i in set2 )
            print
                
        n = 0
        for key in keys:
            values = [ item[key] for item in all_data ]
            #print key, values
            
            #Skip where no difference or three or more things seen
            if len(set(values)) != 2:
                continue
            
            #Skip if *any* ambiguity
            if any( value in AMBIGUOUS for value in values ):
                continue
            
            if any( values[i] != values[set1[0]] for i in set1[1:] ) or \
               any( values[i] != values[set2[0]] for i in set2[1:] ):
                continue
            
            #print key[0], key[1], key[2], values[set1[0]] + ' vs ' + values[set2[0]]
            
            if make_list:
                change_type = key[2]
                if change_type == 'substitution' and '-' in values: change_type = 'deletion'
                print '\t%s\t%s\t%s\t' % (key[0],key[1],change_type) + '\t'.join(values) 
            
            n += 1

        total_differences += n

        if not make_list:
            print ', '.join( all_data_names[i] for i in set1 ) + '   vs   ' + \
                  ', '.join( all_data_names[i] for i in set2 ) + \
                  ': %d differences' %n            

    if not make_list:
        print
        print 'Total: %d' % total_differences


    if make_list:
        print
        print 'Ignored'
        print
    
    n_multiway = 0
    n_ambiguous = 0    
    for key in keys:
        values = [ item[key] for item in all_data ]
        #print key, values
        
        confusing = False
        if any( value in AMBIGUOUS for value in values ):
            n_ambiguous += 1
            confusing = True
        elif len(set(values)) > 2:
            n_multiway += 1
            confusing = True
        
        if make_list and confusing:
            print '\t%s\t%s\t%s\t' % (key[0],key[1],key[2]) + '\t'.join(values)

    if not make_list:
        print
        print 'Ambiguities ignored: %d' % n_ambiguous
        print 'Multi-way changes ignored: %d' % n_multiway
    
    return 0

if __name__ == '__main__':
    sys.exit( main(sys.argv) )

