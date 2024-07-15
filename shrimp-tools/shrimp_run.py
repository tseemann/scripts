#!/usr/bin/python

"""

Run shrimp, set up directory for shrimp_consensus

0.2 - added option to support SOLiD reads

0.4 - autodetect FASTA/FASTQ
      output file of unmapped reads

0.5 - no longer compress,
      so shrimp_consensus can seek in hit file

"""

VERSION = '0.5'

import sys, os, pprint, gzip
from Bio import SeqIO

def how_many_cpus():
    """Detects the number of effective CPUs in the system,
    
       Function nicked from Parallel Python."""
    #for Linux, Unix and MacOS
    try:
        if hasattr(os, "sysconf"):
            if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"): 
                #Linux and Unix
                ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
                if isinstance(ncpus, int) and ncpus > 0:
                    return ncpus
            else: 
                #MacOS X
                return int(os.popen2("sysctl -n hw.ncpu")[1].read())
        #for Windows
        if os.environ.has_key("NUMBER_OF_PROCESSORS"):
            ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
            if ncpus > 0:
                return ncpus
    except:
        print >> sys.stderr, 'Attempt to determine number of CPUs failed, defaulting to 1'
    #return the default value
    return 1


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


def open_possibly_compressed_file(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename,'rb')
    else:
        return open(filename,'rb')


def read_solid(filename):
    reads_file = open(filename,'rU')
    
    while True:
        line1 = reads_file.readline()
        while line1.startswith('#'):
            line1 = reads_file.readline()
        if not line1: break
        assert line1.startswith('>'), 'Not a SOLiD CSFASTA file?'
        line2 = reads_file.readline()

        read_name = line1.rstrip('\n')[1:]
        read_seq = line2.rstrip('\n')
        yield read_name, read_seq

def read_illumina(filename):
    reads_file = open(filename,'rU')
    
    while True:
        line1 = reads_file.readline()
        if not line1: break
        line2 = reads_file.readline()
        line3 = reads_file.readline()
        line4 = reads_file.readline()
            
        assert line1.startswith('@'), 'Not an Illumina FASTQ file?'
        assert line3.startswith('+'), 'Not an Illumina FASTQ file?'
            
        read_name = line1.rstrip('\n')[1:]
        read_seq = line2.rstrip('\n')
        yield read_name, read_seq

def read_fasta(filename):
    reads_file = open(filename,'rU')
    
    line = reads_file.readline()
    while line:
        line = line.rstrip()
        assert line.startswith('>'), 'Not a FASTA file?'
        read_name = line[1:].split()[0]
        
        line = reads_file.readline()
        parts = [ ]
        while line and not line.startswith('>'):
            parts.append(line.rstrip())
            line = reads_file.readline()
        
        yield read_name, ''.join(parts)
        
def read_fasta_or_illumina(filename):
    reads_file = open(filename,'rU')    
    line = reads_file.readline()
    if line.startswith('>'):
        return read_fasta(filename)
    else:
        return read_illumina(filename)
        
   
    

argv = sys.argv[1:]

n_cpus = how_many_cpus()


#illumina, argv = get_flag(argv, '--illumina')
#fasta, argv = get_flag(argv, '--fasta')
solid, argv = get_flag(argv, '--solid')
verbose, argv = get_flag(argv, '--verbose')

limit, argv = get_option_value(argv, '--limit', int, None)
max_shrimps, argv = get_option_value(argv, '--cpus', int, n_cpus)
batch_size, argv = get_option_value(argv, '--batch-size', int, 30000000)

#Default to Illumina reads
#if not illumina and not fasta and not solid:
#    illumina = True

#assert illumina+fasta+solid == 1, 'One format only, please'


if len(argv) < 2:
    print >> sys.stderr, """\
Usage:

    shrimp_run.py output_directory [options] \\
        reference.fa [...] \\
        --reads s_X_1_sequence.txt [...] \\
        [--shrimp-options ...options to pass directly to rmapper-xx... ]

Options:

  --solid         - Reads are from a SOLiD sequencer, ie colorspace
                    (Default is FASTA or FASTQ format)
  
  --verbose       - Show output from SHRiMP
  --cpus N        - How many SHRiMPs to run in parallel 
                    (default: the number of CPUs in your system, %d)
  --batch-size N  - How many bases worth of reads to use in each SHRiMP invocation 
                    (default: 30000000)
  
  --limit N       - Only read N reads from each file

Quality information will not be used.

You can supply as many reference and read files as you like.

Files may be gzipped.
""" % n_cpus
    sys.exit(0)

output_dir = argv[0]
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

input_reference_filenames = [ ]
reads_filenames = [ ]

shrimp_options = [ ]

mode = 'reference'
for item in argv[1:]:
    if item == '--reads':
        mode = 'reads'
    elif item == '--shrimp-options':
        mode = 'shrimp_options'
    elif mode == 'reference':
        input_reference_filenames.append(os.path.abspath(item))
    elif mode == 'reads':
        reads_filenames.append(os.path.abspath(item))
    else:
        shrimp_options.append(item)

assert input_reference_filenames, 'No reference files given'
assert reads_filenames, 'No read files given'



#shrimp = '~/tmp/SHRiMP_1_1_0/bin/rmapper-ls'
#shrimp = 'rmapper-ls'
if solid:
    shrimp = 'rmapper-cs'
else:
    shrimp = 'rmapper-ls'


reference_filename = os.path.join(output_dir,'reference.fa')
reference_file = open(reference_filename,'wb')
for input_reference_filename in input_reference_filenames:
    f = open_possibly_compressed_file(input_reference_filename)
    reader = SeqIO.parse(f,'fasta')
    for record in reader:
        record.description = ''
        SeqIO.write([record], reference_file, 'fasta')
reference_file.close()

config = {
    'reads' : reads_filenames,
    'solid': solid
}
config_file = open(os.path.join(output_dir, 'config.txt'), 'wb')
pprint.pprint(config, config_file)
config_file.close()

#output_file = gzip.open(os.path.join(output_dir, 'shrimp_hits.txt.gz'), 'wb')
output_file = open(os.path.join(output_dir, 'shrimp_hits.txt'), 'wb')
unmapped_file = open(os.path.join(output_dir, 'unmapped.fa'), 'wb')

N = 0
def do_shrimp(read_set):
    global N
    my_number = N
    N += 1
    tempname = os.path.join(output_dir,'temp%d-%d.fa' % (os.getpid(),my_number))
    tempname_out = os.path.join(output_dir,'temp%d-%d.txt' % (os.getpid(),my_number))
    f = open(tempname,'wb')
    for read_name, read_seq in read_set:
        print >> f, '>' + read_name
        print >> f, read_seq
    f.close()

    command = shrimp + ' ' + ' '.join(shrimp_options) + ' ' + \
              tempname + ' ' + reference_filename + ' >' + tempname_out
    if not verbose:
        command += ' 2>/dev/null'
    #f = os.popen(command, 'r')
    child_pid = os.spawnl(os.P_NOWAIT,'/bin/sh','/bin/sh','-c',command)
    print 'SHRiMP %d running' % my_number
    
    def finalize():
        os.waitpid(child_pid, 0)
        
        reads_seen = { }
        
        for line in open(tempname_out,'rb'):
            if line.startswith('>'):
                read_name = line.split(None,1)[0][1:]
                reads_seen[read_name] = True
            output_file.write(line)
        f.close()
        output_file.flush()
        
        for read_name, read_seq in read_set:
            if read_name not in reads_seen:
                print >> unmapped_file, '>' + read_name
                print >> unmapped_file, read_seq
        unmapped_file.flush()

        os.unlink(tempname)
        os.unlink(tempname_out)
        print 'SHRiMP %d finished' % my_number
    return finalize


shrimps = [ ]
for reads_filename in reads_filenames:
    #if illumina:
    #    reader = read_illumina(reads_filename)
    #elif fasta:
    #    reader = read_fasta(reads_filename)
    #else:
    #    reader = read_solid(reads_filename)
    
    if solid:
        reader = read_solid(reads_filename)
    else:
        reader = read_fasta_or_illumina(reads_filename)
    
    read_count = 0

    while True:
        read_set = [ ]
        read_set_bases = 0
        for read_name, read_seq in reader:
            if limit != None and read_count >= limit: break
              
            read_set.append((read_name, read_seq))
            read_count += 1
            read_set_bases += len(read_seq)
            if read_set_bases >= batch_size: break
        if not read_set: break
    
        if len(shrimps) >= max_shrimps:
            shrimps.pop(0)()
        shrimps.append( do_shrimp(read_set) )

while shrimps:
    shrimps.pop(0)()

output_file.close()
unmapped_file.close()
