#!/usr/bin/env python

# ==============================================================================
# PanPhlAn v1.3: PANgenome-based PHyLogenomic ANalysis
# Detecting and characterizing strains in metagenomic samples
#
# Matthias Scholz, Doyle V. Ward, Edoardo Pasolli, Thomas Tolio, Moreno Zolfo,
# Francesco Asnicar, Duy Tin Truong, Adrian Tett, Ardythe L. Morrow, and Nicola Segata.
# Strain-level microbial epidemiology and population genomics from shotgun metagenomics.
# Nature Methods, 13, 435-438, 2016.
#
# For help type "./panphlan_map.py -h"
#
# https://bitbucket.org/CibioCM/panphlan
# ==============================================================================

from argparse import ArgumentParser
from collections import defaultdict
from shutil import copyfileobj
from .utils import end_program, show_interruption_message, show_error_message, time_message, find
import bz2, fnmatch, multiprocessing, operator, os, subprocess, sys, tempfile, time
from distutils.version import LooseVersion

__author__  = 'Matthias Scholz, Thomas Tolio, Leonard Dubois and Nicola Segata (panphlan-users@googlegroups.com)'
__version__ = '1.3'
__date__    = '15 November 2019'


DEFAULT_READ_LENGTH     = 70  # min=70 to account for read-length in the MetaHIT project
PANPHLAN                = 'panphlan_'
TEMPORARY_FILE          = 't'
USER_DEFINED            = 'd'

# Paths
DEFAULT_BOWTIE2_INDEXES = 'BOWTIE2_INDEXES/'

# Pangenome CSV file constants
FAMILY_INDEX            = 0
GENE_INDEX              = 1
CONTIG_INDEX            = 3
FROM_INDEX              = 4
TO_INDEX                = 5

# Messages
WAIT_MESSAGE            = '[W] Please wait. The computation may take several minutes...'
INTERRUPTION_MESSAGE    = '[E] Execution has been manually halted.\n'

# Operating systems
LINUX                   = 'lin'
WINDOWS                 = 'win'

# File extensions
CSV                     = 'csv'
FASTQ                   = 'fastq'
BAM                     = 'bam'
SAM                     = 'sam'
GZ                      = 'gz'
BZ2                     = 'bz2'
TAR_BZ2                 = 'tar.bz2'
TAR_GZ                  = 'tar.gz'
SRA                     = 'sra'

COMPRESSED_FORMATS      = [BZ2, GZ]
ARCHIVE_FORMATS         = [TAR_BZ2, TAR_GZ, SRA]
FILE_EXTENSION          = ''

# input file endings, expected fasta or fastq format, decompression command
INPUT_FORMAT={}
INPUT_FORMAT['tar.bz2']  =('fastq',['tar', '-xOf'])
INPUT_FORMAT['tar.gz']   =('fastq',['tar', '-xOf'])
INPUT_FORMAT['sra']      =('fastq',['fastq-dump', '-Z', '--split-spot', '--minReadLen',str(DEFAULT_READ_LENGTH)])
INPUT_FORMAT['bz2']      =('fastq',['bzcat'])
INPUT_FORMAT['gz']       =('fastq',['gunzip', '-c']) # gzcat, zcat not always available
INPUT_FORMAT['fq.bz2']   =('fastq',['bzcat'])
INPUT_FORMAT['fa.bz2']   =('fasta',['bzcat'])
INPUT_FORMAT['fastq.bz2']=('fastq',['bzcat'])
INPUT_FORMAT['fasta.bz2']=('fasta',['bzcat'])
INPUT_FORMAT['fq.gz']    =('fastq',['gunzip', '-c'])
INPUT_FORMAT['fa.gz']    =('fasta',['gunzip', '-c'])
INPUT_FORMAT['fastq.gz'] =('fastq',['gunzip', '-c'])
INPUT_FORMAT['fasta.gz'] =('fasta',['gunzip', '-c'])
INPUT_FORMAT['fastq']    =('fastq',['cat'])
INPUT_FORMAT['fasta']    =('fasta',['cat'])
INPUT_FORMAT['fq']       =('fastq',['cat'])
INPUT_FORMAT['fa']       =('fasta',['cat'])
INPUT_FORMAT['bam']      =('bam',[''])
# only archive file (tar,sra) are decompressed and piped to bowtie2,
# other files are directly given as bowtie2 input, but might be piped in future

# Error codes
INEXISTENCE_ERROR_CODE      =  1 # File or folder does not exist
UNINSTALLED_ERROR_CODE      =  2 # Software is not installed
FILEFORMAT_ERROR_CODE       =  3 # Format is not acceptable
ARCHIVE_ERROR_CODE          =  4 # Some problem with the zip file (archive)
BOWTIE2_ERROR_CODE          =  5 # Some problem with Bowtie2
SAMTOOLS_ERROR_CODE         =  6 # Some problem with Samtools
INTERRUPTION_ERROR_CODE     =  7 # Computation has been manually halted
INDEXES_NOT_FOUND_ERROR     =  8 # Bowtie2 indexes are not found
PANGENOME_ERROR_CODE        =  9 # Pangenome .csv file cannot be found
PARAMETER_ERROR_CODE        = 10 #


# ------------------------------------------------------------------------------
# ARGUMENT PARSER
# ------------------------------------------------------------------------------

class PanPhlAnParser(ArgumentParser):
    '''
    Subclass of ArgumentParser for parsing command inputs for panphlan.py
    '''
    def __init__(self):
        ArgumentParser.__init__(self)
        self.add_argument('-i','--input',           metavar='INPUT_FILE',                   type=str,                                                   help='Short read input files. If no file is specified, panphlan_map reads from standard input.')
        # self.add_argument('-f', '--input_format',   metavar='INPUT_FORMAT',                 type=str,                                                   help='Old option, will be removed in future version')
        self.add_argument('--i_bowtie2_indexes',    metavar='INPUT_BOWTIE2_INDEXES',        type=str,   default=None,                                   help='Input directory of bowtie2 indexes and pangenome')
        self.add_argument('--fastx',                metavar='FASTX_FORMAT',                 type=str,   default='fastq',choices=['fastq','fasta','bam'],help='Read input format (fasta or fastq), default: fastq, if not fasta recognized by file ending.')
        self.add_argument('-c','--clade',           metavar='CLADE_NAME',                   type=str,                                   required=True,  help='Name of the species or clade: -c ecoli16')
        self.add_argument('-o','--output',          metavar='OUTPUT_FILE',                  type=str,   default='map_results/',                         help='Mapping result output-file: -o path/sampleID_clade.csv')
        self.add_argument('-b','--out_bam',         metavar='OUTPUT_BAM_FILE',              type=str,                                                   help='Get BAM output file: --out_bam sampleID.bam')
        self.add_argument('-p', '--nproc',          metavar='NUMOF_PROCESSORS',             type=int,   default=12,                                     help='Maximum number of processors to use. Default is 12 or a lower number of available processors.')
        self.add_argument('-m', '--mGB',            metavar='MEMORY_GIGABTES_FOR_SAMTOOLS', type=float, default=4.0,                                    help='Maximum amount of memory for Samtools: -mGB 4')
        self.add_argument('--readLength',           metavar='READS_LENGTH',                 type=int,   default=DEFAULT_READ_LENGTH,                    help='Minimum read length, default: --readLength 70')
        self.add_argument('--bt2',                  metavar='BOWTIE2_OPTIONS',              type=str,   default='--very-sensitive',                     help='Additional bowtie2 mapping options, separated by slash: /-D/20/-R/3/, default: -bt2 /--very-sensitive/')
        self.add_argument('--th_mismatches',        metavar='NUMOF_MISMATCHES',             type=int,   default=-1,                                     help='Number of mismatches to filter (bam)')
        self.add_argument('--tmp',                  metavar='TEMP_FOLDER',                  type=str,   default='TMP_panphlan_map',                     help='Directory for temporary files, default: --tmp TMP_panphlan_map')
        self.add_argument('--verbose',              action='store_true',                                                                                help='Display progress information')
        self.add_argument('-v', '--version',        action='version',   version="PanPhlAn version "+__version__+"\t("+__date__+")",                     help='PanPhlAn version')



def detect_input_format(filename):
    '''
    Detect input file format
    Return: file-extension, fasta/fastq format, decompress command
    Error, if unknown file format
    '''
    for extension in reversed(sorted(INPUT_FORMAT.keys(), key=len)):
        if filename.endswith(extension):
            fastx, decompress = INPUT_FORMAT[extension]
            return extension, fastx, decompress
    # Unacceptable format is found. Raise an error and close the program
    show_error_message('Input with unacceptable extension/format.')
    sys.exit(FILEFORMAT_ERROR_CODE)
# ------------------------------------------------------------------------------
def correct_output_name(opath, ipath, panphlan_clade, VERBOSE):
    '''
    complete output filename: add clade, extension, etc..  if missing
    correct filename: ERR260216_ecoli12.csv
    Possible user input:
        a) empty (only directory): use input sampleID and add _ecoli12.csv
        b) -o ERR260216 (adding: _ecoli12.csv(.bz2) )
        c) -o ERR260216_ecoli12.csv.bz2 (avoid double .bz2.bz2)
        d) -o ERR260216_ecoli12.csv (normal case, keep like it was)
        e) -o ERR.260216_ecoli12.csv (allow dot in sample-name)
    '''
    clade=panphlan_clade.replace('panphlan_','')
    orig_path = opath

    # get output directory,  including ending '/'
    outdir=os.path.join(os.path.dirname(opath), '')

    # automatic correction of output ending
    if not opath.endswith('_' + clade + '.csv'): # ending not OK
        # get out-filename without extensions (remove/ignore endings: .csv or .csv.bz2)
        outfilename = os.path.basename(opath).split('.')[0]
        if not outfilename: # empty filename
            if ipath==None: # if stdin (pipe) input
                sys.exit('Please specify output filename: -o Outputpath/sampleID')
            else: # takesample ID from input-filename
                outfilename = os.path.basename(ipath).split('.')[0]

        # add clade if not in filename
        if not outfilename.endswith(clade):
            outfilename = outfilename + '_' + clade
        # add .csv
        outfilename = outfilename + '.csv'
        # get full path
        opath = outdir + outfilename

        if not orig_path==opath:
            print('[I] Extend output filename "' + orig_path + '" to "' + opath + '"')

    # create output directory
    if outdir: # no empty dir string
        if os.path.exists(outdir):
            if VERBOSE: print('[I] Mapping result directory: ' + outdir)
        else:
            if VERBOSE: print('[I] Create mapping result directory: ' + outdir)
            os.makedirs(outdir)
    else:
        if VERBOSE: print('[I] Mapping results are saves in working directory.')

    return opath
# ------------------------------------------------------------------------------
def build_pangenome_dicts(pangenome_file, TIME, VERBOSE):
    '''
    Build the dictionary for contig -> included gene -> location of the gene in the DNA
    '''
    contig2gene = {}
    with open(pangenome_file, mode='r') as f:
        for line in f:
            # line = FAMILY, GENE, CONTIG, FROM, TO
            words = line.strip().split('\t')
            fml, gen, ctg, fr, to = words[FAMILY_INDEX], words[GENE_INDEX], words[CONTIG_INDEX], int(words[FROM_INDEX]), int(words[TO_INDEX])
            if not ctg in contig2gene:
                contig2gene[ctg] = {}
            contig2gene[ctg][gen] = (min(to, fr), max(to, fr))
    # Sort genes in each contig for deterministic results even if pangenome file contains overlaps (and so errors)
    # for ctg in contig2gene:
    #     contig2gene[ctg] = sorted(contig2gene[ctg].items(), key=operator.itemgetter(1))
    if VERBOSE:
        TIME = time_message(TIME, 'Dictionary for {contig:{gene:(from,to)}} has been created.')
    return contig2gene, TIME
# -----------------------------------------------------------------------------
def get_pangenome_file(bowtie2_indexes_dir, clade, VERBOSE):
    '''
    Get the pangenome file for the considered specie
    '''
    pangenome = find(clade + '_pangenome.csv', bowtie2_indexes_dir)
    if len(pangenome) > 1:
        message = '[W] More than one matchable pangenome is found! '
        if VERBOSE:
            message += '\n    Files are: '
            message += ', '.join(pangenome)
            message += '\n    '
        message += 'Choose file "' + pangenome[0] + '". If not good, please resolve manually the problem.'
        print(message)
    elif len(pangenome) < 1:
        sys.stderr.write('[E] Cannot find the pangenome file for ' + clade + ' in directory ' + bowtie2_indexes_dir + '\n')
        sys.exit(PANGENOME_ERROR_CODE)
    else:
        if VERBOSE: print('[I] Pangenome file: ' + pangenome[0])
    return pangenome[0]
# -----------------------------------------------------------------------------
def samtools_sam2bam(in_sam, out_bam, memory, tmp_path, TIME, VERBOSE):
    '''
    Covert a SAM file into BAM file, then sort the BAM
        1.  samtools view -bS <INPUT SAM FILE>
        2.  samtools sort
              samtools version 1.2
                samtools sort <in.bam> <out.prefix>
                cat sample.bam | samtools sort - tmp_sorted
              samtools version 1.3
                samtools sort <in.bam> -o <out.bam>
                cat sample.bam | samtools sort - -o tmp_sorted.bam
    '''
    outcome = (None, None)
    try:
        # get samtools version
        samtools, samtools_version = check_samtools()
        # 1st command: samtools view -bS <INPUT SAM FILE>
        view_cmd = ['samtools', 'view', '-bS', in_sam.name]
        if VERBOSE:
            print('[I] ' + ' '.join(view_cmd))
        p2 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
        if VERBOSE:
            print('[I] Temporary .bam file has been generated')

        try:
            tmp_bam = None
            # 2nd command: samtools sort -m <AMOUNT OF MEMORY> - <OUTPUT BAM FILE>
            sort_cmd = ['samtools', 'sort', '-m', str(int(memory * 1024*1024*1024))]

            if out_bam == None: # .bam file is not saved, only temporary bam file
                if tmp_path == None:
                    tmp_bam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.bam')
                else:
                    tmp_bam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.bam', dir=tmp_path)

                with tmp_bam:
                    if LooseVersion(samtools_version) >= LooseVersion('1.3'): # works also with two dots 1.3.1 (1 of 2 if's)
                        sort_cmd += ['-', '-o', tmp_bam.name]
                    else: # older samtools versions: only prefix, without .bam
                        sort_cmd += ['-', tmp_bam.name[:-4]]
                    if VERBOSE:
                        print('[I] cmd (v'  + samtools_version + '): ' + ' '.join(sort_cmd))
                    p3 = subprocess.Popen(sort_cmd, stdin=p2.stdout, stdout=tmp_bam)
                    p3.wait() # Wait until previous process has finished its computation (otherwise there will be error raised by Samtools)
                    if VERBOSE:
                        print('[I] Temporary .bam file ' + tmp_bam.name + ' has been sorted')
                outcome = (TEMPORARY_FILE, tmp_bam.name)

            else: # .bam file is saved (option -b Bam/sample.bam)
                if LooseVersion(samtools_version) >= LooseVersion('1.3'): # works also with two dots 1.3.1 (2 of 2 if's)
                    sort_cmd += ['-', '-o', out_bam]
                else: # older samtools versions: only prefix, without .bam
                    sort_cmd += ['-', out_bam[:-4]]
                if VERBOSE:
                    print('[I] cmd (v'  + samtools_version + '): ' + ' '.join(sort_cmd))
                with open(out_bam, mode='w') as obam:
                    p3 = subprocess.Popen(sort_cmd, stdin=p2.stdout, stdout=obam)
                    p3.wait() # Wait until previous process has finished its computation (otherwise there will be error raised by Samtools)
                if VERBOSE:
                    print('[I] User-defined .bam file ' + out_bam + ' has been sorted')
                outcome = (USER_DEFINED, out_bam)

            if VERBOSE:
                TIME = time_message(TIME,'Samtools SAM->BAM translation (view+sort) completed.')

        except (KeyboardInterrupt, SystemExit):
            p3.kill()
            show_interruption_message()
            sys.exit(INTERRUPTION_ERROR_CODE)

    except (KeyboardInterrupt, SystemExit):
        p2.kill()
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)
    finally:
        os.unlink(in_sam.name)

    return outcome, TIME
# -----------------------------------------------------------------------------
def mapping(args, bowtie2_indexes, TIME, VERBOSE):
    '''
    Maps the input sample file (.fastq) into a .bam file (passing through a .sam file) using BowTie2 and Samtools commands

        Pipeline:
            1.  bowtie2 --very-sensitive --no-unal -x <SPECIE> -U <INPUT PATH> -p <NUMBER OF PROCESSORS>
            2.  samtools view -bS <INPUT SAM FILE>
            3.  samtools sort -m <AMOUNT OF MEMORY> - <OUTPUT BAM FILE>

            About BowTie2 command:
                --very-sensitive := -D 20 -R 3 -N 0 -L 2 -i S,1,0.5
                -D 20           Up to 20 consecutive seed extension attempts can "fail" before Bowtie 2 moves on, using the alignments found so far. A seed extension "fails" if it does not yield a new best or a new second-best alignment.
                -R 3            3 is the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds. When "re-seeding" Bowtie 2 simply chooses a new set of reads (same length, same number of mismatches allowed) at different offsets and searches for more alignments.
                -N 0            Sets the number of mismatches to allowed in a seed alignment during multiseed alignment.
                -L 2            Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive.
                -i S,1,0.5      Sets a function governing the interval between seed substrings to use during multiseed alignment. This function is f(x) := 1 + 0.5 * sqrt(x)
                --no-unal       Does not create BAM record for unaligned reads.
                -x <SPECIE>     The basename of the index for the reference genome.
                -U <INPUT>      Comma-separated list of files containing unpaired reads to be aligned.
                -S <OUTPUT>     File to write SAM alignments to ("-" == stdout).
            For major details look at http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#usage

            About Samtools commands:
                -bS             Input is in SAM format, output is in BAM format
                -m              Amount of memory it will be used
    '''
    # extracting args variable
    input_set = args['input']
    fastx = args['fastx']
    # if -f FASTQ_TAR_BZ2 or -f FASTQ_TAR_GZ or -f FASTQ_SRA, then the input file is an archive
    is_multi_file = args['input'][1] in ARCHIVE_FORMATS
    clade = bowtie2_indexes + args['clade']
    out_bam = args['out_bam']
    min_length = args['readLength']
    max_numof_mismatches = args['th_mismatches']
    bt2_options = args['bt2']
    memory = args['mGB']
    numof_proc = args['nproc']
    tmp_path = args['tmp']    
    
    outcome = (None, None)
    try:
        XN_FILTER = False if max_numof_mismatches <= -1 else True
        input_path     = input_set[0]
        input_format   = input_set[1]
        decompress_cmd = input_set[3]
        total = 0
        rejected = 0
        outcome = (None, None)
        if VERBOSE: print('[I] Opening ' + input_path)

        if is_multi_file:
            # decompress archive and concatenate all the files inside (both in only one command)
            decompress_cmd.append(input_path)
            if VERBOSE:
                print('[I] ' + ' '.join(decompress_cmd))
            p0 = subprocess.Popen(decompress_cmd, stdout=subprocess.PIPE)

        # bowtie2 --very-sensitive --no-unal -x <SPECIE> -U <INPUT PATH> -p <NUMBER OF PROCESSORS>
        # default: bt2_options = '--very-sensitive'

        bowtie2_cmd = [ 'bowtie2' ] + list(filter(None,bt2_options.split('/'))) + [ '--no-unal', '-x', clade, '-U', '-' if is_multi_file else input_path,
                        ] + ([] if int(numof_proc) < 2 else ['-p', str(numof_proc)])
        if not VERBOSE:
            bowtie2_cmd.append('--quiet')
        else:
            print('[I] ' + ' '.join(bowtie2_cmd))

        if fastx is 'fasta':
            bowtie2_cmd.append('-f') #bowtie2 default is fastq (-q)

        if is_multi_file:
            p1 = subprocess.Popen(bowtie2_cmd, stdin=p0.stdout, stdout=subprocess.PIPE)
        else:
            p1 = subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE)

        if tmp_path == None:
            tmp_sam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sam')
        else:
            tmp_sam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sam', dir=tmp_path)

        if VERBOSE:
            print('[I] Created temporary file ' + tmp_sam.name)
            print(WAIT_MESSAGE)
            print('[I] SAM records filtering: mismatches threshold is at ' + str(max_numof_mismatches) + ', length threshold is at ' + str(min_length))
            
        with tmp_sam:
            # @TODO or line.decode('utf-8')?
            for line in p1.stdout:
                total += 1
                l = line.decode('utf-8')
                if l.startswith('@'):
                    tmp_sam.write(line)
                # @TODO this elif could be useless
                elif line == '':
                    tmp_sam.write(line)
                    break
                else:
                    words = l.strip().split('\t')
                    read_length, numof_snp = len(words[9]), int(words[14].split(':')[-1])
                    # Too short
                    if read_length < min_length:
                        rejected += 1
                        #if VERBOSE:
                        #    print('Filter out read #' + str(total) + ': length is ' + str(readLength))
                    # Too many mismatches
                    elif XN_FILTER:
                        if numof_snp > max_numof_mismatches:
                            rejected += 1
                            #if VERBOSE:
                            #    print('Filter out read #' + str(total) + ': found ' + str(numof_snp) + ' mismatches')
                    # Accept the read
                    else:
                        tmp_sam.write(line)

        if VERBOSE:
            print('[I] Rejected ' + str(rejected) + ' reads over ' + str(total) + ' total')
            TIME = time_message(TIME, 'Bowtie2 mapping and SAM filtering completed.')

        outcome, TIME = samtools_sam2bam(tmp_sam, out_bam, memory, tmp_path, TIME, VERBOSE)

        p1.stdout.close()

    except (KeyboardInterrupt, SystemExit):
        p1.kill()
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)

    return outcome, TIME
# -----------------------------------------------------------------------------
def piling_up(bam_file, isTemp, csv_file, TIME, VERBOSE):
    '''
    Create the indexes and then call the Samtool's mpileup command

        CSV file: meaning of the columns

            Each line consists of 5 (or optionally 6) tab-separated columns:
                1   Sequence identifier
                2   Position in sequence (starting from 1)
                3   Reference nucleotide at that position
                4   Number of aligned reads covering that position (depth of abundance)
                5   Bases at that position from aligned reads
                6   quality of those bases (OPTIONAL)

            Column 5: The bases string
                . (dot)                 means a base that matched the reference on the forward strand
                , (comma)               means a base that matched the reference on the reverse strand
                AGTCN                   denotes a base that did not match the reference on the forward strand
                agtcn                   denotes a base that did not match the reference on the reverse strand
                +[0-9]+[ACGTNacgtn]+    denotes an insertion of one or more bases
                -[0-9]+[ACGTNacgtn]+    denotes a deletion of one or more bases
                ^ (caret)               marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality
                $ (dollar)              marks the end of a read segment
                * (asterisk)            is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation

            (See also at http://samtools.sourceforge.net/pileup.shtml)
    '''

    try:
        # 4th command: samtools index <INPUT BAM FILE>
        index_cmd = ['samtools', 'index', bam_file]
        if VERBOSE:
            print('[I] ' + ' '.join(index_cmd))
        p4 = subprocess.Popen(index_cmd)
        if VERBOSE:
            print('[I] BAM file ' + bam_file + ' has been indexed')

        try:
            with open(csv_file, mode='w') as ocsv:
                # 5th command: samtools mpileup <INPUT BAM FILE> > <OUTPUT CSV FILE>
                mpileup_cmd = ['samtools', 'mpileup', bam_file]
                if VERBOSE:
                    print('[I] ' + ' '.join(mpileup_cmd) + ' > ' + csv_file)
                try:
                    p5 = subprocess.Popen(mpileup_cmd, stdout=ocsv)
                    p5.wait()
                except Exception as err:
                    show_error_message(err)
                    sys.exit(SAMTOOLS_ERROR_CODE)

            # If the .bam file is created as temporary file (not defined by the user), then delete it with the index .bai file
            if isTemp:
                os.unlink(bam_file)
                os.unlink(bam_file + '.bai')
            # Otherwise remove only the index .bai file
            else:
                os.remove(bam_file + '.bai')

            if VERBOSE:
                TIME = time_message(TIME, 'Samtools piling up (view+mpileup) completed.')
            return TIME

        except (KeyboardInterrupt, SystemExit):
            p5.kill()
            if isTemp:
                os.unlink(bam_file)
            show_interruption_message()
            sys.exit(INTERRUPTION_ERROR_CODE)

    except (KeyboardInterrupt, SystemExit):
        p4.kill()
        if isTemp:
            os.unlink(bam_file)
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)
# -----------------------------------------------------------------------------
def remapping(input_pair, out_bam, max_numof_mismatches, memory, tmp_path, TIME, PLATFORM, VERBOSE):
    '''
    Convert a BAM file back to a SAM file, then re-filter it, then re-bam it
    '''
    total, rejected = 0, 0
    input_path = input_pair[0]
    if VERBOSE:
        print('[I] Input BAM file is ' + input_path + '. Now convert into SAM...')
    convertback_cmd = ['samtools', 'view', '-h', '-o', '-', input_path] # samtools view -h -o out.sam in.bam
    if VERBOSE:
        print('[I] ' + ' '.join(convertback_cmd))
    if tmp_path == None:
        tmp_sam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sam')
    else:
        tmp_sam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sam', dir=tmp_path)

    p1 = subprocess.Popen(convertback_cmd, stdout=subprocess.PIPE)

    if VERBOSE:
        print('[I] SAM records filtering: mismatches threshold is at ' + str(max_numof_mismatches))
    with tmp_sam:
        for line in p1.stdout:
            total += 1
            l = line.decode('utf-8')
            if l.startswith('@'):
                tmp_sam.write(line)
            # @TODO this elif could be useless
            elif line == '':
                tmp_sam.write(line)
                break
            else:
                words = l.strip().split('\t')
                numof_snp = int(words[14].split(':')[-1])
                # Accept the read
                if numof_snp <= max_numof_mismatches:
                    tmp_sam.write(line)
                # Too many mismatches
                else:
                    rejected += 1
                    #if VERBOSE:
                    #    print('Filter out read #' + str(total) + ': found ' + str(numof_snp) + ' mismatches')
    if VERBOSE:
        print('[I] Rejected ' + str(rejected) + ' reads over ' + str(total) + ' total')
        TIME = time_message(TIME, 'BAM->SAM reconversion and SAM filtering completed.')

    outcome, TIME = samtools_sam2bam(tmp_sam, out_bam, memory, tmp_path, TIME, VERBOSE)
    p1.stdout.close()

    return outcome, TIME
# ------------------------------------------------------------------------------
def genes_abundances(reads_file, contig2gene, gene_covs_outchannel, TIME, VERBOSE):
    '''
    Compute the abundance for each gene
        Workflow:
            1-  for each read R:
                    take the contig C the read belongs
                    take the position P
                    take its abundance V
            2-  take the list of genes in C
            3-  if P is inside G, then the abundance of G increases by V
    '''
    try:
        genes_abundances = defaultdict(int)
        current_contig = ''
        pos2gene = defaultdict(list)
        if VERBOSE: print(WAIT_MESSAGE)

        with open(reads_file, mode='r') as IN:
            for line in IN:
                # words = CONTIG, POSITION, REFERENCE BASE, COVERAGE, READ BASE, QUALITY
                words = line.strip().split('\t')
                contig, position, abundance = words[0], int(words[1]), int(words[3])
                # File is organised by contigs
                # If we have just passed from a contig to a new one, ...
                if contig != current_contig:
                    # Clear the {position:gene} dictionary (we allocate only a subdictionary for the interesting contig)
                    pos2gene = defaultdict(list)
                    if contig in contig2gene:
                        # For each gene in the contig, create {position:gene}
                        for gene, (fr,to) in contig2gene[contig].items():
                            for a in range(fr, to + 1):
                                pos2gene[a].append(gene)
                        current_contig = contig
                        if VERBOSE:
                            print('[I] Analyzing contig ' + current_contig + '...')

                # Add abundance for each gene covering the position
                if position in pos2gene:
                    genes = pos2gene[position]
                    for g in genes:
                        genes_abundances[g] += abundance

        # Print the gene abundances in stdout
        if gene_covs_outchannel == '-':
            for g in genes_abundances:
                if genes_abundances[g] > 0:
                    sys.stdout.write(str(g) + '\t' + str(genes_abundances[g]) + '\n')
        # Or in the output file (if user-defined)
        else:
            try:
                ocsv = open(gene_covs_outchannel, mode='w')
                for g in genes_abundances:
                    if genes_abundances[g] > 0:
                        ocsv.write(str(g) + '\t' + str(genes_abundances[g]) + '\n')
                ocsv.close()
                # with open(gene_covs_outchannel, 'rb') as icsv:
                try:
                    icsv = open(gene_covs_outchannel, 'rb')
                        # with bz2.BZ2File(gene_covs_outchannel + '.bz2', 'wb', compresslevel=9) as obz2:
                    obz2 = bz2.BZ2File(gene_covs_outchannel + '.bz2', 'wb', compresslevel=9)
                    copyfileobj(icsv, obz2)
                    obz2.close()
                    icsv.close()
                    os.remove(gene_covs_outchannel)

                except (KeyboardInterrupt, SystemExit):
                    os.remove(gene_covs_outchannel + '.bz2')
                    show_interruption_message()
                    sys.exit(INTERRUPTION_ERROR_CODE)

            except (KeyboardInterrupt, SystemExit):
                os.remove(gene_covs_outchannel)
                show_interruption_message()
                sys.exit(INTERRUPTION_ERROR_CODE)

    except (KeyboardInterrupt, SystemExit):
        os.unlink(reads_file)
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)

    if VERBOSE:
        TIME = time_message(TIME, 'Gene abundances computing has just been completed.')
    return genes_abundances

# ------------------------------------------------------------------------------
# CHECKS
# ------------------------------------------------------------------------------

def check_fastqdump(VERBOSE, PLATFORM):
    '''
    If input is a SRA file: check if SRA-toolkit (fastq-dump tool) is installed
    '''
    try:
        fastqdump = ''
        if PLATFORM == WINDOWS:
            fastqdump = subprocess.Popen(['where', 'fastq-dump'], stdout=subprocess.PIPE).communicate()[0]
        else: # Linux, Mac, ...
            fastqdump = subprocess.Popen(['which', 'fastq-dump'], stdout=subprocess.PIPE).communicate()[0]
        if VERBOSE:
            print('[I] SRA-toolkit is installed in the system in path ' + str(fastqdump))
    except Exception as err:
        show_error_message(err)
        print('[W] SRA-toolkit is not installed. SRA sample files cannot be processed.')
        sys.exit(UNINSTALLED_ERROR_CODE)
# -----------------------------------------------------------------------------
def check_samtools(VERBOSE = False, PLATFORM = 'lin'):
    '''
    Check if Samtools is installed
    '''
    try:
        samtools = ''
        if PLATFORM == WINDOWS:
            samtools = subprocess.Popen(['where', 'samtools'], stdout=subprocess.PIPE).communicate()[0].decode()
        else: # Linux, Mac, ...
            samtools = subprocess.Popen(['which', 'samtools'], stdout=subprocess.PIPE).communicate()[0].decode()
        # samtools_version = subprocess.Popen(['samtools', '--version'], stdout=subprocess.PIPE).communicate()[0]
        # samtools_version = samtools_version.decode().split(os.linesep)[0].split()[1]
        ### update: check also for older 0.1.19 version
        ### ('samtools --version' does not exist in 0.1.19, alternatively screen stderr of only 'samtools')
        # Version: 0.1.19-44428cd
        # Version: 1.1 (using htslib 1.1)
        # Version: 1.2 (using htslib 1.2.1)
        # Version: 1.3.1 (using htslib 1.3.1)
        samtools_stdout,samtools_stderr = subprocess.Popen(['samtools'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        samtool_lines = samtools_stderr.decode().split(os.linesep)
        samtools_version_line = [s for s in samtool_lines if 'Version' in s]
        samtools_version = samtools_version_line[0].split()[1].split('-')[0]
        if VERBOSE:
            print('[I] Samtools version ' + str(samtools_version) + ';  path: ' + str(samtools.strip()) )
    except Exception as err:
        # show_error_message(err)
        print('\n[E] Error: Cannot find Samtools, please install from http://www.htslib.org/ \n')
        sys.exit(UNINSTALLED_ERROR_CODE)
    return samtools, samtools_version
# -----------------------------------------------------------------------------
def check_bowtie2(clade, bowtie2_indexes, VERBOSE=False, PLATFORM='lin'):
    '''
    Check if Bowtie2 is alread installed
    '''
    try: # bowtie2 installed?
        if PLATFORM == WINDOWS:
            bowtie2 = subprocess.Popen(['where', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0].decode()
        else: # Linux, Mac, ...
            bowtie2 = subprocess.Popen(['which', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0].decode()
        bowtie2_version = subprocess.Popen(['bowtie2', '--version'], stdout=subprocess.PIPE).communicate()[0]
        bowtie2_version = bowtie2_version.decode().split()[2]
        if VERBOSE:
            print('[I] Bowtie2 version ' + str(bowtie2_version) + ';  path: ' + str(bowtie2).strip())
    except Exception as err:
        show_error_message(err)
        print('\n[E] Please install Bowtie2.\n')
        sys.exit(UNINSTALLED_ERROR_CODE)

    bt2_indexes = []
    if bowtie2: # check for bowtie2 index directory BOWTIE2_INDEXES
        try:
            # use the bowtie2_indexes option value, if set
            if bowtie2_indexes:
                bowtie2_indexes_dir = os.path.abspath(bowtie2_indexes)
            else:
                bowtie2_indexes_dir = '.'
            bt2_indexes = find(clade + '.[1-4].bt2', bowtie2_indexes_dir)
            bt2_indexes.extend( find(clade + '.rev.[1-2].bt2', bowtie2_indexes_dir) )
            if not len(bt2_indexes) == 6:
                # $BOWTIE2_INDEXES not defined in os.environment or indexes files (.bt2) for clade are not found
                bowtie2_indexes_dir = os.environ['BOWTIE2_INDEXES']
                bt2_indexes = find(clade + '.[1-4].bt2', bowtie2_indexes_dir)
                bt2_indexes.extend( find(clade + '.rev.[1-2].bt2', bowtie2_indexes_dir) )
                if not len(bt2_indexes) == 6:
                    raise IOError
        except KeyError:
            print('[E] Unable to find the bowtie2 indexes and pangenome! Use the option --i_bowtie2_indexes or set the environment variable BOWTIE2_INDEXES!')
            sys.exit(INDEXES_NOT_FOUND_ERROR)
        except IOError:
            print('[E] Bowtie2 index files (*.bt2) are not found! Use the option --i_bowtie2_indexes to set the location!')
            sys.exit(INDEXES_NOT_FOUND_ERROR)

        bowtie2_indexes_dir = os.path.join(os.path.dirname(bt2_indexes[0]),'') # '' to get ending '/'

        if VERBOSE:
            print('[I] BOWTIE2_INDEXES in ' + str(bowtie2_indexes_dir))

    return bowtie2, bowtie2_indexes_dir
# ------------------------------------------------------------------------------
def check_args():
    '''
    Check if the input arguments respect the rules of usage

        Usage examples:
            complete bowtie2 and abundance processing (single sample file)
             panphlan.py -c ecoli -i sample.fastq > sample_pangenome.csv
             panphlan.py -c ecoli -i sample.fastq --out_bam sample.bam > sample_pangenome.csv
             panphlan.py -c ecoli -i sample.fastq --out_bam sample.bam -o sample_pangenome.csv
             panphlan.py -c ecoli -i sample.fastq --out_bam sample.bam -o sample_pangenome.csv --input_format fastq

            same for already created .bam files
             panphlan.py -c ecoli -i sample.bam > sample_pangenome.csv

            same but using standard input and pipe
             tar -jxOf SAMPLE.tar.bz2 | panphlan.py -c ecoli > sample_pangenome.csv
    '''
    parser = PanPhlAnParser()
    args_set = vars(parser.parse_args())

    VERBOSE = args_set['verbose']

    if VERBOSE:
        print('\nPanPhlAn map version '+__version__)
        print('Python version: ' + sys.version.split()[0])
        print('System: ' + sys.platform)
        print(' '.join(sys.argv))

    # Check: INPUT_FILE -------------------------------------------------------
    is_compressed = False
    ipath = args_set['input']

    if ipath == None:
        print('[I] Input file is not specified. It will be used the standard input (stdin). Input format is "fastq" by default.')
        # Input format "fastq" is set by default because the stdin is a fasta/q stream
        # i.e. the program has not to decompress anything, the user does
        args_set['input'] = ('-', args_set['fastx'],'','')
        # FILE_EXTENSION = FASTQ
        FILE_EXTENSION = args_set['fastx']
    else:
        # If the file does not exist, halts the program
        if not os.path.exists(ipath):
            show_error_message('Specified input file does not exist.')
            sys.exit(INEXISTENCE_ERROR_CODE)
        else:
            # detect sample file format
            iextension, ifastx, idecompress = detect_input_format(ipath)
            FILE_EXTENSION = iextension

            if VERBOSE:
                print('[I] Input file: ' + ipath + '. Detected extension: ' + iextension + '; format: ' + ifastx)
            args_set['input'] = (ipath, iextension, ifastx, idecompress)
            if ifastx is 'fasta':         # only overwrite for clearly detected 'fasta',
                args_set['fastx']='fasta' # tar.bz2 can be both, user needs to specify if not default 'fastq'

    # Check: CLADE-NAME -------------------------------------------------------
    # replace "_" with "-" if present in species-name: "E-coli"
    args_set['clade'] = args_set['clade'].replace('panphlan_','') # remove panphlan_ prefix (to do: added later only for bowtie2)
    args_set['clade'] = args_set['clade'].replace('_','-') # convert underscore '_' to dash '-' in species-name (underscore is used as separator in _map)
    args_set['clade'] = PANPHLAN + args_set['clade'] # add again panphlan_ prefix as it is still expected in some functions of panphlan_map
    if VERBOSE: print('[I] Clade/Species: ' + args_set['clade'].replace(PANPHLAN,'') )

    # Check: OUTPUT_FILE ------------------------------------------------------
    opath = args_set['output']
    args_set['output'] = correct_output_name(opath, ipath, args_set['clade'], VERBOSE)
    if VERBOSE: print('[I] Output file name: ' + args_set['output'])

    # Check: OUTPUT_BAM_FILE --------------------------------------------------
    bpath = args_set['out_bam']
    # If --out_bam is defined, and input formt is different from BAM...
    if not bpath == None:
        if FILE_EXTENSION != BAM:

            # Create the path if not exists
            folders = os.path.dirname(bpath)
            if not os.path.exists(folders):
                os.makedirs(folders)
                if VERBOSE:
                    print('[I] Created path for (non-mandatory) output BAM file: ' + folders)

            # If out_bal defines a non-BAM file, then change automatically the file extension
            if not bpath.endswith(BAM):
                bpath = os.path.splitext(bpath)[0] + '.' + BAM
                args_set['out_bam'] = bpath
            if VERBOSE:
                print('[I] Output BAM file: ' + args_set['out_bam'])
    else:
        if VERBOSE:
            print('[I] BAM file will not be saved.')

    # Check: NUMOF_MISMATCHES -------------------------------------------------
    if VERBOSE:
        if args_set['th_mismatches'] >= 0:
            print('[I] Maximum number of mismatches: ' + str(args_set['th_mismatches']))
        # else:
        #    print('[I] Maximum number of mismatches: infinite.')

    # Check: NUMOF_PROCESSORS -------------------------------------------------
    # If we set less than 1 processors or more than the processors we really have, then set the minimum default
      # max_available_numof_processors = min(MAX_NUMOF_PROCESSORS, multiprocessing.cpu_count())
    max_available_numof_processors = multiprocessing.cpu_count()
    if args_set['nproc'] < 1 or args_set['nproc'] > max_available_numof_processors:
        args_set['nproc'] = max_available_numof_processors
        if VERBOSE: print('[W] Set number of processors to the maximal number on your machine: ' + str(args_set['nproc']))
    else:
        if VERBOSE: print('[I] Number of processors: ' + str(args_set['nproc']))

    # Check: MEMORY_GIGABTES_FOR_SAMTOOLS -------------------------------------
    # if args_set['mGB'] == None:
    #     args_set['mGB'] = 0.5 # new: default directly in options
    if VERBOSE:
        print('[I] GigaBytes for Samtools memory: -mGB ' + str(args_set['mGB']))

    # Check: READS_LENGTH -----------------------------------------------------
    if args_set['readLength'] < DEFAULT_READ_LENGTH:
        args_set['readLength'] = DEFAULT_READ_LENGTH
    if VERBOSE:
        print('[I] Minimum length threshold of reads: ' + str(args_set['readLength']))

    # Check: TEMP_FOLDER ------------------------------------------------------
    # default TMP_panphlan_map, since system folder /tmp can have space limits
    tmp_path = os.path.join(args_set['tmp'],'')
    if not os.path.exists(os.path.dirname(tmp_path)):
        os.makedirs(tmp_path)
    args_set['tmp'] = tmp_path
    if VERBOSE: print('[I] Temporary folder: ' + args_set['tmp'])
    # TMP folder is not removed after finishing, since it might be used by a parallel _map run.

    return args_set

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

def main():
    if sys.hexversion < 0x02060000:
        print('Python version: ' + sys.version)
        sys.exit('Python versions older than 2.6 are not supported.')
    args = check_args()
    TIME = time.time()
    TOTAL_TIME = time.time()
    VERBOSE = args['verbose']
    if not VERBOSE: print('\nUse option --verbose to display progress information.\n')
    MUST_REMAP = False if args['th_mismatches'] <= -1 else True
    PLATFORM = sys.platform.lower()[0:3]
    
    # --------------------------------------------------------------------------
    
    if VERBOSE: print('\nSTEP 1. Checking software...')
    bowtie2, bowtie2_indexes   = check_bowtie2(args['clade'], args['i_bowtie2_indexes'], VERBOSE, PLATFORM)
    samtools, samtools_version = check_samtools(VERBOSE, PLATFORM)
    pangenome_file = get_pangenome_file(bowtie2_indexes, args['clade'], VERBOSE)
    if FILE_EXTENSION == SRA:
        check_fastqdump(VERBOSE, PLATFORM)
        
    # --------------------------------------------------------------------------
    
    if VERBOSE: print('\nSTEP 2. Mapping short reads...')
    # If the input is a BAM file...
    if args['input'][1] == BAM:
        sorted_bam_file = args['input'][0]
        if MUST_REMAP:
            if VERBOSE: print('[I] BAM file as input argument. The BAM file needs to be remapped\.')
            # Convert BAM back to SAM, filter SAM basing on args['th_mismatches'], and reconvert into BAM and sort
            mapping_outcome, TIME = remapping(args['input'], args['out_bam'], args['th_mismatches'], args['mGB'], args['tmp'], TIME, PLATFORM, VERBOSE)
            isTemp = True if mapping_outcome[0] == TEMPORARY_FILE else False
            sorted_bam_file = mapping_outcome[1]
        else:
            if VERBOSE: print('[I] BAM file as input argument. Bowtie2 and Samtools will NOT be run to produce BAM file.')
    else:    
        mapping_outcome, TIME =  mapping(args, bowtie2_indexes, TIME, VERBOSE)
        sorted_bam_file = mapping_outcome[1]
        isTemp = True if mapping_outcome[0] == TEMPORARY_FILE else False
        
    # --------------------------------------------------------------------------
    
    if VERBOSE: print('\nSTEP 3. Piling up...')
    if args['tmp'] == None:
        tmp_csv__readsfile = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.csv')
    else:
        tmp_csv__readsfile = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.csv', dir=args['tmp'])
    TIME = piling_up(sorted_bam_file, isTemp, tmp_csv__readsfile.name, TIME, VERBOSE)

    # --------------------------------------------------------------------------
    
    if VERBOSE: print('\nSTEP 4. Exporting results...')
    try:
        contig2gene, TIME = build_pangenome_dicts(pangenome_file, TIME, VERBOSE)
        gene_abundances = genes_abundances(tmp_csv__readsfile.name, contig2gene, args['output'], TIME, VERBOSE)
        #os.unlink(tmp_csv__readsfile.name)
    except (KeyboardInterrupt, SystemExit):
        os.unlink(tmp_csv__readsfile.name)
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)
    finally:
        os.unlink(tmp_csv__readsfile.name)

    end_program(time.time() - TOTAL_TIME)
    
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
