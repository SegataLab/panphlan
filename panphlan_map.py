#!/usr/bin/env python

"""
panphlan_map.py
    Using bowtie2 indexes generated using panphlan_prepare_indexes.py script to map a metagenome sample to a pangenome.
"""

import os, subprocess, sys, time, bz2, tempfile
import argparse as ap
from collections import defaultdict
from shutil import copyfileobj

from misc import check_bowtie2

__author__ = 'Leonard Dubois, Matthias Scholz, Thomas Tolio and Nicola Segata (contact on https://forum.biobakery.org/)'
__version__ = '3.0'
__date__ = '20 April 2020'


DEFAULT_MIN_READ_LENGTH = 70

# ------------------------------------------------------------------------------
"""
Reads and parses the command line arguments of the script.
:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-i', '--input', type = str, default=sys.stdin,
                   help='Metagenomic sample to map')
    p.add_argument('--indexes', type = str, 
                   help='Bowtie2 indexes path and file prefix')         
    p.add_argument('-p', '--pangenome', type = str,
                   help='Path to pangenome tsv file exported from ChocoPhlAn')
    p.add_argument('-o', '--output', type = str, default=None,
                   help='Path to output file')
    p.add_argument('--bt2', type=str, default='--very-sensitive',
                   help='Additional bowtie2 mapping options, separated by slash: /-D/20/-R/3/, default: -bt2 /--very-sensitive/')
    p.add_argument('-b','--out_bam', type=str, default=None,
                   help='Get BAM output file')
    p.add_argument('--nproc', type=int, default=12,
                   help='Maximum number of processors to use. Default is 12 or a lower number of available processors.')
    p.add_argument('--min_read_length', type=int, default=DEFAULT_MIN_READ_LENGTH,
                   help='Minimum read length, default 70')
    p.add_argument('--th_mismatches', type=int, default=-1, 
                   help='Number of mismatches to filter (bam)')
    p.add_argument('-m', '--sam_memory', type=float, default=4.0,
                   help='Maximum amount of memory for Samtools (in Gb). Default 4')
    p.add_argument('--fasta', action='store_true',
                   help='Read are fasta format. By default considered as fastq')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='Show progress information')
    return p.parse_args()


"""Check arguments consistency"""
def check_args(args):
    
    if args.input:
        if not os.path.exists(args.input):
            sys.exit('[E] Sample file (' + args.input + ') not found\n')
    else:
        sys.exit('[E] Please provide a valid sample file (argument -i or --input).\n')

    if args.pangenome:
        if not os.path.exists(args.pangenome):
            sys.exit('[E] Pangenome file (' + args.pangenome + ') not found\n')
    else:
        sys.exit('[E] Please provide a valid pangenome file (argument -p or --pangenome).\n')
    
# ------------------------------------------------------------------------------
#   STEP 1 
# ------------------------------------------------------------------------------
"""Check if Samtools is installed. Stops programm if not"""
def check_samtools():
    platform = sys.platform.lower()[0:3]
    try:
        samtools = ''
        if platform == 'win':
            samtools = subprocess.Popen(['where', 'samtools'], stdout=subprocess.PIPE).communicate()[0].decode()
        else: # Linux, Mac, ...
            samtools = subprocess.Popen(['which', 'samtools'], stdout=subprocess.PIPE).communicate()[0].decode()
        samtools_stdout,samtools_stderr = subprocess.Popen(['samtools'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        samtool_lines = samtools_stderr.decode().split(os.linesep)
        samtools_version_line = [s for s in samtool_lines if 'Version' in s]
        samtools_version = samtools_version_line[0].split()[1].split('-')[0]
        print('[I] Samtools version ' + str(samtools_version) + ';  path: ' + str(samtools.strip()) )
    except Exception as err:
        print('\n[E] Please, install Samtools.\n')
        print('     Cannot find Samtools, please install from http://www.htslib.org/\n')
        sys.exit()
    return samtools_version

# ------------------------------------------------------------------------------
#   STEP 2 
# ------------------------------------------------------------------------------
"""Get sample file name and extension"""
def check_input(input_path):
    decompress_cmd = {}
    decompress_cmd['tar.bz2']  = ['tar', '-jxOf']
    decompress_cmd['tar.gz']   = ['tar', '-zxOf']
    decompress_cmd['gz']       = ['gunzip', '-c']
    decompress_cmd['bz2']      = ['bzcat']
    decompress_cmd['sra']      = ['fastq-dump', '-Z', '--split-spot', '--minReadLen', str(DEFAULT_MIN_READ_LENGTH)]
    for extension in decompress_cmd.keys():
        if input_path.endswith(extension):
            to_do = decompress_cmd[extension]
            to_do.append(input_path)
            print('[I] ' + 'input_path')
            return to_do
    return None
    

"""Convert a SAM file into BAM file, then sort the BAM"""
def samtools_sam2bam(in_sam, args):
    """samtools sort
          samtools version 1.2
            samtools sort <in.bam> <out.prefix>
            cat sample.bam | samtools sort - tmp_sorted
          samtools version 1.3
            samtools sort <in.bam> -o <out.bam>
            cat sample.bam | samtools sort - -o tmp_sorted.bam
        About Samtools commands:
            -bS             Input is in SAM format, output is in BAM format
            -m              Amount of memory it will be used
    """
    outcome = (None, None)
    try:
        samtools_version = check_samtools()
        # 1st command: samtools view -bS <INPUT SAM FILE>
        view_cmd = ['samtools', 'view', '-bS', in_sam.name]
        print('[I] ' + ' '.join(view_cmd))
        p2 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
        if args.verbose: print('[I] Temporary .bam file has been generated')

        try:
            # 2nd command: samtools sort -m <AMOUNT OF MEMORY> - <OUTPUT BAM FILE>
            sort_cmd = ['samtools', 'sort', '-m', str(int(args.sam_memory * 1024*1024*1024))]

            if args.out_bam == None: # .bam file is not saved, only temporary bam file
                tmp_bam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.bam')
                is_tmp = True
                out_bam = tmp_bam.name
            else:
                is_tmp = False
            
            sort_cmd += ['-', '-o', out_bam]
            if args.verbose: print('[I] cmd (v'  + samtools_version + '): ' + ' '.join(sort_cmd))
            with open(out_bam, mode='w') as OUT:
                p3 = subprocess.Popen(sort_cmd, stdin=p2.stdout, stdout=OUT)
                p3.wait() # Wait until previous process has finished its computation (otherwise there will be error raised by Samtools)
            if args.verbose:
                print('[I] Temporary .bam file ' + tmp_bam.name + ' has been sorted')
                print('[I] User-defined .bam file ' + out_bam + ' has been sorted')
                print('Samtools SAM->BAM translation (view+sort) completed.')        
            outcome = (is_tmp, out_bam)
                
        except (KeyboardInterrupt, SystemExit):
            p3.kill()
            sys.stderr.flush()
            sys.stderr.write('\r')
            sys.exit('[E] Execution has been manually halted.\n')
    except (KeyboardInterrupt, SystemExit):
        p2.kill()
        sys.stderr.flush()
        sys.stderr.write('\r')
        sys.exit('[E] Execution has been manually halted.\n')
    finally:
        os.unlink(in_sam.name)
    return outcome


"""Maps the input sample file (.fastq) into a .sam file using BowTie2 """
def mapping(args):
    """Pipeline:
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
    """
    
    bt2_options = args.bt2
    try:
        preprocess_cmd = check_input(args.input)
        if preprocess_cmd:
            p0 = subprocess.Popen(preprocess_cmd, stdout=subprocess.PIPE)
        # bowtie2 --very-sensitive --no-unal -x <SPECIE> -U <INPUT PATH> -p <NUMBER OF PROCESSORS>
        # default: bt2_options = '--very-sensitive'
        bowtie2_cmd = ([ 'bowtie2' ] + 
                    list(filter(None, bt2_options.split('/'))) +
                    [ '--no-unal', '-x', args.indexes, '-U', '-' if preprocess_cmd else args.input] +
                    ([] if int(args.nproc) < 2 else ['-p', str(args.nproc)]))
        if not args.verbose: bowtie2_cmd.append('--quiet')
        if args.fasta: bowtie2_cmd.append('-f') #bowtie2 default is fastq (-q)
        print('[I] ' + ' '.join(bowtie2_cmd))
        if preprocess_cmd:
            p1 = subprocess.Popen(bowtie2_cmd, stdin=p0.stdout, stdout=subprocess.PIPE)
        else:
            p1 = subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE)
        tmp_sam = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sam')
        if args.verbose:
            print('[I] Created temporary file ' + tmp_sam.name)
            print('[W] Please wait. The computation may take several minutes...')
            print('[I] SAM records filtering: mismatches threshold is at ' +
                  str(args.th_mismatches) + ', length threshold is at ' + str(args.min_read_length))    
        # Now, filter SAM
        total = 0
        rejected = 0
        with tmp_sam:
            for line in p1.stdout:
                total += 1
                l = line.decode('utf-8')
                if l.startswith('@'): tmp_sam.write(line)
                elif line == '':
                    tmp_sam.write(line)
                    break
                else:
                    words = l.strip().split('\t')
                    read_length, numof_snp = len(words[9]), int(words[14].split(':')[-1])
                    if read_length < args.min_read_length: # Too short
                        rejected += 1
                        if args.verbose: print('Filter out read #' + str(total) + ': length is ' + str(readLength))
                    elif args.th_mismatches > -1: # Too many mismatches
                        if numof_snp > args.th_mismatches:
                            rejected += 1
                            if args.verbose: print('Filter out read #' + str(total) + ': found ' + str(numof_snp) + ' mismatches')
                    else: # Accept the read
                        tmp_sam.write(line)
        print('[I] Rejected ' + str(rejected) + ' reads over ' + str(total) + ' total')
        print('Bowtie2 mapping and SAM filtering completed.')
        p1.stdout.close()
    except (KeyboardInterrupt, SystemExit):
        p1.kill()
        sys.stderr.flush()
        sys.stderr.write('\r')
        sys.exit('[E] Execution has been manually halted.\n')
    return tmp_sam

# ------------------------------------------------------------------------------
#   STEP 3
# ------------------------------------------------------------------------------
def piling_up(bam_file, is_tmp, csv_file, args):
    """Create the indexes and then call the Samtool's mpileup command
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
    """
    try:
        # command: samtools index <INPUT BAM FILE>
        index_cmd = ['samtools', 'index', bam_file]
        print('[I] ' + ' '.join(index_cmd))
        p4 = subprocess.Popen(index_cmd)
        if args.verbose: print('[I] BAM file ' + bam_file + ' has been indexed')
        try:
            with open(csv_file, mode='w') as ocsv:
                # command: samtools mpileup <INPUT BAM FILE> > <OUTPUT CSV FILE>
                mpileup_cmd = ['samtools', 'mpileup', bam_file]
                print('[I] ' + ' '.join(mpileup_cmd) + ' > ' + csv_file)
                try:
                    p5 = subprocess.Popen(mpileup_cmd, stdout=ocsv)
                    p5.wait()
                except Exception as err:
                    show_error_message(err)
                    sys.stderr.flush()
                    sys.stderr.write('\r')
                    sys.exit('[E] Samtools encountered some error.\n')   
            # delete tmp file
            if is_tmp: os.unlink(bam_file)
            os.unlink(bam_file + '.bai')
            if args.verbose: print('Samtools piling up (view+mpileup) completed.')
        except (KeyboardInterrupt, SystemExit):
            p5.kill()
            if isTemp: os.unlink(bam_file)
            sys.stderr.flush()
            sys.stderr.write('\r')
            sys.exit('[E] Execution has been manually halted.\n')
    except (KeyboardInterrupt, SystemExit):
        p4.kill()
        if isTemp: os.unlink(bam_file)
        sys.stderr.flush()
        sys.stderr.write('\r')
        sys.exit('[E] Execution has been manually halted.\n')
        
# ------------------------------------------------------------------------------
#   STEP 4
# ------------------------------------------------------------------------------
"""Build the dictionary for contig -> included gene -> location of the gene in the DNA"""
def build_pangenome_dicts(args):
    contig2gene = {}
    with open(args.pangenome, mode='r') as IN:
        for line in IN:
            # line = FAMILY, GENE_NAME, GENE_UID, CONTIG, FROM, TO
            words = line.strip().split('\t')
            fml, gen, ctg, fr, to = words[0], words[1], words[3], int(words[4]), int(words[5])
            if not ctg in contig2gene:
                contig2gene[ctg] = {}
            contig2gene[ctg][gen] = (min(to, fr), max(to, fr))
    if args.verbose: print('Dictionary for {contig:{gene:(from,to)}} has been created.')
    return contig2gene


"""Compute the abundance for each gene"""
def genes_abundances(reads_file, contig2gene, args):
    try:
        genes_abundances = defaultdict(int)
        if args.verbose: print('[W] Please wait. The computation may take several minutes...')
        # READ
        with open(reads_file, mode='r') as IN:
            for line in IN:
                words = line.strip().split('\t')
                # words = CONTIG, POSITION, REFERENCE BASE, COVERAGE, READ BASE, QUALITY
                contig, position, abundance = words[0], int(words[1]), int(words[3])
                # For each gene in the contig, if position in range of gene, increase its abundance
                if contig in contig2gene.keys():
                    for gene, (fr,to) in contig2gene[contig].items():
                        if position in range(fr, to+1):
                            genes_abundances[gene] += abundance
        # WRITE 
        if args.output == None:
            for g in genes_abundances:
                if genes_abundances[g] > 0:
                    sys.stdout.write(str(g) + '\t' + str(genes_abundances[g]) + '\n')
        else: 
            
            # WRITE AND THEN COMPRESS WITH copyobj()
            with bz2.open(args.output + '.bz2', 'wt', compresslevel=9) as OUT:
                for g in genes_abundances:
                    if genes_abundances[g] > 0:
                        OUT.write(str(g) + '\t' + str(genes_abundances[g]) + '\n')
    except (KeyboardInterrupt, SystemExit):
        os.unlink(reads_file)
        sys.stderr.flush()
        sys.stderr.write('\r')
        sys.exit('[E] Execution has been manually halted.\n')
    if args.verbose: print('Gene abundances computing has just been completed.')

# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------

def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python3, please update Python')       
    
    args = read_params()
    check_args(args)
    
    if args.verbose: print('\nSTEP 1. Checking software...')
    check_bowtie2()
    samtools_version = check_samtools()
    
    if args.verbose: print('\nSTEP 2.  Mapping the reads...')
    tmp_sam =  mapping(args)
    is_tmp, out_bam = samtools_sam2bam(tmp_sam, args)
    
    if args.verbose: print('\nSTEP 3. Piling up...')
    tmp_csv = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.csv')
    piling_up(out_bam, is_tmp, tmp_csv.name, args)
    
    if args.verbose: print('\nSTEP 4. Exporting results...')
    contig2gene = build_pangenome_dicts(args)
    genes_abundances(tmp_csv.name, contig2gene, args)
    os.unlink(tmp_csv.name)
    
    
if __name__ == '__main__':
    start_time = time.time() 
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
    