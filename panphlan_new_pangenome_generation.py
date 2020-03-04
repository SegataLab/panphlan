#!/usr/bin/env python

from argparse import ArgumentParser
import os, subprocess, sys, tempfile, time
from fnmatch import fnmatch
from shutil import copyfileobj

__author__  = 'Leonard Dubois (panphlan-users@googlegroups.com)'
__version__ = '1.3'
__date__    = '11 October 2019'

# Error codes
UNINSTALLED_ERROR_CODE      = 2 # Software is not installed
INTERRUPTION_ERROR_CODE     = 7 # Computation has been manually halted
INTERRUPTION_MESSAGE    = '[E] Execution has been manually halted.\n'
# ------------------------------------------------------------------------------
# INTERNAL CLASSES
# ------------------------------------------------------------------------------
class PanPhlAnGenParser(ArgumentParser):
    '''
    Subclass of ArgumentParser for parsing command inputs for panphlan_pangenome_generation2.py
    '''
    def __init__(self):
        ArgumentParser.__init__(self)
        self.add_argument('-i', '--i_fna',      metavar='INPUT_FNA_FOLDER',     type=str,   default=False,      help='Folder containing the .fna genome sequence files')
        self.add_argument('-c','--clade',    metavar='CLADE_NAME',           type=str,   default=False,      help='Name of the species pangenome database, for example: -c ecoli17')
        self.add_argument('-o','--output',      metavar='OUTPUT_FOLDER',        type=str,   default='database', help='Result folder for all database files')
        self.add_argument('--verbose',          action='store_true',                                            help='Show progress information')
# ------------------------------------------------------------------------------
# MINOR FUNCTIONS
# ------------------------------------------------------------------------------
def show_interruption_message():
	sys.stderr.flush()
	sys.stderr.write('\r')
	sys.stderr.write(INTERRUPTION_MESSAGE)

def time_message(start_time, message):
    current_time = time.time()
    print('[I] ' + message + ' Execution time: ' + str(round(current_time - start_time, 2)) + ' seconds.')
    return current_time

# ------------------------------------------------------------------------------
# MAJOR FUNCTIONS
# ------------------------------------------------------------------------------
def check_bowtie2(VERBOSE):
    '''
    Check if Bowtie2 is installed
    '''
    PLATFORM = sys.platform.lower()[0:3]
    try:
        if PLATFORM == 'win':
            bowtie2 = subprocess.Popen(['where', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        else: # Linux, Mac, ...
            bowtie2 = subprocess.Popen(['which', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        bowtie2_version = subprocess.Popen(['bowtie2', '--version'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        bowtie2_version = bowtie2_version.split()[2]
        if VERBOSE:
            print('[I] Bowtie2 is installed, version: ' + str(bowtie2_version) + ', path: ' + str(bowtie2).strip())

    except OSError as err:
        sys.stderr.write('\n[E] Execution has encountered an error!\n')
        sys.stderr.write('    ' + str(err) + '\n')
        print('\n[E] Please, install Bowtie2.\n')
        if VERBOSE: print('    Bowtie2 is used to generate the .bt2 index files required in panphlan_map.py\n')
        sys.exit(UNINSTALLED_ERROR_CODE)
# ------------------------------------------------------------------------------
def create_bt2_indexes(args, TIME):

    try:
        genomes_files = [f for f in os.listdir(args['i_fna']) if fnmatch(f,'*.fna')]
        tmp_cat_filename = args['output'] + args['clade'] + '_ref_genomes.fna'
        out_all_genomes = open(tmp_cat_filename, "a+")
        for f in genomes_files:
            IN = open(os.path.join(args['i_fna'], f), "r")
            copyfileobj(IN, out_all_genomes)
            IN.close()
        out_all_genomes.close()

        try:
            ind_prefix = args['output'] + 'panphlan_' + args['clade']
            build_cmd = ['bowtie2-build', tmp_cat_filename, ind_prefix]
            if args['verbose']: print('[C] ' + ' '.join(build_cmd))
            p1 = subprocess.Popen(build_cmd)
            p1.wait()

            try:
                # Check generated files
                inspect_cmd = ['bowtie2-inspect', '-n', ind_prefix]
                if not args['verbose']: inspect_cmd.append('--verbose')
                else:
                    print('[I] ' + ' '.join(inspect_cmd))
                p2 = subprocess.Popen(inspect_cmd)
                p2.wait()
                if args['verbose']: TIME = time_message(TIME, 'Bowtie2 indexes have been created.')
                return TIME

            except (KeyboardInterrupt, SystemExit):
                p2.kill()
                show_interruption_message()
                sys.exit(INTERRUPTION_ERROR_CODE)
        except (KeyboardInterrupt, SystemExit):
            p1.kill()
            show_interruption_message()
            sys.exit(INTERRUPTION_ERROR_CODE)
    except (KeyboardInterrupt, SystemExit):
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)
# ------------------------------------------------------------------------------
def check_args():

    parser = PanPhlAnGenParser()
    args = vars(parser.parse_args())
    if not args['verbose']: print('\nUse option --verbose to display progress information.\n')

    if args['clade']:
        pangenome_filename = ''.join((args['i_fna'], 'panphlan_', args['clade'], "_pangenome.csv"))
        if not os.path.exists(pangenome_filename):
            sys.exit('\n Error: Pangenome file ('+ pangenome_filename +') not found\n')
    else:
        sys.exit('\n Error: Please provide clade name.\n')

    if args['i_fna']:
        in_path = os.path.join(args['i_fna'], '')
        if not os.path.exists(in_path):
            sys.exit('\n Error: Input folder not found.\n')
    else:
        sys.exit('\n Error: Please provide input folder.\n')

    if args['output']:
        args['output'] = os.path.join(args['output'], '')
        if not os.path.exists(args['output']):
            os.makedirs(args['output'])

    return args
# ------------------------------------------------------------------------------
def main():

    if not sys.version_info.major == 3:
        print('Python version: ' + sys.version)
        sys.exit('This software uses Python3, please update Python')

    args = check_args()
    start_time = time.time()
    TIME = time.time()
    check_bowtie2(args['verbose'])
    TIME = create_bt2_indexes(args, TIME)

    print('[TERMINATING...] ' + __file__ + ', ' + str(round((time.time() - start_time) / 60.0, 2)) + ' minutes.')
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
