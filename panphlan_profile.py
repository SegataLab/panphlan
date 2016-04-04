#!/usr/bin/env python

from __future__ import with_statement 

# ==============================================================================
# PanPhlAn v1.2: PANgenome-based PHyLogenomic ANalysis
#                for detecting and characterizing strains in metagenomic samples
#
# Authors:  Matthias Scholz, algorithm design
#           Thomas Tolio, programmer
#           Nicola Segata, principal investigator
#
# PanPhlAn is a project of the Computational Metagenomics Lab at CIBIO,
# University of Trento, Italy
#
# For help type "./panphlan_map.py -h"
#
# https://bitbucket.org/CibioCM/panphlan
# ==============================================================================

__author__  = 'Matthias Scholz, Thomas Tolio, Nicola Segata (panphlan-users@googlegroups.com)'
__version__ = '1.2.0.1'
__date__    = '4 April 2016'

# Imports
from argparse import ArgumentParser
from collections import defaultdict
from random import randint
# from math import round 
import fnmatch, operator, os, subprocess, sys, time

try:
    import numpy
except ImportError as err:
    print('\n[E]',err) 
    print('\n[E] Please install the numpy module of Python\n')
    sys.exit(2)


# Pangenome CSV file constants
FAMILY_INDEX    = 0
GENE_INDEX      = 1
GENOME_INDEX    = 2
CONTIG_INDEX    = 3
FROM_INDEX      = 4
TO_INDEX        = 5

# Thresholds
MIN_NONPRESENT_TH   = 0.05
PRESENT_TH          = 0.5
MIN_PRESENT_TH      = 0.10
MIN_MULTICOPY_TH    = 0.15
LEFT_TH             = 1.25 # v1.0: 1.18 strain presence/absence filter (plateau curve)
RIGHT_TH            = 0.75 # v1.0: 0.82
COVERAGE_TH         = 2.0  # v1.0: 5.0
ZERO_NON_PLATEAU_TH = 0.20 # multistrain detection
RNA_MAX_ZERO_TH     = 10.0
SIMILARITY_TH       = 50.0

# Default tokens
DEFAULT_NP          = '-'
UNACCEPTABLE_NP     = ['NA', 'NaN', '1']
DEFAULT_NAN         = 'NA'
UNACCEPTABLE_NAN    = ['-', '1']

# File extensions
CSV     = 'csv'
TXT     = 'txt'
BZ2     = 'bz2'
GZ      = 'gz'
ZIP     = 'zip'
EXTENSIONS = [CSV, TXT, BZ2, GZ, ZIP]

# Error codes
INEXISTENCE_ERROR_CODE  =  1 # File or folder does not exist
PARAMETER_ERROR_CODE    =  2 # Dependent options are missed (e.g. If we define --sample_pairs, we MUST also define both --i_dna and --i_rna)
FILEFORMAT_ERROR_CODE   =  3 # file of DNA RNA sample pairs

# Strings
PANPHLAN        = 'panphlan_'
NO_RNA_FILE_KEY = '# NA #'
INTERRUPTION_MESSAGE    = '[E] Execution has been manually halted.\n'

# Plot's colors
# See also http://en.wikipedia.org/wiki/Html_color
PLOT_COLORS = [
        '#ff0000', '#800000', '#ffff00', '#808000', '#00ff00',
        '#008000', '#00ffff', '#008080', '#008080', '#0000ff',
        '#000080', '#ff00ff', '#800080', '#fa8072', '#ffa07a',
        '#dc143c', '#b22222', '#8b0000', '#ff69b4', '#ff1493',
        '#c71585', '#ff7f50', '#ff4500', '#ffa500', '#ffd700',
        '#bdb76b', '#9400d3', '#4b0082', '#483d8b', '#6a5acd',
        '#7fff00', '#32cd32', '#00fa9a', '#2e8b57', '#006400',
        '#20b2aa', '#4682b4', '#4169e1', '#ffdead', '#f4a460',
        '#d2691e', '#a52a2a', '#a0522d', '#b8860b', '#000000']
COLOR_GREY  = '#c0c0c0'

# ------------------------------------------------------------------------------
# INTERNAL CLASSES
# ------------------------------------------------------------------------------

class PanPhlAnJoinParser(ArgumentParser):
    '''
    Subclass of ArgumentParser for parsing command inputs for panphlan.py
    '''
    def __init__(self):
        ArgumentParser.__init__(self)
        self.add_argument('-i','--i_dna',               metavar='INPUT_DNA_FOLDER',             type=str,   default='',                     help='Input directory of panphlan_map.py results, containing SAMPLE.csv.bz2 files')
        self.add_argument('-c','--clade',               metavar='CLADE_NAME',                   type=str,   required=True,                  help='Panphlan species/clade database (e.g.: ecoli16)')
        self.add_argument('-o','--o_dna',               metavar='OUTPUT_FILE',                  type=str,                                   help='Write gene family presence/absence matrix: gene_presence_absence.csv')
        self.add_argument('--i_rna',                    metavar='INPUT_RNA_FOLDER',             type=str,                                   help='RNA-seq: input directory of RNA mapping results SAMPLE_RNA.csv.bz2')
        self.add_argument('--sample_pairs',             metavar='DNA_RNA_MAPPING',              type=str,                                   help='RNA-seq: list of DNA-RNA sequencing pairs from the same biological sample.')
        self.add_argument('--th_zero',                  metavar='MINIMUM_THRESHOLD',            type=float, default=None,                   help='Gene family presence/absence threshold: lower are non-present gene families.')
        self.add_argument('--th_present',               metavar='MEDIUM_THRESHOLD',             type=float, default=None,                   help='Gene family presence/absence threshold: higher are present gene families.')
        self.add_argument('--th_multicopy',             metavar='MAXIMUM_THRESHOLD',            type=float, default=None,                   help='Gene family presence/absence threshold: higher are multicopy gene families.')
        self.add_argument('--min_coverage',             metavar='MIN_COVERAGE_MEDIAN',          type=float, default=None,                   help='Minimum coverage threshold, default: 2X')
        self.add_argument('--left_max',                 metavar='LEFT_MAX',                     type=float, default=None,                   help='Strain presence/absence plateau curve threshold: left max [1.25]')
        self.add_argument('--right_min',                metavar='RIGHT_MIN',                    type=float, default=None,                   help='Strain presence/absence plateau curve threshold: right min [0.75]')
        self.add_argument('--rna_max_zeros',            metavar='RNA_MAX_ZEROES',               type=float, default=RNA_MAX_ZERO_TH,        help='Max accepted percent of zero coveraged gene-families (default: <10 %%).')          
        self.add_argument('--rna_norm_percentile',      metavar='RNA_NORM_PERCENTILE',          type=float, default=50,                     help='Percentile for normalizing RNA/DNA ratios')
        self.add_argument('--strain_similarity_perc',   metavar='SIMILARITY_PERCENTAGE',        type=float, default=SIMILARITY_TH,          help='Minimum threshold (percentage) for genome length to add a reference genome to presence/absence matrix (default: 50).')
        self.add_argument('--np',                       metavar='NON_PRESENCE_TOKEN',           type=str,   default='NP',                   help='User-defined string to mark non-present genes. [NP]')
        self.add_argument('--nan',                      metavar='NOT_A_NUMBER_TOKEN',           type=str,   default='NaN',                  help='User-defined string to mark multicopy and undefined genes. [NaN]')
        self.add_argument('--o_covplot',                metavar='COV_PLOT_NAME',                type=str,                                   help='Filename for gene-family coverage plot.')
        self.add_argument('--o_covplot_normed',         metavar='NOR_PLOT_NAME',                type=str,                                   help='Filename for normalized gene-family coverage plot.')
        self.add_argument('--o_cov',                    metavar='PANCOVERAGE_FILE',             type=str,                                   help='Write raw gene-family coverage matrix.')
        self.add_argument('--o_idx',                    metavar='DNA_INDEX_FILE',               type=str,                                   help='Write gene-family plateau definitions (1, -1, -2, -3)')
        self.add_argument('--o_rna',                    metavar='RNA_EXPRS_FILE',               type=str,                                   help='Write normalized gene-family transcription values (RNA-seq).')
        self.add_argument('--strain_hit_genes_perc',    metavar='GENEHIT_PERC_PER_STRAIN',      type=str,                                   help='Write overlap of gene-families between samples-strains and reference genomes.')
        self.add_argument('--add_strains',              action='store_true',                                                                help='Add reference genomes to gene-family presence/absence matrix.')
        self.add_argument('--interactive',              action='store_true',                                                                help='Plot coverage curves to screen, and not to a file.')
        self.add_argument('--verbose',                  action='store_true',                                                                help='Display progress information.')
        self.add_argument('-v', '--version',            action='version',   version="PanPhlAn version "+__version__+"\t("+__date__+")",     help='Prints the current PanPhlAn version and exits.')


# ------------------------------------------------------------------------------
# MINOR FUNCTIONS
# ------------------------------------------------------------------------------

def end_program(total_time):
    print('\n[TERMINATING...] ' + __file__ + ', ' + str(round(total_time / 60.0, 2)) + ' minutes.\n')

def show_interruption_message():
    sys.stderr.flush()
    sys.stderr.write('\r')
    sys.stderr.write(INTERRUPTION_MESSAGE)

def show_error_message(error):
    sys.stderr.write('[E] Execution has encountered an error!\n')
    sys.stderr.write('    ' + str(error) + '\n')

def time_message(start_time, message):
    current_time = time.time()
    print(' [I] ' + message + ' Execution time: ' + str(round(current_time - start_time, 2)) + ' seconds.')
    return current_time

def find(pattern, path):
    '''
    Find all the files in the path whose name matches with the specified pattern
    '''
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                target = os.path.join(root, name)
                result.append(target)
    return result

def check_output(opath, odefault, goal, VERBOSE):
    '''
    Check and, whenever necessary, create de novo the path (folders and/or file) for execution's outcome
    '''
    # If path does not exist, then create it
    if opath == None:
        # # no message needed about files, we don't save (confusing), only tell what panphlan saves
        # if VERBOSE:
        #    # print('[I] Output file is not specified for ' + goal + '. It will be used the standard value ' + odefault)
        #    print('[I] No output saved for ' + goal)
        return odefault

    else:
        # Create the path
        if not os.path.exists(os.path.dirname(opath)):
            try:
                folder = os.path.dirname(opath)
                if not folder == '':
                    os.makedirs(folder)
                    if VERBOSE:
                        print('[I] Created output directory: ' + folder)
            except FileNotFoundError as err:
                show_error_message(err)
                sys.exit(INEXISTENCE_ERROR_CODE)
        if VERBOSE:
            print('[I] Output file for ' + goal + ': ' + opath)
        return opath

def get_sampleID_from_path(sample_path, clade):
    # example: "path/to/mapping/result/ERR54632_ecoli14.csv.bz2" -> "ERR54632"
    # if no path, return same ID: "ERR54632" -> "ERR54632"
    filename = os.path.basename(sample_path)
    filename = filename.replace('.csv.bz2','')
    c=clade.replace('panphlan_','')
    sampleID = filename.replace( '_' + c ,'')
    return sampleID

def random_color(used):
    '''
    Get a random color for plotting coverage curves
    '''
    reset = False
    total = PLOT_COLORS
    available = [c for c in total if c not in used]
    # If we have no other available colors, than repeat the picking
    if len(available) == 0:
        available = total
        reset = True
    return (available[randint(0, len(available) - 1)], reset)

# ------------------------------------------------------------------------------
# MAJOR FUNCTIONS
# ------------------------------------------------------------------------------

def is_present(family, sample2family2presence):
    for s in sample2family2presence:
        if sample2family2presence[s][family]:
            return True
    return False
# -----------------------------------------------------------------------------
def print_presence_absence_profiles(sample_order, sample2family2presence, families, fileout, TIME, VERBOSE):
    '''
    Print presence/absence (1/0) family profiles (--o_dna option)
    a) used for printing database reference genomes only
    b) to do: used for printing sample profiles
    c) to do: used for printing sample + addStrains (ref. genomes) profiles
    sample_order - can be subset of samples or reference genomes
    families - all families present in at least one sample
    sample2family2presence - 2D dict dataset
    '''
    with open(fileout, mode='w') as csv:
        csv.write('\t' + '\t'.join(sample_order) + '\n')
        for f in families:
            csv.write(f)
            for s in sample_order:
                p = '\t1' if sample2family2presence[s][f] else '\t0'
                csv.write(p)
            csv.write('\n')
    TIME = time_message(TIME, 'Presence/absence profiles printed to ' + fileout)
# -----------------------------------------------------------------------------
def build_strain2family2presence(ref_genomes, families, genome2families, TIME, VERBOSE):
    '''
    Build the dictionary from strain to gene family to presence
    { STRAIN : { GENE FAMILY : PRESENCE(True or False) } }
    '''
    # Build the dictionary strain2family2presence
    strain2family2presence = defaultdict(dict)
    numof_strains = len(ref_genomes)
    i = 1
    for s in ref_genomes:
        if VERBOSE:
            print('[I] [' + str(i) + '/' + str(numof_strains) + '] Analysing reference strain ' + s + '...')
            i += 1
        for f in families:
            if f in genome2families[s]:
                strain2family2presence[s][f] = True
            else:
                strain2family2presence[s][f] = False
    if VERBOSE:
        TIME = time_message(TIME, 'Gene families presence/absence in strain reference genomes computed.')
    return TIME, strain2family2presence
# -----------------------------------------------------------------------------
def select_related_ref_genomes(ref_genomes, strain2family2presence, similarity, samples_panfamilies, families, TIME, VERBOSE):
    '''
    Select reference genomes similar to strains detected in samples
    to add them in the gene-family presence/absence matrix.
    '''
    # Rejection 2 (vertical filtering)
    numof_strains = len(ref_genomes)
    rejected_strains = []
    i = 1
    for s in ref_genomes:
        if VERBOSE:
            print('[I] [' + str(i) + '/' + str(numof_strains) + '] Analysing ref. genome ' + s + '...')
            i += 1
        f2p = strain2family2presence[s]
        # Get number of present gene families in strain
        strain_families = [f for f in f2p if f2p[f]]
        strain_length = len(strain_families)

        # # Too short genomes
        # lb, ub = genome_length + (genome_length / 10), genome_length - (genome_length / 10)
        # if lb <= strain_length and strain_length <= ub:
        #     rejected_strains.append(s)
        #     if VERBOSE:
        #         print('[W] Strain ' + s + ' is rejected because its genome is too short (size: ' + str(strain_length) + ')')
        # # This part of code is superfluous!

        # Check that at least half of the strain families is present in families
        half = int(strain_length * similarity / 100)
        numof_ss_families = 0
        for f in strain_families:
            if f in samples_panfamilies:
                numof_ss_families += 1
                if numof_ss_families >= half:
                    break # Exit from the loop to improve performances
        if numof_ss_families < half:
            rejected_strains.append(s)
            if VERBOSE:
                print('[W] Strain ' + s + ' is rejected because only ' + str(numof_ss_families) + ' families are present in the samples.')


    # Delete entries in the dictionary
    for s in rejected_strains:
        del(strain2family2presence[s])

    if VERBOSE:
        print('[I] ' + str(len(rejected_strains)) + ' strain genomes filtered out. Strains are: ' + ', '.join(rejected_strains))
    selected_strains = sorted([s for s in strain2family2presence if s not in rejected_strains])
    if VERBOSE:
        print('[I] Selected strains are: ' + ', '.join(selected_strains))

    # Rejections 3 (horizontal filtering): always-zero gene families
    never_present_families = []
    for f in families:
        always_zero = True
        for g in strain2family2presence:
            if strain2family2presence[g][f]: # == True == 1
                always_zero = False
                break
        if always_zero:
            never_present_families.append(f)
    if VERBOSE:
        TIME = time_message(TIME, str(len(never_present_families)) + ' never present gene families filtered out.')

    return TIME, selected_strains, never_present_families
# -----------------------------------------------------------------------------
def get_samples_panfamilies(families, sample2family2presence, TIME, VERBOSE):
    '''
    Get the sorted list of all the families present in the samples
    Can be a subset of the pangenome's set of families
    '''
    panfamilies = set()
    for f in families:
        for s in sample2family2presence:
            if sample2family2presence[s][f]:
                panfamilies.add(f)
                break
    if VERBOSE:
        TIME = time_message(TIME, 'Extracted ' + str(len(panfamilies)) + ' gene families present in the samples.')
    return TIME, sorted(panfamilies)
# -----------------------------------------------------------------------------
def rna_seq(out_channel, sample2family2dnaidx, dna_sample2family2cov, dna_accepted_samples, rna_samples_list, rna_sample2family2cov, rna_max_zeroes, dna2rna, dna_file2id, rna_id2file, families, np_symbol, nan_symbol, clade, rna_norm_percentile, TIME, VERBOSE):
    '''
    DESCRIPTION
        1.  convert DNA samples to get (1,-1,-2,-3) DNA index matrix and DNA coverage values
        2.  convert RNA samples to coverage values only
        3.  select DNA/RNA sample pair
                a.  reject all RNA values not in corresponding DNA plateau area "1" 
                b.  Normalization ...
                c.  set gene-families "-1" and "-2" as missing (NaN or NA?)
                d.  set non-present genes to "-" by default or by symbol specified at command line -np 0 or -np NP (the "np_symbol" given in input)
        4.  merge all RNA samples in a single matrix


    find DNA RNA sample pairs
        select gene-families present in both DNA and RNA datasets
        divide all coverage values: RNA/DNA  (for each samples and each gene-family)
        if DNA is zero, set RNA/DNA also zero (not as undefined)

    pre-result part1: data_sepidermidis_RNAseq_part1_RNAdivDNA.csv
        still RNA/DNA values for all gene-families (including non-plateau gene-families)
        zero lines are removed
        non-plateau samples removed (only selected DNA samples, in which species is present)
        bad RNA samples (low coverage) still not removed


    '''
    # Data from Step 1 are given in input (sample2family2dnaidx, dna_samples_covs_path)
    # Data from Step 2 are given in input (rna_samples_covs_path)    

    sample2family2rna_div_dna = defaultdict(dict)
    rna_samples = []
    rna_ids = []
    dna_accepted_samples = sorted([s for s in dna_accepted_samples if dna_accepted_samples[s]])
    
    for dna_sample in dna_accepted_samples:
        rna_sample = rna_id2file[dna2rna[dna_file2id[dna_sample]]]
        if not rna_sample == NO_RNA_FILE_KEY:
            rna_samples.append(dna_sample)
            rna_ids.append(dna_file2id[dna_sample])
            # For each family present in at least one sample, divide RNA coverage for the correlative DNA coverage
            for f in families:
                dna_cov = dna_sample2family2cov[dna_sample][f]
                if dna_cov == 0.0: # We avoid a division by zero :)
                    sample2family2rna_div_dna[dna_sample][f] = 0.0
                else:
                    rna_cov = rna_sample2family2cov[rna_sample][f]
                    sample2family2rna_div_dna[dna_sample][f] = rna_cov / dna_cov
                    
    # Reverse dictionaries
    rna2dna = dict((v,k) for (k,v) in dna2rna.items())
    rna_file2id = dict((v,k) for (k,v) in rna_id2file.items())
    dna_id2file = dict((v,k) for (k,v) in dna_file2id.items())

    # Define
    sample2family2median_norm = defaultdict(dict)
    rna_samples.sort()
    sample2zeroes = defaultdict(tuple)
    median = defaultdict(float)


    # Step 3.4) Normalization + filtering
    for dna_sample in rna_samples:
        sample2zeroes[dna_sample] = (0,0)

        # Take all the gene families belonging to the plateau and calculte the median of their RNA/DNA values
        plateau_rna_div_dna = [sample2family2rna_div_dna[dna_sample][f] for f in sample2family2rna_div_dna[dna_sample] if sample2family2dnaidx[dna_sample][f] == 1]
        # median[dna_sample] = numpy.median(plateau_rna_div_dna)
        median[dna_sample] = numpy.percentile(plateau_rna_div_dna, rna_norm_percentile) # default: 50

        if VERBOSE:
            print('[I] Median of plateau gene families RNA/DNA values: ' + str(median[dna_sample]))
        
        for f in families:
            # If the family is in the plateau, calculate median normalized RNA/DNA value
            if sample2family2dnaidx[dna_sample][f] == 1:
                sample2family2median_norm[dna_sample][f] = sample2family2rna_div_dna[dna_sample][f] / median[dna_sample]
                # Update the number of zeroes over the total families (belonging to the plateau)
                numof_zeroes, numof_families = sample2zeroes[dna_sample]
                numof_families += 1
                if sample2family2median_norm[dna_sample][f] == 0.0:
                    numof_zeroes += 1
                sample2zeroes[dna_sample] = (numof_zeroes, numof_families)
            # If not in the plateau, set to NaN
            elif sample2family2dnaidx[dna_sample][f] == -3:
                sample2family2median_norm[dna_sample][f] = np_symbol
            else:
                sample2family2median_norm[dna_sample][f] = nan_symbol
        sample2zeroes[dna_sample] = float(sample2zeroes[dna_sample][0]) / sample2zeroes[dna_sample][1]

    # Reject samples with too many zeros
    rnaseq_accepted_samples = []
    for s in sample2zeroes:
        perc = sample2zeroes[s] * 100.0
        if VERBOSE:
            print('[I] Percentage of zero values for sample ' + s + ': ' + str(perc) + '%')
        if perc <= rna_max_zeroes:
            rnaseq_accepted_samples.append(s)
            print('    Sample is accepted.')
        else:
            print('    Sample is rejected.')

    # Log nomalization
    sample2family2log_norm = defaultdict(dict)
    for dna_sample in rnaseq_accepted_samples:
        for f in families:
            v = sample2family2median_norm[dna_sample][f]
            if type(v) is str:
                sample2family2log_norm[dna_sample][f] = sample2family2median_norm[dna_sample][f]
            else:
                sample2family2log_norm[dna_sample][f] = 0.0 if v == 0.0 else (numpy.log2(v) / 10) + 1.0

    # Print
    rnaseq_accepted_samples.sort()
    rnaseq_accepted_ids = [get_sampleID_from_path(s, clade) for s in rnaseq_accepted_samples]
    if not out_channel == '':
        with open(out_channel, mode='w') as csv:
            csv.write('\t' + '\t'.join(rnaseq_accepted_ids) + '\n')
            for f in families:
                # Skip the never present gene families
                all_null = True
                for s in rnaseq_accepted_samples:
                    if not sample2family2log_norm[s][f] == np_symbol:
                        if not sample2family2log_norm[s][f] == nan_symbol:
                            all_null = False
                            break
                if not all_null:
                    csv.write(f)
                    for s in rnaseq_accepted_samples:
                        v = sample2family2log_norm[s][f]
                        if type(v) is float or type(v) is numpy.float64:
                            csv.write('\t' + str(format(v, '.3f')))
                        else:
                            csv.write('\t' + v)
                    csv.write('\n')

    if VERBOSE:
        TIME = time_message(TIME, 'RNA indexing executed.')
    return TIME
# ------------------------------------------------------------------------------
def strains_gene_hit_percentage(ss_presence, genome2families, accepted_samples, out_channel, clade, TIME, VERBOSE):
    '''
    TODO
    '''
    # NB. File's lines must match this pattern: STRAIN z list_of(x's)
    #     where x = percentage, z = total number of gene families in the strain
    # ss_presence = { STRAIN or SAMPLE : { GENE FAMILY : PRESENCE } }
    strain2sample2hit = {}
    ref_genomes = sorted(genome2families.keys())
    samples_list = sorted([get_sampleID_from_path(s, clade) for s in accepted_samples.keys() if accepted_samples[s]])

    if len(samples_list) > 0:
        # Populate { STRAIN : { SAMPLE : HIT PERCENTAGE } } 
        for strain in genome2families:
            strain2sample2hit[strain] = {}
            for sample in samples_list:
                strain2sample2hit[strain][sample] = 0
                for family in genome2families[strain]:
                    # Add 1 if the gene family of the strain is also in the sample
                    if ss_presence[sample][family]:
                        strain2sample2hit[strain][sample] += 1
                # Divide the number of hit genes by the total number of gene families in the strain
                numof_hits = strain2sample2hit[strain][sample]
                strain_len = len(genome2families[strain])
                # strain2sample2hit[strain][sample] = (z,x)
                strain2sample2hit[strain][sample] = float(numof_hits) / strain_len * 100.0

        # Write into a file
        try:
            with open(out_channel, mode='w') as ocsv:
                ocsv.write('strainID\tnumber_of_genes\t' + '\t'.join(samples_list) + '\n')
                for strain in ref_genomes:
                    numof_families = len(genome2families[strain])
                    ocsv.write(strain + '\t' + str(numof_families))
                    for sample in samples_list:
                        perc = strain2sample2hit[strain][sample]
                        ocsv.write('\t' + str(format(perc, '.1f')))
                    ocsv.write('\n')

        except (KeyboardInterrupt, SystemExit):
            os.remove(out_channel)

        if VERBOSE:
            TIME = time_message(TIME, 'Strains hit gene families percentages computed.')
    
    else:
        print('[W] No file has been written for strains gene hit percentages because there is no accepted samples.')
    return strain2sample2hit, TIME
# ------------------------------------------------------------------------------
def samples_strains_presences(sample2family2presence, ref_genomes, strain2family2presence, genome_length, similarity, out_channel, families, clade, TO_BE_PRINTED, TIME, VERBOSE):
    '''
    Compute gene families presence/absence for reference genomes and merge them with detected sample strain profiles

    Combine sample presence/absence matrix with strain presence/absence matrix
        1. merge: first samples columns, then strains columns (still keep all gene-families present in any strain)
        2. reject all strains which have less than 50% of its gene-families in common with the sample matrix. As total number of gene-families, we can use genome_length=2616 (saureus) for all strains. Means 50% = 1308 (saureus) gene-families of a strain have to be present in the sample set, otherwise strain is excluded.
        3. reject completely (all zero) non-present gene-families (present only in rejected strains) 
    NB.
        Some gene-families can be present in samples, but not in the selected (>50%) strains.
        Some gene-families can be present in selected strains, but not in samples (if a strain is selected, we show all of it's gene-families).
    '''
    sample2family2presence = dict((get_sampleID_from_path(k, clade), v) for (k,v) in sample2family2presence.items())

    # Get all present (in at least one sample) families
    TIME, samples_panfamilies = get_samples_panfamilies(families, sample2family2presence, TIME, VERBOSE)
    TIME, selected_strains, never_present_families = select_related_ref_genomes(ref_genomes, strain2family2presence, similarity, samples_panfamilies, families, TIME, VERBOSE)

    # Create sorted list of sample/strains names
    sample_and_strain_sorted_list = []
    sample_and_strain_sorted_list = sorted(sample2family2presence.keys())
    sample_and_strain_sorted_list.extend(selected_strains)

    # Merge the two dictionaries
    # sample_and_strain_presences = dict(sample2family2presence, **strain2family2presence)
    sample_and_strain_presences = {}
    for s in sample2family2presence:
        sample_and_strain_presences[s] = sample2family2presence[s]
    for s in strain2family2presence:
        sample_and_strain_presences[s] = strain2family2presence[s]

    if TO_BE_PRINTED and len(sample_and_strain_sorted_list) > 0:
        if VERBOSE: print(' [I] Print gene-family presence/absence matrix of samples and ref. genomes to: ' + out_channel)
        with open(out_channel, mode='w') as csv:
            # [get_sampleID_from_path(s, clade) for s in sample_and_strain_sorted_list]
            csv.write('\t' + '\t'.join(sample_and_strain_sorted_list) + '\n')
            for f in families:
                if f not in never_present_families:
                    csv.write(f)
                    for s in sample_and_strain_sorted_list:
                        p = '\t1' if sample_and_strain_presences[s][f] else '\t0'
                        csv.write(p)
                    csv.write('\n')

    if VERBOSE:
        TIME = time_message(TIME, 'Samples/strains gene families presence/absence matrix computed.')
    return sample_and_strain_presences, TIME
# ------------------------------------------------------------------------------
def presence_of(dna_index):
    return dna_index >= -1

def presence_to_str(presence):
    return '1' if presence else '0'

def get_genefamily_presence_absence(accepted_samples, dna_files_list, dna_file2id, sample2family2dnaidx, out_channel, families, clade, avg_genome_length, sample_stats, TIME, VERBOSE):
    '''
    Build the gene families presence/absence matrix.
        Take the DNA indexing matrix:
        gene family in sample has DNA index 1 or -1 ==> present (1)
        gene family in sample has DNA index -2 or -3 ==> NOT present (0)
    '''
    sample2family2presence = defaultdict(dict)
    
    dna_sample_ids = [dna_file2id[s] for s in dna_files_list if accepted_samples[s]]
    dna_files_list = [s for s in dna_files_list if accepted_samples[s]]

    if not out_channel == '' and len(dna_files_list) > 0:
        csv = open(out_channel, mode='w')
    
    if not out_channel == '' and len(dna_sample_ids) > 0:
        csv.write('\t' + '\t'.join([get_sampleID_from_path(s, clade) for s in dna_sample_ids]) + '\n')
    for f in families:
        #if sum(sample2family2presence[s][f] for s in sample2family2presence) > 0:
        line = f
        total_presence = False
        for s in dna_files_list:
            presence = presence_of(sample2family2dnaidx[s][f])
            sample2family2presence[s][f] = presence
            total_presence = total_presence or presence
            if accepted_samples[s]:
                if f in sample2family2dnaidx[s]:
                    line = line + '\t' + presence_to_str(presence)
                else:
                    line = line + '\t0'
        if not out_channel == '' and len(dna_files_list) > 0:
            if total_presence:
                csv.write(line + '\n')

    # get number of gene-families per sample (add to dict sample_stats)
    for s in sample2family2presence.keys():
        sampleID = get_sampleID_from_path(s, clade)
        numGeneFamilies = sum( sample2family2presence[s][f] for f in sample2family2presence[s] )
        sample_stats[sampleID].update({'numberGeneFamilies' : numGeneFamilies})
    if VERBOSE:
        print(' [I] Number of gene families per sample-specific strain:')
        for sampleID in sorted(sample_stats.keys()):
            if 'numberGeneFamilies' in sample_stats[sampleID]:  
                print('      ' + sampleID + '\t' + str(sample_stats[sampleID]['numberGeneFamilies']))
        print('      Average number of gene-families in reference genomes: ' + str(avg_genome_length))        
        
    if len(dna_files_list) > 0:
        if VERBOSE:
            print(' [I] Gene family presence/absence matrix is printed to ' + out_channel)
            TIME = time_message(TIME, 'Presence/absence matrix finished.')
    else:
        print('[W] No file has been written for gene families presence/absence because there is no accpeted samples.')
        print('    You can try very sensitive options:  --min_coverage 1 --left_max 1.70 --right_min 0.30')
        print('    Read more: https://bitbucket.org/CibioCM/panphlan/wiki/panphlan_profile_strain_detection')
    return sample2family2presence, sample_stats, TIME
# -----------------------------------------------------------------------------
def print_multistrain_warning(sample_stats, avg_genome_length, VERBOSE):
    print(' ')
    for sampleID in sorted(sample_stats.keys()): # check result of out-plateau zero threshold
        if sample_stats[sampleID]['Multistrain']:
            print('QUALITY WARNING: sample ' + sampleID + ' may contain multiple strains, PanPhlAn extracts the dominant strain')
    for sampleID in sorted(sample_stats.keys()): # check number of gene-families
        if sample_stats[sampleID]['strainIsPresent']:
            n = sample_stats[sampleID]['numberGeneFamilies']
            if n > 1.1 * avg_genome_length:
                print('QUALITY WARNING: gene-families of sample ' + sampleID + ' may come from multiple strains \n  number of gene-families: '+ str(n) +' is 10% higher than expected number (average of ref. genomes): ' + str(avg_genome_length))
            if n < 0.75 * avg_genome_length:
                print('QUALITY WARNING: sample ' + sampleID + ' shows too low number of gene-families, due to low coverage or multiple strains  \n  number of gene-families: '+ str(n) +' is 25% lower than expected number (average of ref. genomes): ' + str(avg_genome_length))
    print('\n')    
# -----------------------------------------------------------------------------
def index_of(min_thresh, med_thresh, max_thresh, normalized_coverage):
    '''
    Return the DNA index for the given median-normalized coverage
    '''
    if normalized_coverage < min_thresh:
        return -3
    elif normalized_coverage <= med_thresh:
        return -2
    elif normalized_coverage <= max_thresh:
        return  1
    else:
        return -1
# -----------------------------------------------------------------------------
def get_idx123_plateau_definitions(accepted_samples, sample2family2normcov, min_thresh, med_thresh, max_thresh, index_file, families, clade, TIME, VERBOSE=False):
    '''
    -o_idx HMP_saureus_DNAindex.csv

    To use later also in RNA-seq, we need an DNA index matrix containing 4 levels (1, -1, -2, -3)

    Take samples that passed plateau criteria and define index based on coverage level of gene-families
         1 means plateau area of gene-families
        -1 means multicopy core genes (left from plateau), present also in other species
        -2 means undefined gene-families between plateau-level and zero
        -3 means "clearly" non-present gene-families

    Settings
        th1=0.30  (lower are non-present genes)
        th2=0.70  (higher are plateau or multi-copy genes)
        th3=1.30  (higher are only multi-copy genes)

    Get DNA index
        X = median normalized coverage values
        DNAindex set to "-3" if (X < th1)
        DNAindex set to "-2" if (X >=th1) & (X <= th2)
        DNAindex set to  "1" if (X > th2) & (X <= th3)
        DNAindex set to "-1" if (X > th3)
    '''
    sample2family2dnaidx = defaultdict(dict)

    accepted_ids = sorted([get_sampleID_from_path(s, clade) for s in accepted_samples])
    id2sample = dict((get_sampleID_from_path(s, clade),s) for s in accepted_samples)

    for sample in accepted_samples:
        if VERBOSE:
            print(' [I] Get DNA -1,1,2,3 levels for sample ' + sample)
        sample_id = get_sampleID_from_path(sample, clade)
        for family in families:
            sample2family2dnaidx[sample][family] = index_of(min_thresh, med_thresh, max_thresh, sample2family2normcov[sample][family])

    if not index_file == '' and len(accepted_ids) > 0:
        with open(index_file, mode='w') as csv:
            csv.write('\t' + '\t'.join(accepted_ids) + '\n')
            for family in families:
                csv.write(family)
                for sample_id in accepted_ids:
                    sample = id2sample[sample_id]
                    csv.write('\t' + str(sample2family2dnaidx[sample][family]))
                csv.write('\n')

    elif len(accepted_ids) == 0:
        print('[W] No file has been written for DNA indexing because there is no accpeted samples.')

    if VERBOSE:
        TIME = time_message(TIME, 'DNA indexing executed.')

    return sample2family2dnaidx, TIME 
# -----------------------------------------------------------------------------
def strain_presence_plateau_filter(samples_coverages, num_ref_genomes, avg_genome_length, th_min_coverage, th_plateau_left_max, th_plateau_right_min, families, clade, TIME, VERBOSE=False):
    '''
    Plateau quality filter based on genes coverage curves

    - take each sample individually
    - sort gene-families by abundance (coverage values)
    - Filter 1) curve needs to have a median coverage higher than 2X ("min_coverage")
    - Filter 2) left plateau side needs to be lower than threshold "left_max",
                right plateau side needs to be higher than "right_min"
    
    filter settings
    genome-length = 2300  (saureus)  (different for each species)

    position_median        = 0.5  ( x genome length) of sorted gene-family vector
    position_plateau_left  = 0.30 (  x genome length)
    position_plateau_right = 0.70 (  x genome length)

    threshold_min_coverage      = 2     (at position_median)
    threshold_plateau_left_max  = 1.18  (at position_plateau_left)
    threshold_plateau_right_min = 0.82  (at position_plateau_right)    
    '''
    sample2accepted = {}
    sample2famcovlist = {} # { SAMPLE NAME : ( [ COVERAGE ], [ GENE FAMILY ] ) }
    median = {}
    sample2color = {}
    median_normalized_covs = defaultdict(list)
    norm_samples_coverages = defaultdict(dict)
    sample_stats = defaultdict(dict)

    # reduce expected number of gene-families, in case of only 1,2 or 3 ref. genomes in DB
    orig_avg_genome_length=avg_genome_length
    if   num_ref_genomes == 3:
        avg_genome_length=int(round(0.90 * avg_genome_length))
        print('[I] Decrease expected gene-families per sample strain to: ' + str(avg_genome_length) + ' (0.9*' + str(orig_avg_genome_length) +') due to only 3 ref. genomes in DB')
    elif num_ref_genomes == 2:
        avg_genome_length=int(round(0.85 * avg_genome_length))
        print('[I] Decrease expected gene-families per sample strain to: ' + str(avg_genome_length) + ' (0.85*' + str(orig_avg_genome_length) +') due to only 2 ref. genomes in DB)')
    elif num_ref_genomes == 1:
        avg_genome_length=int(round(0.75 * avg_genome_length))
        print('[I] Decrease expected gene-families per sample strain to: ' + str(avg_genome_length) + ' (0.75*' + str(orig_avg_genome_length) +') due to only 1 ref. genomes in DB)')

    # set default filter th's
    if th_plateau_left_max is None:
        th_plateau_left_max=LEFT_TH
    if th_plateau_right_min is None:
        th_plateau_right_min=RIGHT_TH    
    if th_min_coverage is None:
        th_min_coverage=COVERAGE_TH
    th_max_zero=ZERO_NON_PLATEAU_TH    
        
    if VERBOSE:
        print(' [I] Minimum median coverage threshold: '                          + str(th_min_coverage))
        print(' [I] Left maximum plateau threshold: '                             + str(th_plateau_left_max))
        print(' [I] Right minimum plateau threshold: '                            + str(th_plateau_right_min))
        print(' [I] Maximum zero non-plateau threshold (multistrain detection): ' + str(th_max_zero))
        
    # Take one sample a time
    for sample in sorted(samples_coverages.keys()):

        sample_id = get_sampleID_from_path(sample, clade)
        sampleID  = get_sampleID_from_path(sample, clade) # try to change keys
        d = samples_coverages[sample]

        # Take families coverage from sample and sort descendently by value (coverage)
        families_covs = sorted(d.items(), key=lambda x: x[1])
        families_covs = families_covs[::-1] # reverse
        del(d)
        # Compute the median
        median[sampleID] = numpy.median([p[1] for p in families_covs][:avg_genome_length])
        sample2famcovlist[sample] = ([p[1] for p in families_covs], [p[0] for p in families_covs])
        # Median-normalization
        # median_normalized_covs[sample] = [c / median[sample] for c in sample2famcovlist[sample][0]]
        for cov in sample2famcovlist[sample][0]:
            normed_cov = 0.0
            if not median[sampleID] == 0:
                normed_cov = cov / median[sampleID]
            median_normalized_covs[sampleID].append(normed_cov)    
        # 
        for f in families:
            normed_cov = 0.0
            if not median[sampleID] == 0:
                normed_cov = samples_coverages[sample][f] / median[sampleID]
            norm_samples_coverages[sample][f] = normed_cov
        # samples_coverages[sample] = {f : samples_coverages[sample][f] / median[sample] for f in samples_coverages[sample]}

        # min coverage & plateau filter
        mediancov = median[sampleID] # self-defined median func: median of avg_genome_length, see above
        leftcov   = median_normalized_covs[sampleID][int(avg_genome_length * 0.3)]
        rightcov  = median_normalized_covs[sampleID][int(avg_genome_length * 0.7)]
        zerocov   = median_normalized_covs[sampleID][int(avg_genome_length * 1.25)]
        sample_stats[sample_id] = {'strainCoverage' : mediancov}
        if VERBOSE:
            print(' [I] ' + sample_id + ' median coverage: ' + str(round(mediancov,2)) +
                  '; left-side cov: ' + str(round(leftcov,2)) +
                  '; right-side cov: ' + str(round(rightcov,2)) +
                  '; out-plateau cov: ' + str(round(zerocov,2)) )
        sample2accepted[sample] = True if mediancov >= th_min_coverage else False # min coverage filter
        if not sample2accepted[sample]:
            print('     ' + sample_id + ': no strain detected, sample below MIN COVERAGE threshold')    
        if sample2accepted[sample]: # check left right plateau coverage
            if leftcov > th_plateau_left_max:
                sample2accepted[sample] = False
                if VERBOSE:
                    print('     ' + sample_id + ': no strain detected, sample does not pass LEFT-side coverage threshold.')
            elif rightcov < th_plateau_right_min:
                sample2accepted[sample] = False
                if VERBOSE:
                    print('     ' + sample_id + ': no strain detected, sample does not pass RIGHT-side coverage threshold.')        
        if sample2accepted[sample]:
            sample_stats[sample_id].update({'strainIsPresent' : True})
            sample_stats[sample_id].update({'Multistrain' : zerocov > th_max_zero})
            if VERBOSE: print('     ' + sample_id + ' OK - strain detected')
            if zerocov > th_max_zero:
                # to do: add to dict
                if VERBOSE: print('     ' + sample_id + ' WARNING: sample may contain multiple strains')
        else:
            sample_stats[sample_id].update({'strainIsPresent' : False})
            sample_stats[sample_id].update({'Multistrain' : False})
    accepted_samples_list = sorted([s for s in sample2accepted if sample2accepted[s]])
    return sample2accepted, accepted_samples_list, norm_samples_coverages, sample2famcovlist, sample2color, median_normalized_covs, median, sample_stats
# ----------------------------------------------------------------------------------------------------
def plot_dna_coverage(sample2accepted, samples_coverages, sample2famcovlist, sample2color, median_normalized_covs, genome_length, clade, plot1_name, plot2_name, INTERACTIVE, TIME, VERBOSE=False):
    '''
    Plot gene-family coverage plots as pdf file.
    a) absolute coverage
    b) median normalized coverage
    Accepted samples are plotted in colors, rejected samples in gray.
    '''
    try:
        if not INTERACTIVE:        # save to file
            import matplotlib      # for non-interactive plots on server without X11
            matplotlib.use('Agg')  # set 'Agg' before import pylab
        import matplotlib.pyplot as plt    
        try:
            from pylab import legend, savefig

            samples = sorted(samples_coverages.keys())
            accepted2samples = defaultdict(list)
            for s in samples:
                if sample2accepted[s]:
                    accepted2samples[True].append(s)
                else:
                    accepted2samples[False].append(s)
            sorted_samples = accepted2samples[False]
            sorted_samples.extend(accepted2samples[True])
            num_accepted=len(accepted2samples[True])

            # Plotting...
            if not plot1_name == '' or not plot2_name == '':
                if VERBOSE:
                    print(' [I] Plot gene-family coverage curves')
                used_colors = []
                for sample in accepted2samples[True]:
                    color, reset = random_color(used_colors)
                    sample2color[sample] = color
                    if reset:
                        used_colors = [color]
                    else:
                        used_colors.append(color)

                # Absolute coverage plot
                fig1 = None
                if not plot1_name == '':
                    plt.suptitle('Gene families coverages')
                    plt.xlabel('Gene families')
                    plt.ylabel('Coverage')
    
                    for sample in sorted_samples:
                        sample_id = get_sampleID_from_path(sample, clade)
                        covs = sample2famcovlist[sample][0]
                        if sample2accepted[sample]:
                            plt.plot(range(1, len(covs) + 1), covs, sample2color[sample], label=sample_id)
                        else:
                            plt.plot(range(1, len(covs) + 1), covs, COLOR_GREY)
                    plt.axis([0.0, genome_length * 1.5, 0.0, 1000.0])
                    try:
                        if num_accepted > 0: plt.legend(loc='upper right', fontsize='xx-small')
                    except TypeError:
                        print(' [W] pylab.legend fontsize does not work (please update your "pylab" module version).')  
                    savefig(plot1_name) 
                    if INTERACTIVE:
                        fig1 = plt.figure(0)
                    plt.close()

                # Median-normalized coverage plot
                fig2 = None
                if not plot2_name == '':
                    plt.suptitle('Gene families normalized coverages')
                    plt.xlabel('Gene families')
                    plt.ylabel('Normalized coverage')
                    for sample in sorted_samples:
                        sample_id = get_sampleID_from_path(sample, clade)
                        covs = median_normalized_covs[sample_id]
                        if sample2accepted[sample]:
                            plt.plot(range(1, len(covs) + 1), covs, sample2color[sample], label=sample_id)
                        else:
                            plt.plot(range(1, len(covs) + 1), covs, COLOR_GREY)
                    plt.axis([0.0, genome_length * 1.5, 0.0, 9.0])
                    
                    # thresholds need to adapted from dna_sample_filtering
                    # plt.plot((0.0, genome_length * 1.5), (th_present, th_present), 'k--') # th_present horizontal
                    # plt.plot((0.9 * genome_length, 0.9 * genome_length), (0.0, 9.0), 'k--') # genome length lowerbound vertical
                    # plt.plot((1.1 * genome_length, 1.1 * genome_length), (0.0, 9.0), 'k--') # genome length upperbound vertical
                    # plt.plot([genome_length * 0.3], [left_max], 'ro') # left_max intersected with genome length * 0.3
                    # plt.plot([genome_length * 0.7], [right_min], 'ro') # right_min intersected with genome length * 0.7

                    try:
                        if num_accepted > 0: plt.legend(loc='upper right', fontsize='xx-small')
                    except TypeError:
                        print(' [W] pylab.legend fontsize does not work (please update your "pylab" module version).')
                    savefig(plot2_name)
                    if INTERACTIVE:
                        fig2 = plt.plot()

            del(samples)
            del(accepted2samples)
            return True

        except ImportError:
            print(' [W] "pylab" module is not installed.')
            print('     To visualize and save charts, you need both "matplotlib" and "pylab" modules.')
    except ImportError:
        print(' [W] "matplotlib" module is not installed.')    
    return False
# -----------------------------------------------------------------------------
def print_coverage_matrix(dna_files_list, dna_file2id, dna_samples_covs_path, out_channel, families, clade, TIME, VERBOSE):
    '''
    TODO
    '''
    dna_sample_ids = sorted([dna_file2id[s] for s in dna_files_list])
    id2file = dict((v,k) for (k,v) in dna_file2id.items())
    if not out_channel == '':
        with open(out_channel, mode='w') as csv:
            csv.write('\t' + '\t'.join([get_sampleID_from_path(s, clade) for s in dna_sample_ids]) + '\n')
            for f in families:
                if sum(dna_samples_covs_path[s][f] for s in dna_samples_covs_path) > 0.0:
                    csv.write(f)
                    for s in dna_sample_ids:
                        csv.write('\t' + str(format(dna_samples_covs_path[id2file[s]][f], '.3f')))
                    csv.write('\n')
        if VERBOSE:
            TIME = time_message(TIME, 'Gene families coverage matrix has been printed in ' + out_channel + ' -')
    return TIME
# -----------------------------------------------------------------------------
def families_coverages(gene2cov, gene2family, lengths, VERBOSE):
    '''
    Compute the gene families coverages clustering the genes coverages
    '''
    # fami_covs = { GENE FAMILY : ( SUM OF FAMILY'S GENE UNNORMALIZED COVERAGES , [ GENE'S LENGTH ] ) }
    family2cov = defaultdict(list)

    for g in gene2family:
        if g in gene2cov:
            family2cov[gene2family[g]].append((gene2cov[g], lengths[g]))
        else:
            family2cov[gene2family[g]].append((0, lengths[g]))

    for f in family2cov:
        # family2cov[f] := [(cov1, len1), (cov2, len2), ...]
        sum_of_covs = float(sum(e[0] for e in family2cov[f]))
        sum_of_lens = float(sum(e[1] for e in family2cov[f]))
        cov = sum_of_covs / (sum_of_lens / len(family2cov[f]))
        family2cov[f] = cov

    return family2cov
# -----------------------------------------------------------------------------
def read_pangenome(panphlan_clade, VERBOSE):
    '''
    1) Search pangenome file, exit if not present
    2) Build the following data structures:
     - (dict) length for each gene
     - (dict) family for each gene_coverages
     - (list) sorted list of family
     - (int) average length of the genomes (in terms of # of gene)
    '''
    gene_lengths = {}
    gene2family = {}
    families = set()
    genome_lengths = defaultdict(int)
    genome2families = defaultdict(set)
    
    # search pangenome file: first in local pwd; second in $BOWTIE2_INDEXES
    clade=panphlan_clade.replace('panphlan_','')
    filename = 'panphlan_' + clade + '_pangenome.csv'
    path_local  = os.path.join(os.getcwd(),filename)
    path_bowtie = os.path.join(os.environ['BOWTIE2_INDEXES'],filename)
    if os.path.exists(path_local):
        pangenome_file=path_local
        if VERBOSE: print(' [I] Pangenome file: ./' + filename)
    elif os.path.exists(path_bowtie):
        pangenome_file = path_bowtie
        if VERBOSE: print(' [I] Pangenome file: ' + pangenome_file)
    else:
        show_error_message('Error: pangenome file "' + filename + '" not found.')
        sys.exit(INEXISTENCE_ERROR_CODE)
     
    # read pangnome file
    with open(pangenome_file, mode='r') as f:
        for line in f:
            words = line.strip().split('\t')
            fml, gene, genome, ctg, fr, to = words[FAMILY_INDEX], words[GENE_INDEX], words[GENOME_INDEX], words[CONTIG_INDEX], int(words[FROM_INDEX]), int(words[TO_INDEX])
            gene_lengths[gene] = abs(to - fr) + 1
            gene2family[gene] = fml
            families.add(fml)
            genome2families[genome].add(fml)
    
    # Get expected median genome length (number of gene families)
    genome_lengths    = dict((g, len(genome2families[g])) for g in genome2families)
    num_ref_genomes   = len(genome_lengths)
    avg_genome_length = int(numpy.median(list(genome_lengths.values())))
    ref_genomes      = sorted(genome2families.keys())
    
    if VERBOSE:
        print('     Number of reference genomes: '                + str(num_ref_genomes))
        print('     Average number of gene-families per genome: ' + str(avg_genome_length))
        print('     Total number of pangenome gene-families '     + str(len(families)))
    return gene_lengths, gene2family, sorted(list(families)), num_ref_genomes, avg_genome_length, genome2families, ref_genomes
# -----------------------------------------------------------------------------
def read_gene_cov_file(input_file):
    '''
    Put the information contained in a file into a dictionary data structure
    '''
    import bz2
    d = {}
    f = bz2.BZ2File(input_file, mode='r')
    for line in f:
        words = line.decode('utf-8').strip().split('\t')
        gene, coverage = words[0], int(words[1])
        d[gene] = coverage
    f.close()
    return d
# -----------------------------------------------------------------------------
def read_map_results(i_dna, i_rna, clade, RNASEQ, VERBOSE):
    '''
    Read results from panphlan_map.py
    '''
    dna_samples_covs = {}
    dna_samples_covs_path = {} # old Thomas version, path as key
    rna_samples_covs = {}
    rna_samples_covs_path = {} # old Thomas version, path as key
    dna_files_list = []
    rna_id_list    = []
    if not i_dna == None:
        dna_files_list = sorted(i_dna.keys())
        for dna_covs_file in dna_files_list:
            sample_id = get_sampleID_from_path(dna_covs_file, clade)
            dna_samples_covs[sample_id] = read_gene_cov_file(dna_covs_file) # new dict
            dna_samples_covs_path[dna_covs_file] = dna_samples_covs[sample_id] # old dict (path as key)
        if RNASEQ:
            rna_id_list = sorted(i_rna.keys())
            for rna_covs_id in rna_id_list:
                rna_covs_file = i_rna[rna_covs_id]
                if not rna_covs_file == NO_RNA_FILE_KEY:
                    rna_samples_covs_path[rna_covs_file] = read_gene_cov_file(rna_covs_file)
    return dna_samples_covs, dna_samples_covs_path, dna_files_list, rna_samples_covs, rna_samples_covs_path, rna_id_list
# -----------------------------------------------------------------------------
def check_args():
    '''
    Check if the input arguments respect the rules of usage

        Usage examples:
            panphlan_join.py -c ecoli -i sampleCSVdiectory  > projectName_saureus_pangenome_coverage.csv
            panphlan_join.py -c ecoli -i sampleCSVdiectory -o projectName_saureus_pangenome_coverage.csv
            panphlan_join.py -c sepidermidis --i_dna DNAdir/ --i_rna RNAdir/ --sample_pairs DNA_RNA_sampleIDs.csv -o HMP_sepidermidis_RNAseq_gene_expression.csv
    '''
    parser = PanPhlAnJoinParser()
    args = vars(parser.parse_args())

    VERBOSE = args['verbose']

    if VERBOSE:
        print('\nPanPhlAn profile version '+__version__)
        print('Python version: ' + sys.version.split()[0])
        print('System: ' + sys.platform)

    # Check CLADE
    clade = args['clade']
    if not clade.startswith(PANPHLAN):
        args['clade'] = PANPHLAN + clade
    if VERBOSE:
        print('[I] Clade: ' + args['clade'])

    # Check DNA_RNA_MAPPING
    pairs_path = args['sample_pairs']
    idna = args['i_dna']
    irna = args['i_rna']
    # dna2rna := { DNA_ID : RNA_ID }
    dna2rna = {}

    if idna == '': # if no samples, print presence/absence of reference genomes
        args['i_dna'] = None
        args['add_strains'] = True
    else:
        # Normal pipeline
        if not pairs_path == None:
            # --sample_pairs is defined: check if the file exists or not
            if not os.path.exists(pairs_path):
                show_error_message('DNA-RNA mapping file does not exist.')
                sys.exit(INEXISTENCE_ERROR_CODE)
            else:
                if idna == None or irna == None:
                    show_error_message('With option --sample_pairs must be defined also options --i_dna and --i_rna.')
                    sys.exit(PARAMETER_ERROR_CODE)
                else:
                    # --i_dna and --i_rna are defined: check if the folders exist or not
                    if not os.path.exists(idna):
                        show_error_message('Input folder for DNA files does not exist.')
                        sys.exit(INEXISTENCE_ERROR_CODE)
                    if not os.path.exists(irna):
                        show_error_message('Input folder for RNA files does not exist.')
                        sys.exit(INEXISTENCE_ERROR_CODE)

                    # --sample_pairs, --i_dna, --i_rna are all defined
                    try: # read file of DNA RNA sample pairs
                        for line in (l.strip() for l in open(pairs_path) if not l.startswith('#')):
                            words = line.split('\t')
                            dna2rna[words[0]] = words[1]
                    except IndexError as err:
                        show_error_message(err)
                        print(line)
                        print('[E] Cannot read DNA RNA --sample_pairs:')
                        print('    File format needs to follow tab-separated pairs of DNA and RNA sample names.')
                        sys.exit(FILEFORMAT_ERROR_CODE)
                    # Search for files reading the DNA-RNA mapping
                    dna_file2id = {}
                    rna_id2file = {}
                    # NB. The idea is: rna_id2file[dna2rna[dna_file2id[sample_dna_file_name]]]
                    #     Also because if we have not defined --sample_pairs, we still have the same dict structure for i_dna[COVERAGES_KEY]
                    #     (i.e. {DNA file path : DNA sample id} <==> {DNA file path : None}) accesing to its values with i_dna[COVERAGES_KEY].keys()
                    for d in sorted([s for s in dna2rna.keys()]):
                        dna_path = find('*' + d + '*.csv.bz2', idna)
                        if dna_path == []:
                            print('[W] DNA file corresponding to ID ' + d + ' has not been found. Analysis and mapping for this DNA will be skipped.')
                            continue
                        dna_file2id[dna_path[0]] = d

                        rna_path = find('*' + dna2rna[d] + '*.csv.bz2', irna)
                        if rna_path == []:
                            print('[W] RNA file corresponding to ID ' + dna2rna[d] + ' has not been found. Analysis for this RNA will be skipped.')
                            rna_id2file[dna2rna[d]] = NO_RNA_FILE_KEY
                        else:
                            rna_id2file[dna2rna[d]] = rna_path[0]
                    args['sample_pairs'] = dna2rna
                    args['i_dna'] = dna_file2id
                    args['i_rna'] = rna_id2file
                    if VERBOSE:
                        print('[I] Input folder for DNAs: ' + idna)
                        print('[I] Input folder for RNAs: ' + irna)
                        print('[I] Gene coverages files:\n\t' + '\n\t'.join(sorted(list(dna_file2id.keys()))))
                        print('[I] Trascripts coverages files:\n\t' + '\n\t'.join(sorted(list(rna_id2file.values()))))
                        print('[I] DNA-RNA projects mapping:\n\t' + '\n\t'.join(k+' >>> '+v for k,v in sorted(dna2rna.items())))
        else:
            # --sample_pairs is not defined: check only --i_dna
            
            if not irna == None: # If --sample_pairs is NOT defined BUT --i_rna yes, then error
                show_error_message('Option --sample_pairs has not been defined, but --i_rna is defined. You must decide if define both or no one of them.')
                sys.exit(PARAMETER_ERROR_CODE)

            if not os.path.exists(idna):
                show_error_message('Input folder for DNA files does not exist.')
                sys.exit(INEXISTENCE_ERROR_CODE)
            
            # Find coverages file
            covs_file_pattern = '*' + args['clade'].replace('panphlan_', '') + '*.csv.bz2'
            if VERBOSE:
                print('[I] Looking for "' + covs_file_pattern + '"-patterned files...')
            covs_files = find(covs_file_pattern, idna)
            for f in covs_files:
                # In the (remote) case where the pangenome file is zipped (.csv.bz2) and located in the same folder of the DNA abundance files, then delete it from the list
                if 'pangenome' in f:
                    covs_files.pop(covs_files.index(f))
            if covs_files == []:
                show_error_message('Any gene coverages file has not been found.')
                sys.exit(INEXISTENCE_ERROR_CODE)
            if VERBOSE:
                print('[I] Found ' + str(len(covs_files)) + ' abundances files.')
            samples_files = dict((f, get_sampleID_from_path(f, args['clade'])) for f in covs_files)

            # TODO choose only one pangenome file if more than one are found

            args['i_dna'] = samples_files            
            if VERBOSE:
                print('[I] Input folder: ' + idna)
                print('[I] Gene coverages files:\n\t' + '\n\t'.join(sorted(covs_files)))

    # Check OUTPUT_FILE
    args['o_dna'] = check_output(args['o_dna'], '', 'gene families presence/absence matrix', VERBOSE)

    # Check COVERAGE_OUT_CSV
    args['o_cov'] = check_output(args['o_cov'], '', 'gene families normalized coverage', VERBOSE)

    # Check DNA_INDEX_FILE
    args['o_idx'] = check_output(args['o_idx'], '', 'gene families DNA indexing', VERBOSE)
    
    # Check RNA_EXPRS_FILE
    args['o_rna'] = check_output(args['o_rna'], '', 'transcript families normalized coverage', VERBOSE) 

    # Check GENEHIT_PERC_PER_STRAIN
    args['strain_hit_genes_perc'] = check_output(args['strain_hit_genes_perc'], '', 'strains gene hit percentage', VERBOSE) 

    # Check MINIMUM_THRESHOLD, MEDIUM_THRESHOLD and MAXIMUM_THRESHOLD
    a, b, c = args['th_zero'], args['th_present'], args['th_multicopy']
    ok = False
    if b != None:
        if b < MIN_PRESENT_TH:
            ok = False
        elif a != None or b != None:
            if a != None and b != None:
                if a < MIN_NONPRESENT_TH or c < MIN_MULTICOPY_TH:
                    ok = False
                elif a < b and b < c:
                    # Unique case where a and c are not automatically set
                    ok = True
            else:
                a = b / 2.0
                c = 3.0 * b
                ok = True
        else:
            # Both re not set
            a = b / 2.0
            c = 3.0 * b
            ok = True
    else:
        if a == None and b == None:
            b = PRESENT_TH
            a = b / 2.0
            c = 3.0 * b
            ok = True
        else:
            ok = False

    args['th_zero'], args['th_present'], args['th_multicopy'] = a, b, c
    if VERBOSE:
        print('[I] Non-presence threshold: ' + str(args['th_zero']))
        print('[I] Presence threshold: ' + str(args['th_present']))
        print('[I] Multicopy threshold: ' + str(args['th_multicopy']))

    if not ok:
        show_error_message('Thresholds are set to unacceptable values.')
        if VERBOSE:
            print('    Please, follow this usage: [--th_present B [--th_zero A --th_multicopy C]]\n    with A < B < C. Default values are A = 0.25, B = 0.50, C = 1.50')
        sys.exit(PARAMETER_ERROR_CODE)        
        
    # Check strain filter (coverage curve) LEFT_MAX, RIGHT_MIN, MIN_COVERAGE_MEDIAN
    if args['left_max'] is not None:
        if args['left_max'] <= 1:
            show_error_message('Threshold left_max must be greater than 1.')
            sys.exit(PARAMETER_ERROR_CODE)
    if args['right_min'] is not None:
        if args['right_min'] >= 1:
            show_error_message('Threshold right_min must be smaller than 1.')
            sys.exit(PARAMETER_ERROR_CODE)    
    if args['min_coverage'] is not None:
        if args['min_coverage'] < 0.0:
            if VERBOSE:
                print('[W] Unacceptable value for minimum median coverage threshold: ' + str(args['min_coverage']) +'. Set default: ' + str(COVERAGE_TH))
            args['min_coverage'] = COVERAGE_TH

    # Check RNA_MAX_ZEROES
    if args['rna_max_zeros'] < 0.0 or args['rna_max_zeros'] > 100.0:
        args['rna_max_zeros'] = RNA_MAX_ZERO_TH
        if VERBOSE:
            print('[W] Unacceptable value for RNA maximum zeros threshold. Set default.')
    if VERBOSE:
        print('[I] RNA maximum zeros threshold: ' + str(args['rna_max_zeros']))

    # Check SIMILARITY_PERCENTAGE
    if args['strain_similarity_perc'] < 0.0 or args['strain_similarity_perc'] > 100.0:
        args['strain_similarity_perc'] = SIMILARITY_TH
        if VERBOSE:
            print('[W] Unacceptable value for strain similiarity percentage threshold. Set default.')
    if VERBOSE:
        print('[I] Strain similiarity percentage threshold: ' + str(args['strain_similarity_perc']))

    # Check NON_PRESENCE_TOKEN
    if args['np'] in UNACCEPTABLE_NP:
        args['np'] = DEFAULT_NP
        if VERBOSE:
            print('[W] Unacceptable string for non-presence (absence) token. Set default.')
    if VERBOSE:
        print('[I] Non-presence token: ' + str(args['np']))

    # Check NOT_A_NUMBER_TOKEN
    if args['nan'] in UNACCEPTABLE_NAN:
        args['nan'] = DEFAULT_NAN
        if VERBOSE:
            print('[W] Unacceptable string for NaN token. Set default.')
    if VERBOSE:
        print('[I] NaN token: ' + str(args['nan']))
    
    # Check COV_PLOT_NAME
    args['o_covplot'] = check_output(args['o_covplot'], '', 'gene coverage plot', VERBOSE) # Will never be never equals to None, so we don't need a default value

    # Check NOR_PLOT_NAME
    args['o_covplot_normed'] = check_output(args['o_covplot_normed'], '', 'gene normalized coverage plot', VERBOSE) # Will never be never equals to None, so we don't need a default value

    return args
# -----------------------------------------------------------------------------

def main():
    # Check Python version
    if sys.hexversion < 0x02060000:
        print('Python version: ' + sys.version)
        sys.exit('Python versions older than 2.6 are not supported.')

    args = check_args()    
        
    TOTAL_TIME = time.time()
    TIME       = time.time()

    VERBOSE     = args['verbose']
    INTERACTIVE = args['interactive']
    ADD_STRAINS = args['add_strains']
    RNASEQ      = True if args['sample_pairs'] else False

    # Create pangenome dicts: gene->family, genome->families, gene->length
    if VERBOSE: print('\nSTEP 1. Read pangenome data...')
    gene_lenghts, gene2family, families, num_ref_genomes, avg_genome_length, genome2families, ref_genomes = read_pangenome(args['clade'], VERBOSE)

    # read mapping result files
    if VERBOSE: print('\nSTEP 2. Read mapping results ...')
    dna_samples_covs, dna_samples_covs_path, dna_files_list, rna_samples_covs, rna_samples_covs_path, rna_id_list = read_map_results(
        args['i_dna'], args['i_rna'], args['clade'], RNASEQ, VERBOSE)

    # Presence/absence matrix only of reference genomes, no samples
    if ADD_STRAINS or args['strain_hit_genes_perc'] != '':
        if VERBOSE: print('\nSTEP 3a. Extracting reference genomes gene repertoire...')
        TIME, strain2family2presence = build_strain2family2presence(ref_genomes, families, genome2families, TIME, VERBOSE)
        if ADD_STRAINS and args['i_dna'] == None:
            if VERBOSE: print('\nSTEP 3b. Printing presence/absence binary matrix only for reference genomes...')
            TIME = print_presence_absence_profiles(ref_genomes, strain2family2presence, families, args['o_dna'], TIME, VERBOSE)
            end_program(time.time() - TOTAL_TIME)
            sys.exit(0) 

    # Merge gene/transcript abundance into family (normalized) coverage
    if VERBOSE: print('\nSTEP 3. Merge single gene abundances to gene family coverages')
    for sample in dna_files_list:
        if VERBOSE: print(' [I] Normalization for DNA sample ' + get_sampleID_from_path(sample, args['clade']) + '...')
        dna_samples_covs_path[sample] = families_coverages(dna_samples_covs_path[sample], gene2family, gene_lenghts, VERBOSE)
    
    # Get samples list
    # Print coverages in file
    TIME = print_coverage_matrix(dna_files_list, args['i_dna'], dna_samples_covs_path, args['o_cov'], families, args['clade'], TIME, VERBOSE)
    

    # Filter DNA samples according to their median coverage value and plot coverage plateau
    if VERBOSE:
        print('\nSTEP 4: Strain presence/absence filter based on coverage plateau curve...')
    sample2accepted, accepted_samples, norm_dna_samples_covs_path, sample2famcovlist, sample2color, median_normalized_covs, sample2median, sample_stats = strain_presence_plateau_filter(dna_samples_covs_path, num_ref_genomes, avg_genome_length, args['min_coverage'], args['left_max'], args['right_min'], families, args['clade'], TIME, VERBOSE)
    result = plot_dna_coverage(sample2accepted, norm_dna_samples_covs_path, sample2famcovlist, sample2color, median_normalized_covs, avg_genome_length, args['clade'], args['o_covplot'], args['o_covplot_normed'], INTERACTIVE, TIME, VERBOSE)


    # DNA indexing
    if VERBOSE: print('\nSTEP 5a: Define multicopy, strain-specific, and non-present gene-families (1,-1,-2,-3 matrix, option --o_idx)')
    sample2family2dnaidx, TIME = get_idx123_plateau_definitions(accepted_samples, norm_dna_samples_covs_path, args['th_zero'], args['th_present'], args['th_multicopy'], args['o_idx'], families, args['clade'], TIME, VERBOSE)
    if VERBOSE: print('\nSTEP 5b: Get presence/absence of gene-families (1,-1 matrix, option --o_dna)')
    dna_sample2family2presence, sample_stats, TIME = get_genefamily_presence_absence(sample2accepted, dna_files_list, args['i_dna'], sample2family2dnaidx, args['o_dna'], families, args['clade'], avg_genome_length, sample_stats, TIME, VERBOSE)

    if ADD_STRAINS or args['strain_hit_genes_perc'] != '':
        if VERBOSE: print('\nSTEP 5c: Calculate percent of identical gene-families between sample-strains and reference-genomes... (option --strain_hit_genes_perc)')
        ss_presence, TIME = samples_strains_presences(dna_sample2family2presence, ref_genomes, strain2family2presence, avg_genome_length, args['strain_similarity_perc'], args['o_dna'], families, args['clade'], ADD_STRAINS, TIME, VERBOSE)
        # ss_presence = { STRAIN or SAMPLE : { GENE FAMILY : PRESENCE } }
        if args['strain_hit_genes_perc'] != '':
            if VERBOSE: print('\nSTEP 5d: Get percent of sample-strain gene-families present in reference-genomes (option --strain_hit_genes_perc)')
            strain2sample2hit, TIME = strains_gene_hit_percentage(ss_presence, genome2families, sample2accepted, args['strain_hit_genes_perc'], args['clade'], TIME, VERBOSE)


    # RNA-seq: get RNA gene-family coverage
    rna_file_list = []
    if RNASEQ:
        if VERBOSE:
            print('\nSTEP 6. RNA-seq: Merge single gene transcript abundances to gene-family transcript coverages')
        for sample_id in rna_id_list:
            sample = args['i_rna'][sample_id]
            if not sample == NO_RNA_FILE_KEY:
                if VERBOSE:
                    print('[I] Normalization for RNA sample ' + get_sampleID_from_path(sample, args['clade']) + '...')
                rna_file_list.append(sample)
                rna_samples_covs_path[sample] = families_coverages(rna_samples_covs_path[sample], gene2family, gene_lenghts, VERBOSE)


    # DNA (and RNA) indexing
    if RNASEQ:
        if VERBOSE: print('\nSTEP 8. RNA-seq: Get strain-specific gene transcription profiles')
        rna_seq(args['o_rna'], sample2family2dnaidx, dna_samples_covs_path, sample2accepted, rna_id_list, rna_samples_covs_path, args['rna_max_zeros'], args['sample_pairs'], args['i_dna'], args['i_rna'], families, args['np'], args['nan'], args['clade'], args['rna_norm_percentile'], TIME, VERBOSE)

    end_program(time.time() - TOTAL_TIME) 

    # WARNING for presence of multiple strains in a sample 
    print_multistrain_warning(sample_stats, avg_genome_length, VERBOSE)


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
