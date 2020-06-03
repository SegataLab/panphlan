#!/usr/bin/env python

"""
panphlan_profiling.py
    Concatenate multiple output files from panphlan_map.py to create a presence/absence matrix of gene families across metagenomic samples.
"""

import os, subprocess, sys, time, bz2
import numpy
import argparse as ap
from collections import defaultdict
from shutil import copyfileobj
from misc import random_color
from random import randint


__author__ = 'Leonard Dubois, Matthias Scholz, Thomas Tolio and Nicola Segata (contact on https://forum.biobakery.org/)'
__version__ = '3.0'
__date__ = '20 April 2020'


# FIXED THRESHOLDS
ZERO_NON_PLATEAU_TH = 0.20 # multistrain detection
#  + Arguments default values

# ------------------------------------------------------------------------------
"""
Reads and parses the command line arguments of the script.
:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    # INPUT ARGUMENTS
    p.add_argument('-i','--i_dna', type=str, default=None,
                   help='Input directory of panphlan_map.py results, containing SAMPLE.csv.bz2 files')
    p.add_argument('-p', '--pangenome', type = str,
                   help='Path to pangenome tsv file exported from ChocoPhlAn')
    p.add_argument('--i_covmat', type=str, default=None,
                   help='Path to precomputed coverage matrix')               
                 
    # OUTPUT ARGUMENTS
    p.add_argument('--o_matrix', type=str, default=None,
                   help='Path for presence/absence matrix output')
    p.add_argument('--o_covmat', type=str, default=None,
                   help='Write raw gene-family coverage matrix in provided file')
    p.add_argument('--o_covplot_normed', type=str, default=None,
                   help='Filename for normalized gene-family coverage plot.')
    p.add_argument('--o_idx', metavar='DNA_INDEX_FILE', type=str, default= None,
                   help='Write gene-family plateau definitions (1, -1, -2, -3)')

    # THRESHOLDS TUNING ARGUMENTS
    # Strains (samples) presence/absence threshold
    p.add_argument('--min_coverage', type=float, default=2,
                   help='Minimum coverage threshold, default: 2X')
    p.add_argument('--left_max', type=float, default=1.25, # v1.0: 1.18 strain presence/absence filter (plateau curve)
                   help='Strain presence/absence plateau curve threshold: left max [1.25]')
    p.add_argument('--right_min', type=float, default=0.75, # v1.0: 0.82
                   help='Strain presence/absence plateau curve threshold: right min [0.75]')
    # Gene families presence/absence threshold
    p.add_argument('--th_non_present', type=float, default = 0.05,
                   help='Gene families threshold: not present if lower')
    p.add_argument('--th_present', type=float, default = 0.5,
                   help='Gene families threshold: present if higher')
    p.add_argument('--th_multicopy', type=float, default = 0.15,
                   help='Gene families threshold: multicopy if higher')
    # filter similarity
    p.add_argument('--strain_similarity_perc', metavar='SIMILARITY_PERCENTAGE', type=float, default=50.0,
                    help='Minimum threshold (percentage) for genome length to add a reference genome to presence/absence matrix (default: 50).')

    # OPTIONAL ARGUMENTS
    p.add_argument('--add_ref', action='store_true',
                   help='Add reference genomes to gene-family presence/absence matrix.')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='Show progress information')
    # FUNCTIONNAL ANNOTATION ARGUMENTS
    p.add_argument('--func_annot', type=str, default=None,
                   help='Path to file mapping UniRef IDs to GO/KEGG/... annotation for functional characterization')
    p.add_argument('-f', '--field', type=int, default=1,
                   help='Field in the annotation file that must be added to the presence/absence matrix')

    return p.parse_args()


"""Check arguments consistency"""
def check_args(args):
    if args.i_dna:
        if not os.path.exists(args.i_dna):
            sys.exit('[E] Sample file directory (' + args.i_dna + ') not found\n')
    else:
        sys.exit('[E] Please provide a valid sample file (argument -i or --i_dna).\n')



# ------------------------------------------------------------------------------
#   STEP 1
# ------------------------------------------------------------------------------
def read_pangenome(pangenome_file):
    """Build the following data structures:
     - (dict) length and family for each gene
     - (list) sorted list of family
     - (int) families in each genome
    Other informations can be extracted from these
    """
    genes_info = defaultdict(dict)
    families = set()
    genome_lengths = defaultdict(int)
    genome2families = defaultdict(set)

    pangenome_file = pangenome_file
    with open(pangenome_file, mode='r') as f:
        for line in f:
            words = line.strip().split('\t')
            fml, gene, genome, ctg, fr, to = words[0], words[1], words[2], words[3], int(words[4]), int(words[5])
            genes_info[gene] = {'length' : abs(to - fr) + 1, 'family' : fml}
            families.add(fml)
            if not genome.startswith('REF_'):
                genome = 'REF_' + genome
            genome2families[genome].add(fml)

    # Get expected median genome length (number of gene families)
    genome_lengths    = dict((g, len(genome2families[g])) for g in genome2families)
    num_ref_genomes   = len(genome_lengths)
    avg_genome_length = int(numpy.median(list(genome_lengths.values())))
    print('     Number of reference genomes: '                + str(num_ref_genomes))
    print('     Average number of gene-families per genome: ' + str(avg_genome_length))
    print('     Total number of pangenome gene-families '     + str(len(families)))
    return genes_info, sorted(list(families)), genome2families

# ------------------------------------------------------------------------------
#   STEP 1 BIS
# ------------------------------------------------------------------------------
def build_ref2family2presence(families, genome2families, VERBOSE):
    """Build the dictionary from ref to gene family to presence
    { REF : { GENE FAMILY : PRESENCE(True or False) } }
    """
    ref_genomes = sorted(genome2families.keys())
    ref2family2presence = defaultdict(dict)
    numof_ref = len(ref_genomes)
    i = 1
    for s in ref_genomes:
        if VERBOSE:
            print('[I] [' + str(i) + '/' + str(numof_ref) + '] Analysing reference genome ' + s + '...')
            i += 1
        for f in families:
            ref2family2presence[s][f] = f in genome2families[s]
    if VERBOSE:
        print('Gene families presence/absence in reference genomes computed.')
    return  ref2family2presence

# ------------------------------------------------------------------------------
#   MATRIX OUTPUT
# ------------------------------------------------------------------------------
def filter_never_present(presences_dict, args):
    """Remove gene families never present.
    Also remove those which are only present in the samples (because some ref strain has been filtered out )
    """
    sample_and_strains = sorted(presences_dict.keys())
    families = sorted(presences_dict[sample_and_strains[0]].keys())
    never_present_families = []
    for f in families:
        always_zero = True
        for s in presences_dict:
            if presences_dict[s][f]: # == True == 1
                always_zero = False
                break
        if always_zero:
            never_present_families.append(f)
    if args.verbose:
        print(' [I] '+ str(len(never_present_families)) + ' never present gene families filtered out.')
    return never_present_families

def write_presence_absence_matrix(presences_dict, args, family2annot):
    """Function writing the presence/absence matrix in csv file from
    a dict variable containing data about samples, strains or both.
    It can also add a annotation collumn
    """
    sample_and_strains = sorted(presences_dict.keys())
    families = sorted(presences_dict[sample_and_strains[0]].keys())
    never_present_families = filter_never_present(presences_dict, args)

    if len(sample_and_strains) > 0:
        if args.verbose: print(' [I] Print gene-family presence/absence matrix to: ' + args.o_matrix)
        OUT = open(args.o_matrix, mode='w')
        header = '\t'.join(sample_and_strains) + '\n'
        if not family2annot == None:
            header = '\t' + 'annotation' + '\t' + header
        else:
            header = '\t' + header
        OUT.write(header)

        for f in families:
            if f not in never_present_families:
                line = f
                if not family2annot == None:
                    if not str(family2annot[f]) == "" :
                        line = line + '\t' + str(family2annot[f])
                    else :
                        line = line + '\t' + "NA" + '\t'
                for s in sample_and_strains:
                    if presences_dict[s][f]:
                        line = line + '\t1'
                    else:
                        line = line + '\t0'
                OUT.write(line + '\n')
        OUT.close()

# ------------------------------------------------------------------------------
#  FUNCTIONNAL ANNOTATION
# ------------------------------------------------------------------------------
def create_annot_dict(presences_dict, args):
    """Build dict mapping families to annotation before writing presence/abscence matrix
    """
    sample_and_strains = sorted(presences_dict.keys())
    families = sorted(presences_dict[sample_and_strains[0]].keys())

    # if annot file provided is the same as pangenome file
    if args.func_annot == args.pangenome:
        pangenome_file  = os.path.join(os.getcwd(), args.pangenome)
        with file(pangenome_file) as f:
            line = f.readline()
        if len(line.split()) < 7:
            if args.verbose : print(' [I] No annotation data were found.\n No information about families will be added to the presence/absence matrix\n')
            return None
        else:
            # get the annotation from pangenome file
            family2annot = defaultdict(str)
            with file(pangenome_file) as IN:
                for line in IN:
                    ids = line.strip().split('\t')
                    #1st field uniref90 mapped to 6th one, annotation (UniRef50)
                    family2annot[ids[0]] = ids[7]
            return family2annot
    else:
        # read the file provided
        if args.verbose : print(' [I] Mapping families to annotation... This operation can take several minutes')
        family2annot = defaultdict(str)

        # Check file extension bz2
        if (args.func_annot).endswith('.bz2'):
            IN = bz2.BZ2File(args.func_annot, mode='r')
        else:
            IN = open(args.func_annot, 'r')

        line = IN.readline()
        for line in IN:
            ids = line.strip().split('\t')
            uniref90 = ids[0]
            if len(ids) < args.field :
                continue
            annot = ids[args.field - 1]
            if uniref90 in families:
                family2annot[uniref90] = annot
            else:
                family2annot[uniref90] = ""

        IN.close()
    return family2annot

# ------------------------------------------------------------------------------
#  STEP 2 Create coverage matrix
# ------------------------------------------------------------------------------
def get_sampleID_from_path(sample_path):
    # example: "path/to/mapping/result/ERR54632_ecoli14.csv.bz2" -> "ERR54632_ecoli14"
    sampleID = os.path.basename(sample_path)
    sampleID = sampleID.replace('_map.tsv.bz2','')
    return sampleID

def read_gene_cov_file(input_file):
    """Convert coverage mapping file into a dictionary data structure"""
    d = {}
    f = bz2.BZ2File(input_file, mode='r')
    for line in f:
        words = line.decode('utf-8').strip().split('\t')
        gene, coverage = words[0], int(words[1])
        d[gene] = coverage
    f.close()
    return d

def read_map_results(i_dna, VERBOSE):
    """Read results from panphlan_map.py"""
    dna_samples_covs = {}
    dna_files_list =  os.listdir(i_dna)
    for dna_covs_file in dna_files_list: # i_dna: path2id
        dna_sample_id = get_sampleID_from_path(dna_covs_file)
        if VERBOSE: print(' [I] Reading mapping result file: ' + dna_covs_file )
        dna_samples_covs[dna_sample_id] = read_gene_cov_file(os.path.join(i_dna, dna_covs_file))
    return dna_samples_covs

def get_genefamily_coverages(gene2cov, genes_info, VERBOSE):
    """Sum single gene coverage to gene-families coverages based on pangenome clustering"""
    # family2gene_info = { GENE FAMILY : ( SUM OF FAMILY'S GENE UNNORMALIZED COVERAGES , [ GENE'S LENGTH ] ) }
    family2gene_info = defaultdict(list)
    family2cov = defaultdict(float)
    for g in genes_info.keys():
        if g in gene2cov:
            family2gene_info[genes_info[g]['family']].append((gene2cov[g], genes_info[g]['length']))
        else:
            family2gene_info[genes_info[g]['family']].append((0, genes_info[g]['length']))
    for f in family2gene_info:
        # family2gene_info[f] := [(cov1, len1), (cov2, len2), ...]
        sum_of_covs = float(sum(e[0] for e in family2gene_info[f]))
        sum_of_lens = float(sum(e[1] for e in family2gene_info[f]))
        cov = sum_of_covs / (sum_of_lens / len(family2gene_info[f]))
        family2cov[f] = cov
    return family2cov

def print_coverage_matrix(dna_samples_covs, out_channel, families, VERBOSE):
    """Print merged table of gene-family coverage for all samples (option: --o_cov)"""
    dna_sample_ids = sorted(dna_samples_covs.keys())
    with open(out_channel, mode='w') as OUT:
        OUT.write('\t' + '\t'.join(dna_sample_ids) + '\n')
        for f in families:
            tmp = [dna_samples_covs[s][f] for s in dna_samples_covs]
            if len(tmp) > 0:
                if sum(tmp) > 0.0:
                    OUT.write(f)
                    for s in dna_sample_ids:
                        OUT.write('\t' + str(format(dna_samples_covs[s][f], '.3f')))
                    OUT.write('\n')
    if VERBOSE: print('Gene families coverage matrix has been printed in ' + out_channel)

# Or READ EXISTING COVERAGE MATRIX

def read_coverage_matrix(cov_matrix_file):
    """Read coverage matrix (option --o_cov) for re-analysis using other thresholds"""
    
    sample2family2cov = defaultdict(dict)

    if not os.path.exists(cov_matrix_file):
        sys.exit('\nERROR: Could not find --i_covmat input file: ' + cov_matrix_file)

    with open(cov_matrix_file, mode='r') as IN:
        sample_list = IN.readline().strip().split('\t') # get headerline
        for i,line in enumerate(IN):
            cols = line.strip().split('\t')
            genefamilyID = cols[0]
            coverage_values = cols[1:]
            if not len(sample_list)==len(coverage_values):
                print('[E] ERROR while reading --i_cov: coverage lines does not fit number of sampleIDs in headerline')
            for sample, cov_str in zip(sample_list, coverage_values):
                try:
                    cov = float(cov_str)
                except ValueError:
                    print('[E] ERROR while reading --i_covmat: Could not convert coverage value "'+ cov_str +'" to number, line:' + str(i))
                sample2family2cov[sample][genefamilyID] = cov

    return sample2family2cov

# ------------------------------------------------------------------------------
#  STEP 3 Strain presence/absence filter based on coverage plateau curve
# ------------------------------------------------------------------------------
def adjust_genome_length(genome2families):
    """reduce expected number of gene-families, in case of only 1,2 or 3 ref. genomes in DB"""
    num_ref_genomes = len(genome2families.keys())
    avg_genome_length = int(numpy.median([len(genome2families[g]) for g in genome2families]))

    orig_avg_genome_length = avg_genome_length
    if num_ref_genomes == 3:
        avg_genome_length = int(round(0.90 * avg_genome_length))
        print('[I] Decrease expected gene-families per sample strain to: ' + str(avg_genome_length) + ' (0.9*' + str(orig_avg_genome_length) +') due to only 3 ref. genomes in DB')
    elif num_ref_genomes == 2:
        avg_genome_length = int(round(0.85 * avg_genome_length))
        print('[I] Decrease expected gene-families per sample strain to: ' + str(avg_genome_length) + ' (0.85*' + str(orig_avg_genome_length) +') due to only 2 ref. genomes in DB)')
    elif num_ref_genomes == 1:
        avg_genome_length = int(round(0.75 * avg_genome_length))
        print('[I] Decrease expected gene-families per sample strain to: ' + str(avg_genome_length) + ' (0.75*' + str(orig_avg_genome_length) +') due to only 1 ref. genomes in DB)')
    return avg_genome_length

def defining_normalized_coverage(samples_coverages, avg_genome_length, families):
    norm_samples_coverages = defaultdict(dict)
    median_cov = defaultdict()
    for sample in sorted(samples_coverages.keys()):
        d = samples_coverages[sample]
        # Take families coverage from sample and sort descendently by value (coverage)
        # items() turn dict to list of tuples (keys, values)
        families_covs = sorted(d.items(), key=lambda x: x[1])
        families_covs = families_covs[::-1] # reverse
        #del(d)
        median_cov[sample] = numpy.median([p[1] for p in families_covs][:avg_genome_length])
        for f in families: # all gene-families of the pangenome
            normed_cov = 0.0
            if (not median_cov[sample] == 0) and (f in samples_coverages[sample]):
                normed_cov = samples_coverages[sample][f] / median_cov[sample]
            norm_samples_coverages[sample][f] = normed_cov
    return norm_samples_coverages, median_cov

def strain_presence_plateau_filter(norm_samples_coverages, avg_genome_length, median_cov, args):
    """Check if a strain is present in a sample.
        Plateau quality criteria based on genes coverage curve.
        For each sample:
        - sort gene-families by abundance (coverage values)
        - threshold 1: curve needs to have a median coverage higher than 2X ("min_coverage")
        - threshold 2: left  plateau side needs to be lower  than value "left_max",
                       right plateau side needs to be higher than value "right_min"
        position_median        = 0.5  ( x genome length) of sorted gene-family vector
        position_plateau_left  = 0.30 ( x genome length)
        position_plateau_right = 0.70 ( x genome length)
    """
    VERBOSE = args.verbose
    sample_stats = defaultdict(dict)
    th_max_zero = ZERO_NON_PLATEAU_TH

    if VERBOSE:
        print(' [I] Minimum median coverage threshold: '                          + str(args.min_coverage))
        print(' [I] Left maximum plateau threshold: '                             + str(args.left_max))
        print(' [I] Right minimum plateau threshold: '                            + str(args.right_min))
        print(' [I] Maximum zero non-plateau threshold (multistrain detection): ' + str(th_max_zero))

    for sample in sorted(norm_samples_coverages.keys()):
        ordered_cov = [x[1] for x in norm_samples_coverages[sample].items()] # extract all coverage value from sample
        ordered_cov = sorted(ordered_cov, reverse = True)

        leftcov = ordered_cov[int(avg_genome_length * 0.3)]
        rightcov = ordered_cov[int(avg_genome_length * 0.7)]
        loc = int(avg_genome_length * 1.25) # sample may have less gene-families than N*1.25
        zerocov = ordered_cov[loc] if (len(ordered_cov) > loc) else 0

        if VERBOSE:
            print(' [I] ' + sample + ' median coverage: ' + str(round( median_cov[sample],2)) +
                  '; left-side cov: ' + str(round(leftcov, 2)) +
                  '; right-side cov: ' + str(round(rightcov, 2)) +
                  '; out-plateau cov: ' + str(round(zerocov, 2)) )

        sample_stats[sample] = {'strainCoverage' :  median_cov[sample]}
        sample_stats[sample].update({'accepted' : median_cov[sample] >= args.min_coverage})
        if not sample_stats[sample]['accepted']:
            print('\t' + sample + ': no strain detected, sample below MIN COVERAGE threshold')
        else: # check left right plateau coverage
            if leftcov > args.left_max:
                sample_stats[sample]['accepted'] = False
                if VERBOSE:
                    print('\t' + sample + ': no strain detected, sample does not pass LEFT-side coverage threshold.')
            elif rightcov < args.right_min:
                sample_stats[sample]['accepted'] = False
                if VERBOSE:
                    print('\t' + sample + ': no strain detected, sample does not pass RIGHT-side coverage threshold.')
        if sample_stats[sample]['accepted']:
            if VERBOSE: print('\t ' + sample + ' OK - strain detected')
            sample_stats[sample].update({'Multistrain' : zerocov > th_max_zero})
            if zerocov > th_max_zero:
                if VERBOSE: print('\t' + sample + ' WARNING: sample may contain multiple strains')
        else:
            sample_stats[sample].update({'Multistrain' : False})
    # accepted_samples_list = sorted([s for s in sample2accepted if sample2accepted[s]])
    return sample_stats

def plot_dna_coverage(samples_coverages, sample_stats, genome_length, args, normalized ):
    """Plot gene-family coverage plots.
    a) absolute coverage
    b) median normalized coverage
    Accepted samples are plotted in colors, rejected samples in gray.
    """
    sample2color = {}
    try:
        import matplotlib      # for non-interactive plots on server without X11
        matplotlib.use('Agg')  # set 'Agg' before import pylab
        import matplotlib.pyplot as plt
        try:
            from pylab import legend, savefig

            samples = sorted(samples_coverages.keys())
            accepted2samples = defaultdict(list)
            num_accepted = 0
            for s in samples:
                if sample_stats[s]['accepted']:
                    accepted2samples[s] = True
                    num_accepted = num_accepted +1
                else:
                    accepted2samples[s] = False

            if normalized : # find a way to define this kind of of boolean
                plot_name = args.o_covplot_normed
                plt.suptitle('Gene families coverages')
                plt.ylabel('Coverage')
            else:
                plot_name = args.o_covplot
                plt.suptitle('Gene families normalized coverages')
                plt.ylabel('Normalized coverage')

            if not plot_name == '':
                used_colors = []
                for sample in accepted2samples.keys():
                    color, reset = random_color(used_colors)
                    sample2color[sample] = color
                    if reset:
                        used_colors = [color]
                    else:
                        used_colors.append(color)
                
                plt.xlabel('Gene families')
                for s in samples:
                    covs = samples_coverages[s].values() # also finc a way to extract covs from here
                    covs = sorted(list(covs), reverse =True)
                    if accepted2samples[s]:
                        plt.plot(range(1, len(covs) +1), covs, sample2color[s], label=s)
                    #elif not sum(covs) == 0:
                    #    plt.plot(range(1, len(covs) +1), covs, '#c0c0c0')
                plt.axis([0.0, genome_length * 1.5, 0.0, 15])
                try:
                    if num_accepted > 0: plt.legend(loc='upper right', fontsize='xx-small')
                except TypeError:
                    print(' [W] pylab.legend fontsize does not work (please update your "pylab" module version).')
                savefig(plot_name, dpi = 300)
                plt.close()

                del(samples)
                del(accepted2samples)

        except ImportError:
            print(' [W] "pylab" module is not installed.')
            print('     To visualize and save charts, you need both "matplotlib" and "pylab" modules.')
    except ImportError:
        print(' [W] "matplotlib" module is not installed.')

# ------------------------------------------------------------------------------
#  STEP 4 Define strain-specific gene-families presence/absence
# ------------------------------------------------------------------------------
def index_of(th_non_present, th_present, th_multicopy, normalized_coverage):
    if normalized_coverage < th_non_present:
        return -3
    elif normalized_coverage <= th_present:
        return -2
    elif normalized_coverage <= th_multicopy:
        return  1
    else:
        return -1

def get_idx123_plateau_definitions(sample_stats, norm_samples_coverages, families, args):
    """-o_idx HMP_saureus_DNAindex.csv
    To use later also in RNA-seq, we need an DNA index matrix containing 4 levels (1, -1, -2, -3)

    Take samples that passed plateau criteria and define index based on coverage level of gene-families
         1 means plateau area of gene-families
        -1 means multicopy core genes (left from plateau), present also in other species
        -2 means undefined gene-families between plateau-level and zero
        -3 means "clearly" non-present gene-families
    """
    sample2family2dnaidx = defaultdict(dict)
    accepted_samples = []
    for x in sample_stats:
        if sample_stats[x]['accepted']:  accepted_samples.append(x)
    accepted_samples = sorted(accepted_samples)

    for sample in accepted_samples:
        if args.verbose: print(' [I] Get DNA 1,-1,-2,-3 levels for sample ' + sample)
        for family in families:
            sample2family2dnaidx[sample][family] = index_of(args.th_non_present, args.th_present, args.th_multicopy, norm_samples_coverages[sample][family])

    if args.o_idx and len(accepted_samples) > 0:
        with open(args.o_idx, mode='w') as OUT:
            OUT.write('\t' + '\t'.join(accepted_samples) + '\n')
            for family in families:
                OUT.write(family)
                for sample in accepted_samples:
                    OUT.write('\t' + str(sample2family2dnaidx[sample][family]))
                OUT.write('\n')
    elif len(accepted_samples) == 0:
        print('[W] No DNA 1,2,3 index file has been written because no strain was detected.')
    return sample2family2dnaidx

# ------------------------------------------------------------------------------
#  STEP 5 Get presence/absence of gene-families
# ------------------------------------------------------------------------------
def get_genefamily_presence_absence(sample2family2dnaidx, sample_stats, avg_genome_length, args):
    """Get the gene-family presence/absence matrix.
    Convert the 1,2,3 index matrix:
    gene family in sample has DNA index  1 or -1 ==> present (1)
    gene family in sample has DNA index -2 or -3 ==> NOT present (0)
    """
    sample2family2presence = defaultdict(dict)
    dna_samples = sorted(sample2family2dnaidx.keys())
    # get keys (families names) from the sub dict.
    # Since information is already here no need to pass further arg to the function
    families = sample2family2dnaidx[dna_samples[0]].keys()

    for f in families:
        for sample in dna_samples:
            sample2family2presence[sample][f] = sample2family2dnaidx[sample][f] >= -1

    # get number of gene-families per sample (add to dict sample_stats)
    for sample in sample2family2presence.keys():
        numGeneFamilies = sum( sample2family2presence[sample][f] for f in sample2family2presence[sample] )
        sample_stats[sample].update({'numberGeneFamilies' : numGeneFamilies})

    if args.verbose:
        print(' [I] Number of gene families per sample-specific strain:')
        for sample in sorted(sample_stats.keys()):
            if 'numberGeneFamilies' in sample_stats[sample]:
                print('      ' + sample + '\t' + str(sample_stats[sample]['numberGeneFamilies']))
        print('      Average number of gene-families in reference genomes: ' + str(avg_genome_length))

    if len(dna_samples) > 0:
        if args.verbose:
            print(' [I] Gene family presence/absence matrix is printed to ' + args.o_matrix)
    else:
        print('[W] No file has been written for gene-family presence/absence because no strain could be detected in any of your samples.')
        print('    (a) You can try the very sensitive options:  --min_coverage 1 --left_max 1.70 --right_min 0.30')
        print('    (b) You can check the gene-family coverage curves of your samples, using options:  --o_covplot_normed covplot_normed.png')
        print('    (c) If your reads are shorter than 70bp, you can run panphlan_map.py again, using a lower min read length: --readLength 60')

    # no need to return sample_stats for dictionnaries are passed by reference in Python
    return sample2family2presence

# ------------------------------------------------------------------------------
#  STEP 5b
# ------------------------------------------------------------------------------
def select_related_ref_genomes(genome2families, samples_panfamilies, args):
    """Select reference genomes similar to strains detected in samples
    to add them in the gene-family presence/absence matrix.
    """
    numof_strains = len(genome2families)
    rejected_strains = []
    i = 1
    for ref in genome2families.keys():
        if args.verbose:
            print('[I] [' + str(i) + '/' + str(numof_strains) + '] Analysing ref. genome ' + ref + '...')
            i += 1
        ref_gen_length = len(genome2families[ref])

        # Check that at least half of the strain families is present in families
        half = int(ref_gen_length * args.strain_similarity_perc / 100)
        numof_ss_families = 0
        for f in genome2families[ref]:
            if f in samples_panfamilies:
                numof_ss_families += 1
                if numof_ss_families >= half:
                    break # Exit from the loop to improve performances
        if numof_ss_families < half:
            rejected_strains.append(ref)
            if args.verbose:
                print('[W] Strain ' + ref + ' is rejected because only ' + str(numof_ss_families) + ' families are present in the samples.')

    for s in rejected_strains:
        del(genome2families[s])

    if args.verbose:
        print('[I] ' + str(len(rejected_strains)) + ' strain genomes filtered out. Strains are: ' + ', '.join(rejected_strains))
    selected_strains = sorted([s for s in genome2families if s not in rejected_strains])
    if args.verbose:
        print('[I] Selected strains are: ' + ', '.join(selected_strains))


def get_samples_panfamilies(families, sample2family2presence):
    """Get the sorted list of all the families present in the samples
    Can be a subset of the pangenome's set of families"""
    panfamilies = set()
    for f in families:
        for s in sample2family2presence:
            if sample2family2presence[s][f]:
                panfamilies.add(f)
                break
    return sorted(panfamilies)


def merge_samples_strains_presences(sample2family2presence, genome2families, args):
    """Compute gene families presence/absence for reference genomes and merge them with detected sample strain profiles
    1. merge: first samples columns, then strains columns (still keep all gene-families present in any strain)
    2. reject all strains which have less than 50% of its gene-families in common with the sample matrix.
        As total number of gene-families, we can use genome_length=2616 (saureus) for all strains.
        Means 50% = 1308 (saureus) gene-families of a strain have to be present in the sample set, otherwise strain is excluded.
    NB. Some gene-families can be present in samples, but not in the selected (>50%) strains.
        Some gene-families can be present in selected strains, but not in samples (if a strain is selected, we show all of it's gene-families).
    """
    dna_samples = sorted(sample2family2presence.keys())
    families = sample2family2presence[dna_samples[0]].keys()
    ref_genomes = genome2families.keys()
    # Get all present (in at least one sample) families
    samples_panfamilies = get_samples_panfamilies(families, sample2family2presence)
    select_related_ref_genomes(genome2families, samples_panfamilies, args)

    # ADAPT forone is dict of presence (dict of dict) and other list of families present
    # Merge the two dictionaries
    sample_and_strain_presences = {}
    for s in sample2family2presence:
        sample_and_strain_presences[s] = sample2family2presence[s]
    for ref in genome2families:
        sample_and_strain_presences[ref] = {f : f in genome2families[ref] for f in families}
    return sample_and_strain_presences

# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------
def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python3, please update Python')

    args = read_params()
    check_args(args)

    print('\nSTEP 1. Processing genes informations from pangenome file...')
    genes_info, families, genome2families = read_pangenome(args.pangenome)
    if args.add_ref:
        print('\nSTEP 1b. Get genes present in reference genomes...')
        ref2family2presence = build_ref2family2presence(families, genome2families, args.verbose)
        if args.i_dna == None and args.i_covmat == None:
            print('\nSTEP 1c. Print presence/absence binary matrix only for reference genomes...')
            if not args.func_annot is None:
                family2annot = create_annot_dict(ref2family2presence, args)
                write_presence_absence_matrix(ref2family2presence, args, family2annot)
            else:
                write_presence_absence_matrix(ref2family2presence, args, None)
            sys.exit(0)
        
        
    if args.i_covmat == None:
        # no shortcut        
        print('\nSTEP 2. Create coverage matrix')
        dna_samples_covs = read_map_results(args.i_dna, args.verbose)
        # Merge gene/transcript abundance into family (normalized) coverage
        for sample in sorted(dna_samples_covs.keys()):
            if args.verbose: print(' [I] Gene family normalization for DNA sample ' + sample + '...')
            dna_samples_covs[sample] = get_genefamily_coverages(dna_samples_covs[sample], genes_info, args.verbose)
            # dict of samples, for each sample : nested dict with familly and normalized coverage
        if args.o_covmat:
            print_coverage_matrix(dna_samples_covs, args.o_covmat, families, args.verbose)
    else:
        # shortcut possible, precomputed coverage matrix available
        print('\nSTEP 2. Read provided coverage matrix')
        dna_samples_covs = read_coverage_matrix(args.i_covmat)
    

    print('\nSTEP 3: Strain presence/absence filter based on coverage plateau curve...')
    avg_genome_length = adjust_genome_length(genome2families)
    norm_samples_coverages, median_cov = defining_normalized_coverage(dna_samples_covs, avg_genome_length, families)
    sample_stats = strain_presence_plateau_filter(norm_samples_coverages, avg_genome_length, median_cov, args)
    # if not args.o_covplot is None:
    #     plot_dna_coverage(dna_samples_covs, sample_stats, avg_genome_length, normalized = False, args)
    if args.o_covplot_normed:
        plot_dna_coverage(norm_samples_coverages, sample_stats, avg_genome_length, args, normalized = True)


    print('\nSTEP 4: Define strain-specific gene-families presence/absence (1,-1,-2,-3 matrix, option --o_idx)')
    sample2family2dnaidx = get_idx123_plateau_definitions(sample_stats, norm_samples_coverages, families, args)


    print('\nSTEP 5: Get presence/absence of gene-families (1,-1 matrix, option --o_matrix)')
    sample2family2presence = get_genefamily_presence_absence(sample2family2dnaidx, sample_stats, avg_genome_length, args)

    # ADD STRAINS PRESENCE ABSCENCE IF NEEDED
    if args.add_ref:
        print('\nSTEP 5b: Add reference genomes in matrix of presence/absence')
        # ss_presence = { STRAIN or SAMPLE : { GENE FAMILY : PRESENCE T or F } }
        ss_presence = merge_samples_strains_presences(sample2family2presence, genome2families, args)

    if args.func_annot:
        print('\nOPTIONAL STEP: Adding functionnal annotation of genes... (option --func_annot)')
        family2annot = create_annot_dict(sample2family2presence, args)
    else:
        family2annot = None

    if args.o_matrix:
        print('\nFINAL STEP: Writing presence/absence matrix...')
        if args.add_ref:
            write_presence_absence_matrix(ss_presence, args, family2annot)
        else:
            write_presence_absence_matrix(sample2family2presence, args, family2annot)



if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
