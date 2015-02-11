#!/usr/bin/env python

from __future__ import with_statement 

# ==============================================================================
# PanPhlAn v0.9: PANgenome-based PHyLogenomic ANalysis
#                for taxonomic classification of metagenomic data
#
# Authors: Thomas Tolio (thomas.tolio@unitn.it)
#          @TODO future contributors
#
# Please type "./panphlan_pangenome_generation.py -h" for usage help
#
# ==============================================================================

__author__  = 'Thomas Tolio (thomas.tolio@studenti.unitn.it)'
__version__ = '1.0'
__date__    = '3 February 2015'

# Imports
from argparse import ArgumentParser
from collections import defaultdict
import os, subprocess, sys, tempfile, time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
    
# Operating systems
LINUX                   = 'lin'
WINDOWS                 = 'win'
# @TODO macintosh?

# File extensions
FFN                     = 'ffn'
FNA                     = 'fna'
TXT                     = 'txt'
UC                      = 'uc'

INTERRUPTION_MESSAGE    = '[E] Execution has been manually halted.\n'

# Error codes
INEXISTENCE_ERROR_CODE      =  1 # File or folder does not exist
UNINSTALLED_ERROR_CODE      =  2 # Software is not installed
INTERRUPTION_ERROR_CODE     =  7 # Computation has been manually halted

# ------------------------------------------------------------------------------
# INTERNAL CLASSES
# ------------------------------------------------------------------------------

class PanPhlAnGenParser(ArgumentParser):
    '''
    Subclass of ArgumentParser for parsing command inputs for panphlan_pangenome_generation.py
    '''
    def __init__(self):
        ArgumentParser.__init__(self)
        self.add_argument('--i_ffn',        metavar='INPUT_FFN_FOLDER',     type=str,   required=True,  help='Folder containing the .ffn files, i.e. the files of gene sequences of all the genomes for pangenome generation.')
        self.add_argument('--i_fna',        metavar='INPUT_FNA_FOLDER',     type=str,   required=True,  help='Folder containing the .fna files, i.e. the files of genomes sequences for Bowtie2 indexes generation.')
        self.add_argument('-c','--clade',   metavar='CLADE_NAME',           type=str,   required=True,  help='Name of the specie to consider, i.e. the basename of the index for the reference genome used by Bowtie2 to align reads.')
        self.add_argument('-o','--output',  metavar='OUTPUT_FOLDER',        type=str,   required=True,  help='Directory where to store the produced files (six .bt2 files for Bowtie2 indexes, one .csv file for the pangenome).')
        self.add_argument('--th',           metavar='IDENTITY_PERCENATGE',  type=float, default=95.0,   help='Threshold of gene sequence similarity (in percentage). Default value is 95.0 %.')
        self.add_argument('--tmp',          metavar='TEMP_FOLDER',          type=str,                   help='Alternative folder for temporary files.')
        self.add_argument('--uc',           action='store_true',                                        help='Defines if to keep usearch7 output (mainly centroids.ffn and the pangenome-clusters)')
        self.add_argument('--verbose',      action='store_true',                                        help='Defines if the standard output must be verbose or not.')


# ------------------------------------------------------------------------------
# MINOR FUNCTIONS
# ------------------------------------------------------------------------------

def end_program(total_time):
    print('[TERMINATING...] ' + __file__ + ', ' + str(round(total_time / 60.0, 2)) + ' minutes.')



def show_interruption_message():
	sys.stderr.flush()
	sys.stderr.write('\r')
	sys.stderr.write(INTERRUPTION_MESSAGE)



def show_error_message(error):
    sys.stderr.write('[E] Execution has encountered an error!\n')
    sys.stderr.write('    ' + str(error) + '\n')



def time_message(start_time, message):
    current_time = time.time()
    print('[I] ' + message + ' Execution time: ' + str(round(current_time - start_time, 2)) + ' seconds.')
    return current_time

# ------------------------------------------------------------------------------
# MAJOR FUNCTIONS
# ------------------------------------------------------------------------------

def create_bt2_indexes(fna_folder, clade, output_path, tmp_path, TIME, VERBOSE):
    '''
    Call the build function of Bowtie2 to create the indexes for the given specie and check them
    
    Ex.
        cat /banche_dati/sharedCM/genomes/ecoli/*.fna > genomes.fna
        bowtie2-build genomes.fna panphlan_ecoli
        bowtie2-inspect -n panphlan_ecoli
    '''
    # tmp_fna == genomes.fna
    try:

        if tmp_path == None:
            tmp_fna = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.fna')
        else:
            tmp_fna = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.fna', dir=tmp_path)


        cat_cmd = ['cat']
        for root, dirs, files in os.walk(fna_folder):
            for f in files:
                # Check that the extension is '.fna'
                extension = os.path.splitext(f)[1].replace('.', '')
                if extension == FNA:
                    cat_cmd.append(fna_folder + f)
        if VERBOSE:
            print('[C] ' + ' '.join(cat_cmd) + ' > ' + tmp_fna.name)
        p4 = subprocess.Popen(cat_cmd, stdout=tmp_fna)
        p4.wait()

        try:
            clade = output_path + 'panphlan_' + clade
            build_cmd = ['bowtie2-build', tmp_fna.name, clade]
            if VERBOSE:
                print('[C] ' + ' '.join(build_cmd))
            p5 = subprocess.Popen(build_cmd)
            p5.wait()

            try:
                # Check generated files
                inspect_cmd = ['bowtie2-inspect', '-n', clade]
                if not VERBOSE:
                    inspect_cmd.append('--verbose')
                else:
                    print('[I] ' + ' '.join(inspect_cmd))
                p6 = subprocess.Popen(inspect_cmd)
                p6.wait()

                os.unlink(tmp_fna.name)
                if VERBOSE:
                    TIME = time_message(TIME, 'Bowtie2 indexes have been created.')
                return TIME

            except (KeyboardInterrupt, SystemExit):
                p6.kill()
                show_interruption_message()
                sys.exit(INTERRUPTION_ERROR_CODE)

        except (KeyboardInterrupt, SystemExit):
            p5.kill()
            show_interruption_message()
            sys.exit(INTERRUPTION_ERROR_CODE)

    except (KeyboardInterrupt, SystemExit):
        p4.kill()
        os.unlink(tmp_fna.name)
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)


# ------------------------------------------------------------------------------

def combining(gene2loc, gene2family, gene2genome, output_path, clade, TIME, VERBOSE):
    '''
    Create the pangenome combining all the information from gene mappings (location (contig, from, to), family and genome)

        (1) Include gene families
            currently using grep to find geneID in gene-cluster file (line-number becomes genefamilyID)
            grep -nw geneID family-cluster-file.txt | 'awk -F ":" ''{print $1}'
            but very slow, needs optimization!!

        (2) Include genome names
            using contig-info from previous step (get_contigs.py to get contig names for each genome)
        
    Result: final pangenome file (tab-separated), with this format:
        geneFamily | geneID | genomeName(filename) | contigID | start | stop
    '''
    # line := FAMILy GENE GENOME CONTIG FROM TO
    pangenome_csv = output_path + 'panphlan_' + clade + '_pangenome.csv'
    with open(pangenome_csv, mode='w') as ocsv:
        genes_list = sorted(gene2loc.keys())
        for gene in genes_list:
            ocsv.write(gene2family[gene] + '\t' + gene + '\t' + gene2genome[gene] + '\t' + gene2loc[gene][0] + '\t' + str(gene2loc[gene][1]) + '\t' + str(gene2loc[gene][2]) + '\n')
    return TIME
    

# ------------------------------------------------------------------------------

def get_gene_locations(ffn_folder):
    '''
    Read the .ffn files in order to generate a dictionary that map each gene to its location in its genome
    Produce a .txt file with the locations and return a dictionary {gene:(contig,from,to)}

    NB. Need of Biopython module
    '''
    # { GENE : ( CONTIG, FROM, TO ) }
    gene2loc = defaultdict(tuple)
    for root, dirs, files in os.walk(ffn_folder):
        for f in files:
            for r in SeqIO.parse(open(ffn_folder + f, mode='r'), 'fasta'): # A .ffn file is a FASTA file
                # Through SeqIO, automatially parse the .ffn file and extract the information
                pos1 = int(r.id.split(':')[1].split('-')[0].replace('c',''))
                pos2 = int(r.id.split(':')[1].split('-')[-1].replace('c',''))
                contig = r.id.split(':')[0]
                # Automatically assign to 'from' the min value, and to 'to' the max, even if panphlan_map check the order and switch them if inverted
                start, stop = min(pos1, pos2), max(pos1, pos2)
                gene2loc[r.id] = (str(contig), start, stop)
    return gene2loc


# ------------------------------------------------------------------------------

def get_contigs(fna_folder):
    '''
    Map each genome to its own set of contigs
    NB. Use Biopython
    '''
    # { CONTIG : GENOME }
    genome2contigs = defaultdict(set)
    for root, dirs, files in os.walk(fna_folder):
        for f in files:
            # Check that the file extension is 'fna'
            extension = os.path.splitext(f)[1].replace('.', '')
            if extension == FNA:
                genome = f.split('.')[0].split('/')[-1]
                ltot = 0
                s = ''
                for r in SeqIO.parse(open(fna_folder + f, mode='r'), 'fasta'):
                    l = len(r.seq)
                    ltot += l
                    s = s + "\t" + r.id
                    genome2contigs[genome].add(r.id)
    return genome2contigs


# ------------------------------------------------------------------------------

def gene2genome_mapping(ffn_folder, VERBOSE):
    '''
    Map each gene to its own genome
    NB. Use Biopython
    '''
    gene2genome = {}

    for root, dirs, files in os.walk(ffn_folder):
        for f in files:
            # Check that the file extension is 'ffn'
            extension = os.path.splitext(f)[1].replace('.', '')
            if extension == FFN:
                genome = f[:-4]
                print('[I] Genome ' + genome + '...\r')
                for seq_record in SeqIO.parse(open(ffn_folder + f, mode='r'), 'fasta'):
                    # Take the record ID as the name of the gene, and add it to the list of genes for this genome
                    gene2genome[seq_record.id] = genome

    return gene2genome


# ------------------------------------------------------------------------------

def pangenome_generation(ffn_folder, fna_folder, merged_txt, clade, output_path, TIME, VERBOSE):
    '''
    TODO
    
        (1) Extract gene locations from gene-identifier in .ffn files
            currently done by this script: get_gene_locations.py FOLDER/*.ffn > gene_locations.txt 
            NB. A gene can have multiple regions: e.g. "789327-789398,789400-790437"
                In this case we simply take the complete region (ignoring little gaps):
                start:= 789327, stop:= 790437
            NB. Locations can already be sorted such that: start < stop (even though panphlan_map is checking again)

        (2) Get contig names for each genome-name (filename)
            currently done by: get_contigs.py FOLDER/*.fna > panphlan_CLADE_contigs.txt
            (genome-name/file-name is in first column)
            NB. We need this to know which contigs belong to which filename (genomename)

    Result: basic part of the final pangenome-file (tab-separated):
        geneID | start | stop
    '''

    if VERBOSE:
        print('[I] Get gene locations, gene families, contigs and genomes for each gene.')
    gene2family = familydictization(merged_txt, VERBOSE)
    gene2loc = get_gene_locations(ffn_folder)
    gene2genome = gene2genome_mapping(ffn_folder, VERBOSE)
    genome2contigs = get_contigs(fna_folder)
    
    # Create the pangenome
    combining(gene2loc, gene2family, gene2genome, output_path, clade, TIME, VERBOSE)

    if VERBOSE:
        TIME = time_message(TIME, 'Pangenome has been generated.')

    return TIME


# ------------------------------------------------------------------------------

def family_of(index):
    return 'g' + str(format(index, '06d'))

def familydictization(merged_txt, VERBOSE):
    '''
    Return a dictionary mapping genes to their own gene family
    gene2family := { GENE : FAMILY }

    Input: a file with genes are clustered by line

    Gene family is named as 'g<NUM OF LINE IN THE FILE>'
    EX.
        line 1    ==> family 'g000001'
        line 317  ==> family 'g000317'
    '''
    gene2family = {}
    numof_line = 0 
    with open(merged_txt, mode='r') as itxt:
        for line in itxt:
            numof_line += 1
            family = family_of(numof_line)
            if '\t' in line:
                words = line.strip().split('\t')
                for gene in words:
                    gene2family[gene] = family
            else:
                gene2family[line.strip()] = family
    if VERBOSE:
        print('[I] Pangenome contains ' + str(len(gene2family)) + ' genes, clustered in ' + str(numof_line) + ' gene families.')
    return gene2family


# ------------------------------------------------------------------------------

def conversion(merged_uc, merged_txt, TIME, VERBOSE):
    '''
    Convert the UC file into a TXT file
    See also: http://drive5.com/usearch/manual/ucout.html
    '''
    # Dictionarization of the UC file's content
    uc2cl = defaultdict(set)
    with open(merged_uc, mode='r') as iuc:
        for typ, cln, seql, pid, strand, ing1, ign2, aln, query, target in (line.split('\t') for line in iuc):
            if typ == 'H':
                uc2cl[target.strip()].add( query )
            elif typ == 'S' and  query not in uc2cl:
                uc2cl[query] = set()
    
    # Printing in the TXT file
    with open(merged_txt, mode='w') as otxt:
        # k is the centroid gene, v is a gene of the its cluster
        for index, (k,v) in enumerate(sorted(uc2cl.items(), key=lambda x:str(x[0])) ,1): # clusters are sorted by centroid-IDs; added index starting at 1    
            # Each line of the clusters .txt file is composed by: GENE_FAMILY | CENTROID_GENE | GENE | ... | GENE
            otxt.write(family_of(index) + '\t' + '\t'.join([k] + sorted(list(v))) + '\n') # Intra-line sorting
    
    if VERBOSE:
        TIME = time_message(TIME, 'UC --> TXT conversion has been done.')
    return TIME


# ------------------------------------------------------------------------------

def clustering(sorted_merged_ffn, identity, clade, output_path, tmp_path, KEEP_UC, TIME, VERBOSE):
    '''
    Group gene sequence in clusters by similarity
    Default similarity is 95%
    '''
    merged_uc_name = output_path + 'usearch7_' + clade + '_cluster.uc'
    if KEEP_UC:
        merged_uc = open(merged_uc_name, mode='w')
    else:
        if tmp_path == None:
            merged_uc = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_usearch7_', suffix='.uc')
        else:
            merged_uc = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_usearch7_', suffix='.uc', dir=tmp_path)
    centroids_ffn = output_path + 'usearch7_' + clade + '_centroids.ffn'
    # 3rd command: usearch7 -cluster_smallmem merged_file.sorted.ffn -id 0.95 -maxaccepts 32 -maxrejects 128 -wordlength 3 -strand both -uc merged_file.uc
    try:
        clust_cmd = ['usearch7', '--cluster_smallmem', sorted_merged_ffn, '--id', str(identity),
                    '--maxaccepts', '32', '--maxrejects', '128', '--wordlength', '3', '--strand', 'both',
                    '--uc', merged_uc.name, '--centroids', centroids_ffn]
        if not VERBOSE:
            clust_cmd.append('--quiet')
        else:
            print('[I] ' + ' '.join(clust_cmd))
        p3 = subprocess.Popen(clust_cmd)
        p3.wait()
        if VERBOSE:
            print('[I] Clustering has been done.')
        if KEEP_UC:
            merged_uc.close()

    except (KeyboardInterrupt, SystemExit):
        p3.kill()
        # merged_uc.close()
        if KEEP_UC:
            os.remove(merged_uc_name)
        else:
            os.unlink(merged_uc.name)
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)

    finally:
        os.unlink(sorted_merged_ffn)
    
    if VERBOSE:
        TIME = time_message(TIME, 'Clustering with Usearch has been done.')
    return merged_uc, TIME

# ------------------------------------------------------------------------------

def merging(ffn_folder, tmp_path, TIME, VERBOSE):
    '''
    Merge all the FFN files into a unique one, and then sort it by length
    '''
    # NB. Here we are sure that ffn_folder finishes with a '/'
    # Create a temporary file for merging FFNs
    if tmp_path == None:
        tmp_ffn = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.ffn')
        tmp_sorted_ffn = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sorted.ffn')
    else:
        tmp_ffn = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.ffn', dir=tmp_path)
        tmp_sorted_ffn = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sorted.ffn', dir=tmp_path)
    try:
        with tmp_ffn:
            # 1st command: cat INPUT_FOLDER/*.ffn > merged_file.ffn
            cat_cmd = ['cat']
            for root, dirs, files in os.walk(ffn_folder):
                for f in files:
                    # Check that the extension is '.ffn'
                    extension = os.path.splitext(f)[1].replace('.', '')
                    if extension == FFN:
                        cat_cmd.append(ffn_folder + f)
            # NB. It seems that we need to declare explicitly all the ffn files becuase Popen with the '*.ffn' does not work
            if VERBOSE:
                print('[C] ' + ' '.join(cat_cmd) + ' > ' + tmp_ffn.name)
            p1 = subprocess.Popen(cat_cmd, stdout=tmp_ffn)
            p1.wait()
            if VERBOSE:
                print('[I] FFN files have been merged into one.')

        try:
            with tmp_sorted_ffn:
                # 2nd command: usearch7 -sortbylength merged_file.ffn -output merged_file.sorted.ffn -minseqlength 1
                sort_cmd = ['usearch7', '--sortbylength', tmp_ffn.name, '--output', tmp_sorted_ffn.name, '--minseqlength', '1']
                if not VERBOSE:
                    sort_cmd.append('--quiet')
                else:
                    print('[I] ' + ' '.join(sort_cmd))
                p2 = subprocess.Popen(sort_cmd)
                p2.wait()
                if VERBOSE:
                    print('[I] Merged FFN has been sorted.')

        except (KeyboardInterrupt, SystemExit):
            p2.kill()
            os.unlink(tmp_sorted_ffn.name)
            show_interruption_message()
            sys.exit(INTERRUPTION_ERROR_CODE)

    except (KeyboardInterrupt, SystemExit):
        p1.kill()
        show_interruption_message()
        sys.exit(INTERRUPTION_ERROR_CODE)
    finally:
        os.unlink(tmp_ffn.name)

    if VERBOSE:
        TIME = time_message(TIME, 'FFN merging and Usearch sorting has been done.')
    # Get in output the temporary merged and sorted .ffn file
    return TIME, tmp_sorted_ffn


# ------------------------------------------------------------------------------

def gene_families_clustering(ffn_folder, identity_threshold_perc, clade, output_path, tmp_path, KEEP_UC, TIME, VERBOSE):
    '''
    
    NB. If KEEP_UC, then <clusters>.uc is a file written in the output directory.
        Otherwise, <clusters>.uc is a temp file (in /tmp), deleted at the end of the computation
    '''
    # Merge & Sort
    TIME, tmp_sorted_ffn = merging(ffn_folder, tmp_path, TIME, VERBOSE)
    # Cluster
    tmp_uc, TIME = clustering(tmp_sorted_ffn.name, identity_threshold_perc / 100.0, clade, output_path, tmp_path, KEEP_UC, TIME, VERBOSE)
    # Convert
    merged_txt = output_path + 'usearch7_' + clade + '_genefamily_cluster.txt'
    TIME = conversion(tmp_uc.name, merged_txt, TIME, VERBOSE)
    if not KEEP_UC:
        os.unlink(tmp_uc.name)
    return merged_txt, TIME
    

# ------------------------------------------------------------------------------

def check_biopython(VERBOSE):
    '''
    Check if the Bio module (Biopython) is installed
    '''
    try:
        output = __import__('Bio')
        if VERBOSE:
            print('[I] Biopython module is already installedin the system.')
            return output
    except ImportError as err:
        show_error_message(err)
        print('\n[E] Please, install Biopython.\n')
        if VERBOSE:
            print('    The "Bio" module is necessary to efficiently build the mapping between a')
            print('    gene and its location in the genome (location is a tuple with contig,')
            print('    staring position and ending position).')
        sys.exit(UNINSTALLED_ERROR_CODE)


# ------------------------------------------------------------------------------

def check_usearch7(VERBOSE, PLATFORM='lin'):
    '''
    Check if Usearch 7 is already installed in the system
    '''
    try:
        if PLATFORM == LINUX:
            output = subprocess.Popen(['which', 'usearch7'], stdout=subprocess.PIPE).communicate()[0]
        elif PLATFORM == WINDOWS:
            output = subprocess.Popen(['where', 'usearch7'], stdout=subprocess.PIPE).communicate()[0]
    
    except Exception as err:
        show_error_message(err)
        print('\n[E] Please, install Usearch 7.\n')
        if VERBOSE:
            print('    Usearch 7  is necessary to efficiently cluster sequences by similarity, in')
            print('    order to group genes in functional families. Gene families will be useful')
            print('    for further computation in PanPhlAn.')
        sys.exit(UNINSTALLED_ERROR_CODE)

    if VERBOSE:
        print('[I] Usearch v.7 is already installed in the system.')


# ------------------------------------------------------------------------------

def check_bowtie2(VERBOSE, PLATFORM='lin'):
    '''
    Check if Bowtie2 is already installed in the system
    '''
    try:
        if PLATFORM == LINUX:
            bowtie2 = subprocess.Popen(['which', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0]
        elif PLATFORM == WINDOWS:
            bowtie2 = subprocess.Popen(['where', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0]
        bowtie2_version = subprocess.Popen(['bowtie2', '--version'], stdout=subprocess.PIPE).communicate()[0]
        bowtie2_version = bowtie2_version.split()[2]
        if VERBOSE:
            print('[I] Bowtie2 is already installed in the system with version ' + str(bowtie2_version) + ' in path ' + str(bowtie2).strip())
    
    except Exception as err:
        show_error_message(err)
        print('\n[E] Please, install Bowtie2.\n')
        if VERBOSE:
            print('    Bowtie2 is necessary to generate the specie indexes to use in further PanPhlAn')
            print('    computation. Moreover, after having generated the six index files, Bowtie2 checks')
            print('    them extracting information to show.')
        sys.exit(UNINSTALLED_ERROR_CODE)

# ------------------------------------------------------------------------------

def check_args():
    '''
    Check if the input arguments respect the rules of usage
    '''
    parser = PanPhlAnGenParser()
    args = vars(parser.parse_args())
    VERBOSE = args['verbose']

    # Check: FFN_FOLDER --------------------------------------------------------
    ipath = args['i_ffn']
    if not os.path.exists(ipath):
        show_error_message(err)
        sys.exit(INEXISTENCE_ERROR_CODE)
        
    ipath = os.path.abspath(ipath)
    if not ipath[-1] ==  '/':
        args['i_ffn'] = ipath + '/'
    if VERBOSE:
        print('[I] Input FFN folder is ' + args['i_ffn'])

    # Check: FNA_FOLDER --------------------------------------------------------
    ipath = args['i_fna']
    if not os.path.exists(ipath):
        show_error_message(err)
        sys.exit(INEXISTENCE_ERROR_CODE)

    ipath = os.path.abspath(ipath)
    if not ipath[-1] ==  '/':
        args['i_fna'] = ipath + '/'
    if VERBOSE:
        print('[I] Input FNA folder is ' + args['i_fna'])

    # Check: CLADE -------------------------------------------------------------
    if VERBOSE:
        print('[I] Clade is ' + args['clade'])

    # Check: IDENTITY_PERCENATGE -----------------------------------------------
    identity_threshold_perc = args['th']
    if identity_threshold_perc < 0.0 or identity_threshold_perc > 100.0:
        args['th'] = 95.0
        if VERBOSE:
            print('[I] Invalid value for identity threshold percentage. Default value (95.0 %) has been set.')
    else:
        if VERBOSE:
            print('[I] Identity threshold percentage is ' + str(args['th']) + ' %.')

    # Check: OUTPUT_FOLDER -----------------------------------------------------
    opath = args['output']
    if not opath[-1] == '/':
        opath +=  '/'
    if not os.path.exists(os.path.dirname(opath)):
        os.makedirs(opath)
    args['output'] = opath
    if VERBOSE:
        print('[I] Output file: ' + args['output'])

    # Check: TEMP_FOLDER ------------------------------------------------------
    tmp_path = args['tmp']
    if not tmp_path[-1] == '/':
        tmp_path +=  '/'
    if not os.path.exists(os.path.dirname(tmp_path)):
        os.makedirs(tmp_path)
    args['tmp'] = tmp_path
    if VERBOSE:
        print('[I] Temporary folder: ' + args['tmp'])

    return args


# ------------------------------------------------------------------------------

def main():
    print('\nSTEP 0. Initialization...')
    TOTAL_TIME = time.time()

    # Check options correctness
    args = check_args()
    VERBOSE = args['verbose']
    KEEP_UC = args['uc']
    PLATFORM = sys.platform.lower()[0:3]
    TIME = time.time()

    merged_txt = ''
    
    # Check if software is installed
    if VERBOSE:
        print('\nSTEP 1. Checking software...')
    bowtie2 = check_bowtie2(VERBOSE, PLATFORM)
    usearch7 = check_usearch7(VERBOSE, PLATFORM)
    biopython = check_biopython(VERBOSE)

    # Get gene families cluster
    if VERBOSE:
        print('\nSTEP 2. Getting gene families cluster...')
    merged_txt, TIME = gene_families_clustering(args['i_ffn'], args['th'], args['clade'], args['output'], args['tmp'], KEEP_UC, TIME, VERBOSE)
    # end else

    # Get pangenome file
    if VERBOSE:
        print('\nSTEP 3. Getting pangenome file...')
    TIME = pangenome_generation(args['i_ffn'], args['i_fna'], merged_txt, args['clade'], args['output'], TIME, VERBOSE)
    os.remove(merged_txt) # This file is not useful anymore, and if users want to keep these information, they have the .uc file
    # Get Bowtie2 indexes
    TIME = create_bt2_indexes(args['i_fna'], args['clade'], args['output'], args['tmp'], TIME, VERBOSE)

    end_program(time.time() - TOTAL_TIME)

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
