#!/usr/bin/env python

from __future__ import with_statement 

# ==============================================================================
# PanPhlAn v1.0: PANgenome-based PHyLogenomic ANalysis
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

__author__  = 'Thomas Tolio, Matthias Scholz, Nicola Segata (panphlan-users@googlegroups.com)'
__version__ = '1.0.2'
__date__    = '28 May 2015'

# Imports
from argparse import ArgumentParser
from collections import defaultdict
import os, subprocess, sys, tempfile, time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from fnmatch import fnmatch
import re # for gene genome mapping
    
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
INEXISTENCE_ERROR_CODE      = 1 # File or folder does not exist
UNINSTALLED_ERROR_CODE      = 2 # Software is not installed
FILEFORMAT_ERROR_CODE       = 3 # FFN file content not in NCBI format
INTERRUPTION_ERROR_CODE     = 7 # Computation has been manually halted

# ------------------------------------------------------------------------------
# INTERNAL CLASSES
# ------------------------------------------------------------------------------

class PanPhlAnGenParser(ArgumentParser):
    '''
    Subclass of ArgumentParser for parsing command inputs for panphlan_pangenome_generation.py
    '''
    def __init__(self):
        ArgumentParser.__init__(self)
        self.add_argument('--i_ffn',         metavar='INPUT_FFN_FOLDER',     type=str,   required=True,  help='Folder containing the .ffn files, i.e. the files of gene sequences of all the genomes for pangenome generation.')
        self.add_argument('--i_fna',         metavar='INPUT_FNA_FOLDER',     type=str,   required=True,  help='Folder containing the .fna files, i.e. the files of genomes sequences for Bowtie2 indexes generation.')
        self.add_argument('-c','--clade',    metavar='CLADE_NAME',           type=str,   required=True,  help='Name of the specie to consider, i.e. the basename of the index for the reference genome used by Bowtie2 to align reads.')
        self.add_argument('-o','--output',   metavar='OUTPUT_FOLDER',        type=str,   required=True,  help='Directory where to store the produced files (six .bt2 files for Bowtie2 indexes, one .csv file for the pangenome).')
        self.add_argument('--th',            metavar='IDENTITY_PERCENATGE',  type=float, default=95.0,   help='Threshold of gene sequence similarity (in percentage). Default value is 95.0 %%.')
        self.add_argument('--tmp',           metavar='TEMP_FOLDER',          type=str,   default='TMP_panphlan', help='Alternative folder for temporary files.')
        self.add_argument('--uc',            action='store_true',                                        help='Defines if to keep usearch7 output (mainly centroids.ffn and the pangenome-clusters)')
        self.add_argument('--verbose',       action='store_true',                                        help='Defines if the standard output must be verbose or not.')
        self.add_argument('-v', '--version', action='version',   version="PanPhlAn version "+__version__+"\t("+__date__+")", help='Prints the current PanPhlAn version and exits.')


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

def create_bt2_indexes(ffn_folder, fna_folder, clade, output_path, tmp_path, TIME, VERBOSE):
    '''
    Call the build function of Bowtie2 to create the indexes for a given species
    
    Merge all genomes into a single file and run bowtie2-build 
        cat /ecoli_genomes_2014/*.fna > genomes.fna
        bowtie2-build genomes.fna panphlan_ecoli14
        bowtie2-inspect -n panphlan_ecoli14
    '''
    # tmp_fna == genomes.fna
    try:

        if tmp_path == None:
            tmp_fna = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.fna')
        else:
            tmp_fna = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.fna', dir=tmp_path)

        cat_cmd = ['cat']
        genomefiles = [f for f in os.listdir(fna_folder) if fnmatch(f,'*.'+FNA)]
        for f in genomefiles:
            path_genomefile_fna = fna_folder + f
            path_genefile_ffn = ffn_folder + f.replace('.'+FNA,'.'+FFN)
            if not os.path.exists(path_genefile_ffn):
                print('[W] Cannot find gene-file:\n    ' + path_genefile_ffn)
                print('    Excluding genome ' + f + ' from bowtie2 index')
            else:
                cat_cmd.append(path_genomefile_fna)
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

def get_gene_locations(pathgenomefiles, pathgenefiles, VERBOSE):
    '''
    Get gene locations: Read all .ffn files to extract start and stop location from gene-name,
    If location cannot be extracted from gene-name, use a blast-like mapping of genes against their genomes. 
    
    Gene-locations are returned as dictionary {geneID:(contig,start,stop)}
    
    Requires: Biopython module
    '''
    gene2loc = defaultdict(tuple)
    for (genomefile, genefile) in zip(pathgenomefiles,pathgenefiles):
        try: # extract gene-location from geneIDs
            if VERBOSE:
                print('[I] ' + genefile + ': Extract gene-location from geneIDs')
            for r in SeqIO.parse(open(genefile, mode='r'), 'fasta'): 
                # extract gene-locations from gi-gene-IDs, examples
                #   gi|545636471|ref|NC_022443.1|:3480-3965
                #   gi|387779217|ref|NC_017349.1|:789327-789398,789400-790437
                pos1 = int(r.id.split(':')[1].split('-')[0].replace('c',''))
                pos2 = int(r.id.split(':')[1].split('-')[-1].replace('c',''))
                contig = r.id.split(':')[0]
                start, stop = min(pos1, pos2), max(pos1, pos2) # to always have start < stop
                gene2loc[r.id] = (str(contig), start, stop)
        except IndexError as err: # alternatively, run BLAST-like python gene-genome mapping to get locations
            if VERBOSE:
                print('    Extraction from geneID failt, map gene-sequences against genome')
            gene2multiloc = {} # tmp-dict for all hits, including sets of multiple gene locations
            # read genome sequence
            contignames = []
            contigseqs = []
            for c in SeqIO.parse(open(genomefile, mode='r'), 'fasta'):
                contignames.append(c.id)
                contigseqs.append(str(c.seq))
            # loop over gene sequences
            for g in SeqIO.parse(open(genefile, mode='r'), 'fasta'):
                gene2multiloc[g.id] = [] # append, to join hits of all contigs
                for (cn,cseq) in zip(contignames,contigseqs):
                    loc=(m.start() for m in re.finditer(str(g.seq)+'|'+str(g.seq.reverse_complement()),cseq))
                    for s in loc:  # genes can have multiple hits
                        gene2multiloc[g.id].append((cn,s+1,s+len(g))) # geneID:[(contigID, start, stop),(contigID, start, stop),...]
            # convert single hits into final dict
            delKeys=[] 
            for k, v in gene2multiloc.items():
                if len(v)==0:
                    print('gene: ' + k)
                    print('[W] gene sequence does not match genome sequence')
                    delKeys.append(k)
                    # del gene2multiloc[k] # cannot del during iteration: Python3 RuntimeError: dictionary changed size during iteration
                elif len(v)==1: 
                    gene2loc[k]=gene2multiloc[k][0]
                    delKeys.append(k)
                    # del gene2multiloc[k]
            for k in delKeys: # remove empty or single hits, keep multiple hits
                del gene2multiloc[k]
            # handle remaining multiple gene hits (multi-copy genes)
            unique_gene_loc=[] # get unique sets of multi-copy gene locations
            for k, v in gene2multiloc.items():
                if v not in unique_gene_loc:
                    unique_gene_loc.append(v)
            unique_geneIDsets=[] # get all geneIDs that hit to same location-set   
            for i in unique_gene_loc:
                geneIDset=[]
                for k, v in gene2multiloc.items():
                    if v==i:
                        geneIDset.append(k)   
                unique_geneIDsets.append(geneIDset)    
            # add multi-copy genes to gene2loc dictionary (assign to each multi-copy geneID a different location)    
            for geneIDset,locSet in zip(unique_geneIDsets,unique_gene_loc):
                for g,c in zip(geneIDset,locSet):    
                    gene2loc[g]=c
    return gene2loc


# ------------------------------------------------------------------------------

def get_contigs(fna_folder):
    '''
    Map each genome to its own set of contigs
    NB. Use Biopython
    '''
    # { CONTIG : GENOME }
    genome2contigs = defaultdict(set)
    genomefiles = [f for f in os.listdir(fna_folder) if fnmatch(f,'*.'+FNA)]
    # for root, dirs, files in os.walk(fna_folder):
    for f in genomefiles:
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
    genefiles = [f for f in os.listdir(ffn_folder) if fnmatch(f,'*.'+FFN)]
    # for root, dirs, files in os.walk(ffn_folder):
    for f in genefiles:
        # Check that the file extension is 'ffn'
        extension = os.path.splitext(f)[1].replace('.', '')
        if extension == FFN:
            genome = f[:-4]
            # print('[I] Genome ' + genome + '...\r') # shows also excluded files
            for seq_record in SeqIO.parse(open(ffn_folder + f, mode='r'), 'fasta'):
                # Take the record ID as the name of the gene, and add it to the list of genes for this genome
                gene2genome[seq_record.id] = genome
    return gene2genome


# ------------------------------------------------------------------------------

def pangenome_generation(pathgenomefiles, pathgenefiles, ffn_folder, fna_folder, merged_txt, clade, output_path, TIME, VERBOSE):
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
    gene2family    = familydictization(merged_txt, VERBOSE)
    gene2loc       = get_gene_locations(pathgenomefiles, pathgenefiles, VERBOSE)
    gene2genome    = gene2genome_mapping(ffn_folder, VERBOSE)
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
    centroids_ffn = output_path + 'panphlan_' + clade + '_centroids.ffn'
    # 3rd command:  usearch7 -cluster_smallmem merged_file.sorted.ffn -id 0.95 -maxaccepts 32 -maxrejects 128 -wordlength 3 -strand both -uc merged_file.uc
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

def merging(ffn_folder, fna_folder, tmp_path, TIME, VERBOSE):
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
            # 1st command: cat GENE_SEQ_FOLDER/*.ffn > merged_file.ffn
            cat_cmd = ['cat'] # same as for bowtie2, but based on gene.ffn files 
            genefiles = [f for f in os.listdir(ffn_folder) if fnmatch(f,'*.'+FFN)]
            for f in genefiles:
                path_genefile_ffn = ffn_folder + f
                path_genomefile_fna = fna_folder + f.replace('.'+FFN,'.'+FNA)
                if not os.path.exists(path_genomefile_fna):
                    print('[W] Cannot find genome-file:\n    ' + path_genomefile_fna)
                    print('    Excluding corresponding genes of file: ' + f + ' from usearch7 clustering')
                else:
                    cat_cmd.append(path_genefile_ffn)
            if VERBOSE:
                print('[C] ' + ' '.join(cat_cmd) + ' > ' + tmp_ffn.name)
            p1 = subprocess.Popen(cat_cmd, stdout=tmp_ffn)
            p1.wait()
            if VERBOSE:
                print('[I] FFN files have been merged into one.')

            #cat_cmd = ['cat']
            #for root, dirs, files in os.walk(ffn_folder):
            #    for f in files:
            #        # Check that the extension is '.ffn'
            #        extension = os.path.splitext(f)[1].replace('.', '')
            #        if extension == FFN:
            #            cat_cmd.append(ffn_folder + f)
            ## NB. It seems that we need to declare explicitly all the ffn files becuase Popen with the '*.ffn' does not work
            #if VERBOSE:
            #    print('[C] ' + ' '.join(cat_cmd) + ' > ' + tmp_ffn.name)
            #p1 = subprocess.Popen(cat_cmd, stdout=tmp_ffn)
            #p1.wait()
            #if VERBOSE:
            #    print('[I] FFN files have been merged into one.')

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

def gene_families_clustering(ffn_folder, fna_folder, identity_threshold_perc, clade, output_path, tmp_path, KEEP_UC, TIME, VERBOSE):
    '''
    
    NB. If KEEP_UC, then <clusters>.uc is a file written in the output directory.
        Otherwise, <clusters>.uc is a temp file (in /tmp), deleted at the end of the computation
    '''
    # Merge & Sort
    TIME, tmp_sorted_ffn = merging(ffn_folder, fna_folder, tmp_path, TIME, VERBOSE)
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
            print('[I] Biopython module is installed.')
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
    Check if Usearch 7 is already installed
    '''
    try:
        if PLATFORM == LINUX:
            usearch7_path = subprocess.Popen(['which','usearch7'], stdout=subprocess.PIPE).communicate()[0]
        elif PLATFORM == WINDOWS:
            usearch7_path = subprocess.Popen(['where','usearch7'], stdout=subprocess.PIPE).communicate()[0]
        usearch7_version = subprocess.Popen(['usearch7','--version'], stdout=subprocess.PIPE).communicate()[0]
        usearch7_version = usearch7_version.split()[1]
    except OSError as err:
        show_error_message(err)
        print('\n[E] Please, install Usearch 7.\n')
        if VERBOSE:
            print('    Usearch 7  is necessary to efficiently cluster sequences by similarity, in')
            print('    order to group genes in functional families. Gene families will be useful')
            print('    for further computation in PanPhlAn.')
        sys.exit(UNINSTALLED_ERROR_CODE)

    if VERBOSE:
        print('[I] Usearch v.7 is installed, version: ' + str(usearch7_version) + ', path: ' + str(usearch7_path).strip())
    

# ------------------------------------------------------------------------------

def check_bowtie2(VERBOSE, PLATFORM='lin'):
    '''
    Check if Bowtie2 is already installed
    '''
    try:
        if PLATFORM == LINUX:
            bowtie2 = subprocess.Popen(['which', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0]
        elif PLATFORM == WINDOWS:
            bowtie2 = subprocess.Popen(['where', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0]
        bowtie2_version = subprocess.Popen(['bowtie2', '--version'], stdout=subprocess.PIPE).communicate()[0]
        bowtie2_version = bowtie2_version.split()[2]
        if VERBOSE:
            print('[I] Bowtie2 is installed, version: ' + str(bowtie2_version) + ', path: ' + str(bowtie2).strip())
    except OSError as err:
        show_error_message(err)
        print('\n[E] Please, install Bowtie2.\n')
        if VERBOSE:
            print('    Bowtie2 is necessary to generate the specie indexes to use in further PanPhlAn')
            print('    computation. Moreover, after having generated the six index files, Bowtie2 checks')
            print('    them extracting information to show.')
        sys.exit(UNINSTALLED_ERROR_CODE)

# ------------------------------------------------------------------------------

def check_blastn(VERBOSE, PLATFORM='lin'):
    '''
    Check if BLASTn is installed
    '''
    try:
        if PLATFORM == LINUX:
            blastn_path = subprocess.Popen(['which','blastn'], stdout=subprocess.PIPE).communicate()[0]
        elif PLATFORM == WINDOWS:
            blastn_path = subprocess.Popen(['where','blastn'], stdout=subprocess.PIPE).communicate()[0]
        blastn_version = subprocess.Popen(['blastn','-version'], stdout=subprocess.PIPE).communicate()[0]
        blastn_version = blastn_version.split()[1]
    except OSError as err:
        show_error_message(err)
        print('\n[E] Please, install BLASTn.\n')
        if VERBOSE:
            print('    BLASTn is required to get gene-locations by mapping genes against genomes.')
        sys.exit(UNINSTALLED_ERROR_CODE)

    if VERBOSE:
        print('[I] BLASTn is installed, version: ' + str(blastn_version) + ', path: ' + str(blastn_path).strip())

# ------------------------------------------------------------------------------

def check_genomes(ffn_folder, fna_folder, VERBOSE):
    '''
    Check if genome files and .fna .fnn pairs are present and calculate expected runtime
    '''
    genomefiles = [f for f in os.listdir(fna_folder) if fnmatch(f,'*.'+FNA)]
    genefiles   = [f for f in os.listdir(ffn_folder) if fnmatch(f,'*.'+FFN)]

    if len(genomefiles)==0:
        print('\n[E] Cannot find any genome.fna file in folder:\n    ' + fna_folder)
        print('    Genome files need to end with .fna\n')
        sys.exit('Missing genome files')

    # check if all genome-gene file pairs exist, remove single genome or gene-files from list
    for f in genomefiles:
        path_genefile_ffn = os.path.join(ffn_folder, f.replace('.'+FNA,'.'+FFN) )
        if not os.path.exists(path_genefile_ffn):
            print('[W] Cannot find gene-file:\n    ' + path_genefile_ffn)
            print('    Excluding genome ' + f + ' from pangenome database')
            genomefiles.remove(f)
    for f in genefiles:    
        path_genomefile_fna = os.path.join(fna_folder, f.replace('.'+FFN,'.'+FNA) )
        if not os.path.exists(path_genomefile_fna):
            print('[W] Cannot find genome-file:\n    ' + path_genomefile_fna)
            print('    Excluding corresponding genes of file: ' + f + ' from pangenome database')    
            genefiles.remove(f)

    if len(genomefiles)==0:
        print('\n[E] Cannot find any pair of genome-gene files having the same filename: genome.fna gene.ffn')
        sys.exit('Missing genome-gene file pairs')

    # add full path to genefile list        
    pathgenomefiles = [os.path.join(fna_folder,f) for f in genomefiles]
    pathgenefiles   = [os.path.join(ffn_folder,f) for f in genefiles]   

    print('\nExpected runtime: ' + str(len(genomefiles)*20) + ' minutes (20 min per genome)')
    if not VERBOSE:
        print('Use option --verbose to display progress information.')

    return pathgenomefiles, pathgenefiles

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
        show_error_message('Input folder -i_ffn does not exist.')
        sys.exit(INEXISTENCE_ERROR_CODE)
        
    ipath = os.path.abspath(ipath)
    if not ipath[-1] ==  '/':
        args['i_ffn'] = ipath + '/'
    if VERBOSE:
        print('[I] Input FFN folder is ' + args['i_ffn'])

    # Check: FNA_FOLDER --------------------------------------------------------
    ipath = args['i_fna']
    if not os.path.exists(ipath):
        show_error_message('Input folder -i_fna does not exist.')
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
    # Check options correctness
    args = check_args()
    VERBOSE = args['verbose']
    KEEP_UC = args['uc']
    PLATFORM = sys.platform.lower()[0:3]
    
    pathgenomefiles, pathgenefiles = check_genomes(args['i_ffn'], args['i_fna'], VERBOSE)
    if VERBOSE:
        print('\nSTEP 0. Initialization...')
    
    TOTAL_TIME = time.time()
    TIME = time.time()

    merged_txt = ''
    
    
    
    # Check if software is installed
    if VERBOSE:
        print('\nSTEP 1. Checking required software installations...')
    # blastn = check_blastn(VERBOSE, PLATFORM) # to map genes against genomes
    bowtie2 = check_bowtie2(VERBOSE, PLATFORM) # index generation
    usearch7 = check_usearch7(VERBOSE, PLATFORM) # get gene-family cluster
    biopython = check_biopython(VERBOSE) # get geneIDs and contigIDs from ffn/fna files

    # Get gene families cluster
    if VERBOSE:
        print('\nSTEP 2. Getting gene families cluster...')
    merged_txt, TIME = gene_families_clustering(args['i_ffn'], args['i_fna'], args['th'], args['clade'], args['output'], args['tmp'], KEEP_UC, TIME, VERBOSE)
    # end else

    # Get pangenome file
    if VERBOSE:
        print('\nSTEP 3. Getting pangenome file...')
    TIME = pangenome_generation(pathgenomefiles, pathgenefiles, args['i_ffn'], args['i_fna'], merged_txt, args['clade'], args['output'], TIME, VERBOSE)
    os.remove(merged_txt) # This file is not useful anymore, and if users want to keep these information, they have the .uc file
    # Get Bowtie2 indexes
    TIME = create_bt2_indexes(args['i_ffn'], args['i_fna'], args['clade'], args['output'], args['tmp'], TIME, VERBOSE)

    end_program(time.time() - TOTAL_TIME)

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
