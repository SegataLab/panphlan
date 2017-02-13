#!/usr/bin/env python

# ==============================================================================
# PanPhlAn v1.2.2: PANgenome-based PHyLogenomic ANalysis
#                for detecting and characterizing strains in metagenomic samples
#
# Authors:  Matthias Scholz, algorithm design
#           Thomas Tolio, programmer
#           Nicola Segata, principal investigator
#
# PanPhlAn is a project of the Computational Metagenomics Lab at CIBIO,
# University of Trento, Italy
#
# For help type: ./panphlan_map.py -h
#
# https://bitbucket.org/CibioCM/panphlan
# ==============================================================================

from __future__ import with_statement 
from argparse import ArgumentParser
from collections import defaultdict
import os, subprocess, sys, tempfile, time
from fnmatch import fnmatch
import re # for gene genome mapping

__author__  = 'Matthias Scholz, Thomas Tolio, Nicola Segata (panphlan-users@googlegroups.com)'
__version__ = '1.2.2'
__date__    = '30 January 2017'

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError as err:
    print('\n[E]',err) 
    print('\n[E] Please install Biopython.')
    print('    The "Bio" module is required for extracting gene locations')
    print('    by mapping genes against their genome.\n')
    sys.exit(2)
    
# Operating systems
LINUX                   = 'lin'
WINDOWS                 = 'win'

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
NONUNIQUEGENE_ERROR_CODE    = 4 # More than one gene of identical geneID over the complete genome set

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
        self.add_argument('--tmp',           metavar='TEMP_FOLDER',          type=str,   default='TMP_panphlan_db', help='Alternative folder for temporary files.')
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
    sys.stderr.write('\n[E] Execution has encountered an error!\n')
    sys.stderr.write('    ' + str(error) + '\n')


def time_message(start_time, message):
    current_time = time.time()
    print('[I] ' + message + ' Execution time: ' + str(round(current_time - start_time, 2)) + ' seconds.')
    return current_time

# ------------------------------------------------------------------------------
# MAJOR FUNCTIONS
# ------------------------------------------------------------------------------
def create_bt2_indexes(pathgenomefiles, clade, output_path, tmp_path, TIME, VERBOSE):
    '''
    Call the build function of Bowtie2 to create the indexes for a given species
    
    Merge all genomes into a single file and run bowtie2-build 
        cat /ecoli_genomes_2014/*.fna > genomes.fna
        bowtie2-build genomes.fna panphlan_ecoli14
        bowtie2-inspect -n panphlan_ecoli14
    '''
    try:
        tmp_fna = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.fna', dir=tmp_path)

        cat_cmd = ['cat']
        for f in pathgenomefiles:
            cat_cmd.append(f)
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
def write_pangenome(gene2loc, gene2family, gene2genome, output_path, clade, TIME, VERBOSE):
    '''
    Create the pangenome database file combining all the information from
    gene mappings (location (contig, from, to), family and genome)
        
    Result: pangenome file (tab-separated):
        geneFamily | geneID | genomeName(filename) | contigID | start | stop
    '''
    pangenome_csv = output_path + 'panphlan_' + clade + '_pangenome.csv'
    with open(pangenome_csv, mode='w') as ocsv:
        genes_list = sorted(gene2loc.keys())
        for gene in genes_list:
            if gene in gene2family:
                ocsv.write(gene2family[gene] + '\t' + gene + '\t' + gene2genome[gene] + '\t' + gene2loc[gene][0] + '\t' + str(gene2loc[gene][1]) + '\t' + str(gene2loc[gene][2]) + '\n')
            else:
                print('[W] Could not find gene in usearch7 cluster result (dict gene2family): '   + gene)
                print('    Check presence of gene in file: usearch7_species_cluster.uc, option --uc')
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
        if VERBOSE:
            print('[I] genomefile: ' + genomefile)
            print('    genefile: '   + genefile)
        try: # extract gene-location from geneIDs
            if VERBOSE:
                print('    Extract gene-location from geneIDs')
            for r in SeqIO.parse(open(genefile, mode='r'), 'fasta'): 
                # extract gene-locations from gi-gene-IDs, examples
                #   gi|545636471|ref|NC_022443.1|:3480-3965
                #   gi|387779217|ref|NC_017349.1|:789327-789398,789400-790437
                # pos1 = int(r.id.split(':')[1].split('-')[0].replace('c',''))
                # pos2 = int(r.id.split(':')[1].split('-')[-1].replace('c',''))
                # contig = r.id.split(':')[0]
                # new: ffn files with prefix filename:origGeneID
                #   [0] = filename
                #   [1] = contig
                #   [2] = location
                pos1 = int(r.id.split(':')[-1].split('-')[0].replace('c','')) 
                pos2 = int(r.id.split(':')[-1].split('-')[-1].replace('c',''))
                contig = r.id.split(':')[-2]
                start, stop = min(pos1, pos2), max(pos1, pos2) # to always have start < stop
                gene2loc[r.id] = (str(contig), start, stop)
        except (IndexError, ValueError) as err: # alternatively, run BLAST-like python gene-genome mapping to get locations
            if VERBOSE:
                print('    Extraction from geneID failed, map gene-sequences against genome...')
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
                    loc=(m.start() for m in re.finditer(str(g.seq)+'|'+str(g.seq.reverse_complement()), cseq, re.IGNORECASE))
                    for s in loc:  # genes can have multiple hits
                        gene2multiloc[g.id].append((cn,s+1,s+len(g))) # geneID:[(contigID, start, stop),(contigID, start, stop),...]
            # convert single hits into final dict
            delKeys=[] 
            for k, v in gene2multiloc.items():
                if len(v)==0:
                    print('[W] gene sequence does not match genome sequence')
                    print('    geneID: ' + k)
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
                for g,c in zip(sorted(geneIDset),sorted(locSet)): # sorted: to get indentical results in python 2 and 3   
                    gene2loc[g]=c
    return gene2loc

# ------------------------------------------------------------------------------
def get_contigs(pathgenomefiles):
    '''
    Map each genome (filename) to its contig-set
    
    Function is not used!!
        Contig-name is already included in gene2loc (contig,start,stop)
        Genome-name (filename) is included in gene2genome
    
    requires Biopython (Bio module)
    '''
    # { genome : contig-set }
    genome2contigs = defaultdict(set)
    print('[I] Get genome-contigset list (dict)')
    for f in pathgenomefiles:
        # genome = os.path.basename(f).split('.')[0] # get genome filename without extension
        genome = os.path.splitext(os.path.basename(f))[0] # get genome filename without extension (allow dots in genome-name)
        print('    Genome: ' + genome)
        for r in SeqIO.parse(open(f, mode='r'), 'fasta'):
            genome2contigs[genome].add(r.id)
    return genome2contigs

# ------------------------------------------------------------------------------
def gene2genome_mapping(pathgenefiles, VERBOSE):
    '''
    Map each gene to its genome (genome-filename)
    
    requires Biopython (Bio module)
    '''
    gene2genome     = {} # {geneID : genomefilename}
    gene2description = {} # {geneID : description}
    print('[I] Get gene-genome list (dict)')
    for f in pathgenefiles:
        # genome = os.path.basename(f).split('.')[0] # get genome filename without extension
        genome = os.path.splitext(os.path.basename(f))[0] # get genome filename without extension (allow dots in genome-name)
        print('    Genome: ' + genome)
        for seq_record in SeqIO.parse(open(f, mode='r'), 'fasta'):
            # if not in dict, else error
            if seq_record.id not in gene2genome:    # add only if not in dictionary already
                gene2genome[seq_record.id]     = genome
                description = ' '.join(seq_record.description.split()[1:]) # remove .id from .description
                # print(description)
                gene2description[seq_record.id] = description
            else:
                print('\n[E] Error: non-unique gene identifier.')
                print('    following gene-ID appears multiple times in the genome set')
                print('    ' + seq_record.id + '\n')
                sys.exit(NONUNIQUEGENE_ERROR_CODE)
    return gene2genome, gene2description

# ------------------------------------------------------------------------------
def pangenome_generation(pathgenomefiles, pathgenefiles, gene2family, clade, output_path, gene2genome, TIME, VERBOSE):
    '''
    (1) Extract gene locations from gene-identifier or blast-like search
    (2) Get contig names for each genome-name (filename)
        We need this to know which contigs belong to which filename (genome-name)
    Result: basic part of the final pangenome-file (tab-separated):
        geneID | start | stop
    '''
    if VERBOSE: print('[I] Get gene locations and contigs for each gene.')
    gene2loc       = get_gene_locations(pathgenomefiles, pathgenefiles, VERBOSE)
    # genome2contigs = get_contigs(pathgenomefiles) # not used
    
    # Write the pangenome database file: panphlan_clade_pangenome.csv
    write_pangenome(gene2loc, gene2family, gene2genome, output_path, clade, TIME, VERBOSE)

    if VERBOSE:
        TIME = time_message(TIME, 'Pangenome has been generated.')

    return TIME

# ------------------------------------------------------------------------------
# --- usearch7 -----------------------------------------------------------------
def family_of(index):
    return 'g' + str(format(index, '06d'))
def usearch_get_gene2family_dict(merged_txt, VERBOSE):
    '''
    Return a dictionary mapping genes to gene-families, based on usearch7 cluster result
    gene2family := { GENE : FAMILY }

    Input: Text file with gene-families line-wise
           All geneIDs per line belong to one gene-family cluster
           Line 1 represents gene-family g000001

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
def convert_usearch_result(merged_uc, merged_txt, TIME, VERBOSE):
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
def usearch_centroids_add_geneID_prefix(clade, gene2family, output_path, gene2description):
    '''
    Add prefix 'clade:genefamID:' to geneIDs in centroid.ffn sequence file
    1) copy usearch7 result: panphlan_species_centroids.ffn as .._centroids_orig.ffn
    2) read file and add prefix: species:g12345:old_geneID
    3) write new version of panphlan_species_centroids.ffn
    called from --> pangenome_generation()
    '''
    centroids_ffn      = os.path.join(output_path,'panphlan_' + clade + '_centroids.ffn')
    centroids_orig_ffn = os.path.join(output_path,'panphlan_' + clade + '_centroids_orig.ffn')

    if os.path.exists(centroids_ffn):
        # move panphlan_species_centroids.ffn to panphlan_species_centroids_orig.ffn
        os.rename(centroids_ffn,centroids_orig_ffn)
        # add prefix species:g12345:old_geneID (read centroid_orig.ffn, write new centroid.ffn)    
        centroid_sequences = SeqIO.parse(open(centroids_orig_ffn),'fasta')
        with open(centroids_ffn, 'w') as f:
            for seq in centroid_sequences:
                # seq.description=''
                seq.description=gene2description[seq.id] # add description (gene annotation)
                genefamID=gene2family[seq.id] # 'g12345'
                seq.id = clade + ':' + genefamID + ':' + seq.id
                seq.name=''
                r = SeqIO.write(seq, f, 'fasta')
                if r!=1:
                    sys.exit('[E] Error while writing centroid sequence:  ' + seq.id)
        os.remove(centroids_orig_ffn)
# ------------------------------------------------------------------------------
def run_usearch(sorted_merged_ffn, identity, clade, output_path, tmp_path, KEEP_UC, TIME, VERBOSE):
    '''
    Group gene sequence in clusters by similarity
    Default similarity threshold is 95%
    '''
    merged_uc_name = output_path + 'usearch7_' + clade + '_cluster.uc'
    if KEEP_UC:
        merged_uc = open(merged_uc_name, mode='w')
    else:
        if tmp_path == None:
            merged_uc = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_usearch7_', suffix='.uc')
        else:
            merged_uc = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_usearch7_', suffix='.uc', dir=tmp_path)
    centroids_ffn = os.path.join(output_path,'panphlan_' + clade + '_centroids.ffn')
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
def usearch_sortbylength(pathgenefiles, tmp_path, TIME, VERBOSE):
    '''
    Merge all the gene-sequence FFN files into a unique one, then sort by length
    '''
    
    tmp_ffn        = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.ffn',        dir=tmp_path)
    tmp_sorted_ffn = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sorted.ffn', dir=tmp_path)
        
    try:
        with tmp_ffn:
            # cat genefiles.ffn > merged_file.ffn
            cat_cmd = ['cat'] # same as for bowtie2, but based on gene.ffn files 
            for f in pathgenefiles:
                cat_cmd.append(f)
            if VERBOSE:
                print('[C] ' + ' '.join(cat_cmd) + ' > ' + tmp_ffn.name)
            p1 = subprocess.Popen(cat_cmd, stdout=tmp_ffn)
            p1.wait()
            if VERBOSE:
                print('[I] FFN files have been merged into one.')

        try:
            with tmp_sorted_ffn:
                # usearch7 -sortbylength merged_file.ffn -output merged_file.sorted.ffn -minseqlength 1
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
def usearch_clustering(pathgenefiles, identity_threshold_perc, clade, output_path, tmp_path, KEEP_UC, gene2description, TIME, VERBOSE):
    '''
    Note: If KEEP_UC, then <clusters>.uc is a file written in the output directory.
    Otherwise, <clusters>.uc is a temp file (in /tmp), deleted at the end of the computation
    '''
    # Merge and sort genes by length
    TIME, tmp_sorted_ffn = usearch_sortbylength(pathgenefiles, tmp_path, TIME, VERBOSE)
    # usearch7 clustering
    tmp_uc, TIME = run_usearch(tmp_sorted_ffn.name, identity_threshold_perc / 100.0, clade, output_path, tmp_path, KEEP_UC, TIME, VERBOSE)
    # Convert usearch7 result
    merged_txt = output_path + 'usearch7_' + clade + '_genefamily_cluster.txt'
    TIME = convert_usearch_result(tmp_uc.name, merged_txt, TIME, VERBOSE)
    # get dictionary gene2family
    gene2family = usearch_get_gene2family_dict(merged_txt, VERBOSE)
    # Add prefix clade:genefamID: to geneIDs in centroid.ffn sequence file
    #  to do: add also function (gene description)
    usearch_centroids_add_geneID_prefix(clade, gene2family, output_path, gene2description) 
    # clean up tmp files
    if not KEEP_UC:
        os.unlink(tmp_uc.name)
    if VERBOSE: print('[I] Remove usearch7 tmp results')
    os.remove(merged_txt)
        
    return gene2family, TIME
    
# ------------------------------------------------------------------------------
def check_usearch7(VERBOSE, PLATFORM='lin'):
    '''
    Check if Usearch 7 is installed
    '''
    try:
        if PLATFORM == WINDOWS:
            usearch7_path = subprocess.Popen(['where','usearch7'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        else: # Linux, Mac, ...
            usearch7_path = subprocess.Popen(['which','usearch7'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        usearch7_version = subprocess.Popen(['usearch7','--version'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        usearch7_version = usearch7_version.split()[1]
    except OSError as err:
        show_error_message(err)
        print('\n[E] Please, install Usearch 7\n    download from: http://drive5.com/usearch/')
        if VERBOSE:
            print('    Usearch 7 is required to merge gene sequences into gene-family cluster.\n')
        sys.exit(UNINSTALLED_ERROR_CODE)

    if VERBOSE:
        print('[I] Usearch v.7 is installed, version: ' + str(usearch7_version) + ', path: ' + str(usearch7_path).strip())

# ------------------------------------------------------------------------------
def check_bowtie2(VERBOSE, PLATFORM='lin'):
    '''
    Check if Bowtie2 is installed
    '''
    try:
        if PLATFORM == WINDOWS:
            bowtie2 = subprocess.Popen(['where', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        else: # Linux, Mac, ...
            bowtie2 = subprocess.Popen(['which', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        bowtie2_version = subprocess.Popen(['bowtie2', '--version'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        bowtie2_version = bowtie2_version.split()[2]
        if VERBOSE:
            print('[I] Bowtie2 is installed, version: ' + str(bowtie2_version) + ', path: ' + str(bowtie2).strip())
    except OSError as err:
        show_error_message(err)
        print('\n[E] Please, install Bowtie2.\n')
        if VERBOSE:
            print('    Bowtie2 is used to generate the .bt2 index files required in panphlan_map.py\n')
        sys.exit(UNINSTALLED_ERROR_CODE)        

# ------------------------------------------------------------------------------
def add_filename_to_geneIDs(pathgenefiles, tmp_path, VERBOSE):
    '''
    To get unique geneIDs across all genomes:
    1) Copy gene ffn files to TMP
    2) add filename as prefix to geneIDs "Filename:originalGenID"

    requires Biopython (Bio module)
    '''
    if VERBOSE:
        print('[I] To get unique geneIDs across genomes: add filename as prefix to geneIDs')

    # create new folder 'ffn_uniqueGeneIDs' in TMP
    new_ffn_folder = os.path.join(tmp_path,'ffn_uniqueGeneIDs','') # '' to get ending '/'
    os.makedirs(new_ffn_folder)
    
    # create new list of ffn files: new_pathgenefiles
    new_pathgenefiles = [os.path.join(new_ffn_folder,os.path.basename(f)) for f in pathgenefiles]

    for (ffn_in,ffn_out) in zip(pathgenefiles,new_pathgenefiles):
        filename = os.path.splitext(os.path.basename(ffn_in))[0]
        with open(ffn_out, 'w') as f_out:
            for seq in SeqIO.parse(open(ffn_in), 'fasta'):
                if seq.id == seq.name:
                    seq.name=''
                if seq.id == seq.description.split()[0]:    
                    seq.description=' '.join(seq.description.split()[1:])
                seq.id = filename + ':' + seq.id    
                r = SeqIO.write(seq, f_out, 'fasta')
                if r!=1:
                    sys.exit('[E] Error while writing sequence to ffn-file:\n    ' + ffn_out)    

    return new_pathgenefiles

# ------------------------------------------------------------------------------
def check_genomes(ffn_folder, fna_folder, VERBOSE):
    '''
    Check if genome files and .fna .fnn pairs are present and calculate expected runtime
    Result are lists of correct genome and gene file pairs 
    '''
    genomefiles = [f for f in os.listdir(fna_folder) if fnmatch(f,'*.'+FNA)]
    genefiles   = [f for f in os.listdir(ffn_folder) if fnmatch(f,'*.'+FFN)]

    if len(genomefiles)==0:
        print('\n[E] Cannot find any genome.fna file in folder:\n    ' + fna_folder)
        print('    Genome files need to end with .fna\n')
        sys.exit('Missing genome files')

    # check if all genome-gene file pairs exist, remove single genome or gene-files from list
    for f in genomefiles[:]: # need a copy [:], as we remove items from list
        path_genefile_ffn = os.path.join(ffn_folder, f.replace('.'+FNA,'.'+FFN) )
        if not os.path.exists(path_genefile_ffn):
            print('[W] Cannot find gene-file:\n    ' + path_genefile_ffn)
            print('    Excluding genome ' + f + ' from pangenome database')
            genomefiles.remove(f)
    for f in genefiles[:]:  # need a copy [:], as we remove items from list
        path_genomefile_fna = os.path.join(fna_folder, f.replace('.'+FFN,'.'+FNA) )
        if not os.path.exists(path_genomefile_fna):
            print('[W] Cannot find genome-file:\n    ' + path_genomefile_fna)
            print('    Excluding corresponding genes of file: ' + f + ' from pangenome database')    
            genefiles.remove(f)

    if not len(genomefiles) == len(genefiles):
        sys.exit('\nError in function check_genomes of panphlan_pangenome_generation.py\n')

    if VERBOSE:
        print('[I] Total number of genomes (genome.fna gene.ffn file pairs): ' + str(len(genomefiles)))
    if len(genomefiles)==0:
        print('\n[E] Cannot find any genome-gene pair of .fna .fnn files having the same filename: ID.fna ID.ffn')
        sys.exit('Missing genome-gene file pairs')

    # add full path to genefile list        
    pathgenomefiles = sorted([os.path.join(fna_folder,f) for f in genomefiles])
    pathgenefiles   = sorted([os.path.join(ffn_folder,f) for f in genefiles])
    
    # check usearch7 cluster depends on file order? (bvulgatus14) Yes, result depends on gene-file order 
    # pathgenomefiles = [ pathgenomefiles[i] for i in [2,1,0,3]]  # Original file order
    # pathgenefiles   = [ pathgenefiles[i]   for i in [2,1,0,3]]
    
    print('\nExpected runtime: ' + str(len(genomefiles)*20) + ' minutes (start time: ' + time.strftime("%b %d %Y %H:%M") + ')\n')
    if not VERBOSE:
        print('Use option --verbose to display progress information.\n')

    return pathgenomefiles, pathgenefiles

# ------------------------------------------------------------------------------
def clean_up(pathgenefiles, tmp_path, VERBOSE):
    '''
    Remove files not needed anymore
    1) copy of gene-sequence ffn files having prefix to geneIDs "Filename:originalGenID"
    2) TMP folder
    '''
    if VERBOSE:
        print('[I] Remove copies of gene-sequence ffn files from TMP/')
    for f in pathgenefiles:
        if tmp_path in f: # make sure we deleting in the TMP directory
            os.remove(f)
    
    if VERBOSE:
        print('[I] Remove TMP/ directory')
    os.rmdir(os.path.join(tmp_path,'ffn_uniqueGeneIDs'))
    os.rmdir(tmp_path)

# ------------------------------------------------------------------------------
def check_args():
    '''
    Check input arguments
    '''
    parser = PanPhlAnGenParser()
    args = vars(parser.parse_args())
    VERBOSE = args['verbose']

    if VERBOSE:
        print('\nPanPhlAn pangenome generation version '+__version__)
        print('Python version: ' + sys.version.split()[0])
        print('System: ' + sys.platform)

    # Check: FFN_FOLDER --------------------------------------------------------
    ipath = args['i_ffn']
    if not os.path.exists(ipath):
        show_error_message('Input folder -i_ffn does not exist.')
        sys.exit(INEXISTENCE_ERROR_CODE)
        
    ipath = os.path.abspath(ipath)
    ipath = os.path.join(ipath,'')
    args['i_ffn'] = ipath
    # if not ipath[-1] ==  '/':
    #    args['i_ffn'] = ipath + '/'
    if VERBOSE:
        print('[I] Input gene FFN folder: ' + args['i_ffn'])

    # Check: FNA_FOLDER --------------------------------------------------------
    ipath = args['i_fna']
    if not os.path.exists(ipath):
        show_error_message('Input folder -i_fna does not exist.')
        sys.exit(INEXISTENCE_ERROR_CODE)

    ipath = os.path.abspath(ipath)
    ipath = os.path.join(ipath,'')
    args['i_fna'] = ipath
    # if not ipath[-1] ==  '/':
    #    args['i_fna'] = ipath + '/'
    if VERBOSE:
        print('[I] Input genome FNA folder: ' + args['i_fna'])

    # Check: CLADE -------------------------------------------------------------
    if VERBOSE:
        print('[I] Clade or species name: ' + args['clade'])

    # Check: IDENTITY_PERCENATGE -----------------------------------------------
    identity_threshold_perc = args['th']
    if identity_threshold_perc < 0.0 or identity_threshold_perc > 100.0:
        args['th'] = 95.0
        if VERBOSE:
            print('[I] Invalid value for identity threshold percentage. Default value (95.0 %) has been set.')
    else:
        if VERBOSE:
            print('[I] Identity threshold percentage: ' + str(args['th']) + ' %.')

    # Check: OUTPUT_FOLDER -----------------------------------------------------
    opath = os.path.join(args['output'],'')
    # if not opath[-1] == '/':
    #    opath +=  '/'
    if not os.path.exists(os.path.dirname(opath)):
        os.makedirs(opath)
    args['output'] = opath
    if VERBOSE:
        print('[I] Output folder: ' + args['output'])

    # Check: TEMP_FOLDER ------------------------------------------------------
    tmp_path = os.path.join(args['tmp'],'')
    # if not tmp_path[-1] == '/':
    #    tmp_path +=  '/'
    if not os.path.exists(os.path.dirname(tmp_path)):
        os.makedirs(tmp_path)
    args['tmp'] = tmp_path
    if VERBOSE:
        print('[I] Temporary folder: ' + args['tmp'])

    return args

# ------------------------------------------------------------------------------
def main():
    # Check Python version
    if sys.hexversion < 0x02060000:
        print('Python version: ' + sys.version)
        sys.exit('Python versions older than 2.6 are not supported.')

    args = check_args()
    VERBOSE = args['verbose']
    KEEP_UC = args['uc']
    PLATFORM = sys.platform.lower()[0:3]
    
    TOTAL_TIME = time.time()
    TIME = time.time()
    
    # Check if software is installed
    if VERBOSE: print('\nSTEP 1. Checking required software installations...')
    bowtie2   = check_bowtie2(VERBOSE, PLATFORM)  # for generating .bt2 index files
    usearch7  = check_usearch7(VERBOSE, PLATFORM) # for getting gene-family cluster
    
    # check input genome and gene files
    pathgenomefiles, pathgenefiles = check_genomes(args['i_ffn'], args['i_fna'], VERBOSE)
    pathgenefiles = add_filename_to_geneIDs(pathgenefiles, args['tmp'], VERBOSE)
    gene2genome, gene2description = gene2genome_mapping(pathgenefiles, VERBOSE)

    # Get gene families cluster (usearch7)
    if VERBOSE: print('\nSTEP 2. Generating gene families cluster (usearch7) ...')
    gene2family, TIME = usearch_clustering(pathgenefiles,args['th'],args['clade'],
                                           args['output'],args['tmp'],KEEP_UC,gene2description,TIME,VERBOSE)

    # Get pangenome and bowtie2 index file
    if VERBOSE: print('\nSTEP 3. Getting pangenome file...')
    TIME = pangenome_generation(pathgenomefiles, pathgenefiles, gene2family,
                                args['clade'], args['output'], gene2genome, TIME, VERBOSE)
    TIME = create_bt2_indexes(pathgenomefiles, args['clade'], args['output'],
                              args['tmp'], TIME, VERBOSE)
    
    clean_up(pathgenefiles, args['tmp'], VERBOSE)
    end_program(time.time() - TOTAL_TIME)

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
