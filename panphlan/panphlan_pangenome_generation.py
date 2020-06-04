#!/usr/bin/env python

# ==============================================================================
# PanPhlAn v1.2.2: PANgenome-based PHyLogenomic ANalysis
# Detecting and characterizing strains in metagenomic samples
#
# Matthias Scholz, Doyle V. Ward, Edoardo Pasolli, Thomas Tolio, Moreno Zolfo,
# Francesco Asnicar, Duy Tin Truong, Adrian Tett, Ardythe L. Morrow, and Nicola Segata.
# Strain-level microbial epidemiology and population genomics from shotgun metagenomics.
# Nature Methods, 13, 435-438, 2016.
#
# For help type: ./panphlan_pangenome_generation.py -h
#
# Example:
# ./panphlan_pangenome_generation.py --i_gff ncbi_download_gff/ --i_fna ncbi_download_fna/ -c speciesname
#
# https://bitbucket.org/CibioCM/panphlan
# ==============================================================================

from __future__ import print_function # to give equal outputs in python2 and python3
from __future__ import with_statement
from argparse import ArgumentParser
from collections import defaultdict
import os, subprocess, sys, tempfile, time
from fnmatch import fnmatch
import re # for gene genome mapping
import shutil
import gzip

if '--i_gff' in sys.argv:
    try:
        from BCBio import GFF
    except ImportError as err:
        print('\n[Error]',err)
        print('\n[Error] Please install the BCBio package of Python')
        print('  https://pypi.python.org/pypi/bcbio-gff')
        print('  pip install bcbio-gff\n\n')
        sys.exit(2)
    try:
        from Bio import SeqIO
        from Bio.Seq import UnknownSeq
        from Bio.SeqRecord import SeqRecord
    except ImportError as err:
        print('\n[Error]',err)
        print('\n[Error] Please install the BCBio package of Python')
        print('  http://biopython.org\n  https://pypi.python.org/pypi/biopython/')
        print('  pip install biopython\n\n')
        sys.exit(2)

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError as err:
    print('\n[E]',err)
    print('\n[E] Please install Biopython.')
    print('    The "Bio" module is required for extracting gene locations')
    print('    by mapping genes against their genome.\n')
    sys.exit(2)


__author__  = 'Matthias Scholz, Moreno Zolfo, Thomas Tolio, Nicola Segata (panphlan-users@googlegroups.com)'
__version__ = '1.2.3.7'
__date__    = '27 February 2018'


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
        self.add_argument('--i_ffn',         metavar='INPUT_FFN_FOLDER',     type=str,   default=False,      help='Folder containing the .ffn gene sequence files')
        self.add_argument('--i_fna',         metavar='INPUT_FNA_FOLDER',     type=str,   default=False,      help='Folder containing the .fna genome sequence files')
        self.add_argument('--i_gff',         metavar='INPUT_GFF_FOLDER',     type=str,   default=False,      help='Folder containing the .gff gene annotation files')
        self.add_argument('-c','--clade',    metavar='CLADE_NAME',           type=str,   default=False,      help='Name of the species pangenome database, for example: -c ecoli17')
        self.add_argument('-o','--output',   metavar='OUTPUT_FOLDER',        type=str,   default='database', help='Result folder for all database files')
        self.add_argument('--roary_dir',     metavar='ROARY_FOLDER',         type=str,   default=False,      help='Use pre-processed Roary pangenome clustering (instead of usearch): Folder containing gene family cluster results of Roary based on gff')
        self.add_argument('--th',            metavar='IDENTITY_PERCENATGE',  type=float, default=95.0,       help='Threshold of gene sequence similarity (in percentage), default: 95.0 %%.')
        self.add_argument('--tmp',           metavar='TEMP_FOLDER',          type=str,   default='TMP_panphlan_db', help='Folder for temporary files, default: TMP_panphlan_db')
        self.add_argument('--uc',            action='store_true',                                            help='Keep all usearch7 output files')
        self.add_argument('--verbose',       action='store_true',                                            help='Show progress information')
        self.add_argument('-v', '--version', action='version',   version="PanPhlAn version "+__version__+"\t("+__date__+")", help='Prints the current PanPhlAn version and exits')

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

def openread( filename, mode = "r" ):
    # open file for reading, allow both compressed or uncompressed versions
    if filename.endswith('.gz'):
        return gzip.open(filename,mode+"t") #  to force text mode 't' in python3 (read gff.gz files)
    elif filename.endswith('.bz2'):
        import bz2 # bz2 not always available on servers: import only when needed
        return bz2.BZ2File(filename,mode)
    else:
        return open(filename,mode)

# ------------------------------------------------------------------------------
# MAJOR FUNCTIONS
# ------------------------------------------------------------------------------
def create_bt2_indexes(path_genome_fna_files, clade, output_path, tmp_path, TIME, VERBOSE):
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
        for f in path_genome_fna_files:
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
def write_pangenome_file(gene2loc, gene2family, gene2genome, output_path, clade, VERBOSE):
    '''
    Create the pangenome database file combining all the information from
    gene mappings (location (contig, from, to), family and genome)

    Result: pangenome file (tab-separated):
        geneFamily | geneID | genomeName(filename) | contigID | start | stop
    '''
    pangenome_csv = output_path + 'panphlan_' + clade + '_pangenome.csv'
    missing_geneIDs_in_familycluster = []
    n=0
    with open(pangenome_csv, mode='w') as ocsv:
        genes_list = sorted(gene2loc.keys())
        for gene in genes_list:
            if gene in gene2family:
                n=n+1
                ocsv.write(gene2family[gene] + '\t' + gene + '\t' + gene2genome[gene] + '\t' + gene2loc[gene][0] + '\t' + str(gene2loc[gene][1]) + '\t' + str(gene2loc[gene][2]) + '\n')
            else:
                missing_geneIDs_in_familycluster.append(gene)
                # print('[W] Could not find gene in gene-family cluster result (dict gene2family): '   + gene)
                # print('    Check presence of gene in file: usearch7_species_cluster.uc, option --uc')
    maxN=min(5,len(missing_geneIDs_in_familycluster))
    if maxN>0:
        if VERBOSE: print('    ' + str(len(missing_geneIDs_in_familycluster)) + ' of ' + str(len(genes_list)) + ' genes are not present in gene-family clusters: ')
        for g in missing_geneIDs_in_familycluster[0:maxN]:
            if VERBOSE: print('      '+ g)
        if VERBOSE: print('      o o o')
        print('      If many genes are missing, check file: usearch7_species_cluster.uc, option --uc, or Roary gene_presence_absence.csv')
    if VERBOSE: print('[I] Pangenome file has been generated ('+str(n)+' genes):\n    ' + pangenome_csv)
# ------------------------------------------------------------------------------
def gff_add_genome_seq(gff_folder, fna_folder, VERBOSE):
    '''
    Add genome sequence from .fna file to .gff file (gff's from NCBI don't have genome sequence included)
    PanPhlAn (and Roary) work with genome sequence included as in Prokka gff's.
    1) read all .gff, gff.gz, and gff,bz2 files (if both compressed and uncompressed version exist, take uncompressed)
    2) check for missing genome sequence (not missing -> don't add fna, copy orig gff to new folder)
    3) add sequence from .fna file
    '''
    # read .gff or .gff.gz
    path_gff_files = [os.path.join(gff_folder,f) for f in os.listdir(gff_folder) if f.endswith(('.gff','.gff.gz','.gff.bz2'))]
    if len(path_gff_files) < 1: print('\nERROR: Cannot find .gff files in folder\n  '+gff_folder+'\n')
    # if both raw .gff and compressed gff.gz are present, keep only uncompressed .gff in list
    path_gff_files =  [f for f in path_gff_files if not os.path.splitext(f)[0] in path_gff_files]

    # create new gff folder 'gff_added_fna' (genome sequence included in gff files)
    new_gff_folder = os.path.join(os.getcwd(),'gff_added_fna','') # '' to get ending '/'
    if os.path.exists(new_gff_folder):
        print('\nERROR: New gff directory exist already:\n  '+new_gff_folder+'\n  Please rename or remove the directory and try again.\n')
        sys.exit(2)
    else:
        os.makedirs(new_gff_folder)

    for gff_file in path_gff_files:
        # read and check presence of seq
        genomeID = os.path.basename(gff_file).replace('.gff.gz','').replace('.gff.bz2','').replace('.gff','') # to allow dots in genome-name
        # genomeID = os.path.splitext(os.path.basename(gff_file))[0] # not working for compressed files, gff remains
        if VERBOSE: print('[I] ' + genomeID + '\n    Read gff file:\n    ' + gff_file)

        try:
            contig = next(GFF.parse(openread(gff_file)))
        except Exception as err:
            print('BCBio GFF parse error')
            print('Error:',err)
            print('Current Python version: ' + sys.version)
            print('Please use Python3\n\n')
            sys.exit(1)
        # GFF.parse does not work under conda Python version: 2.7.15, Feb 28 2019
        # when applied to gff file, fasta seq included, single contig
        # contig = next(GFF.parse(openread(gff_file)))

        # check if FASTA genome sequence is already present
        gff_version = int(contig.annotations['gff-version'][0])
        if VERBOSE: print('    gff-version: ' + str(gff_version))
        if type(contig.seq)==UnknownSeq: # type check for Bio.Seq.UnknownSeq
            fna_present = False
        else:
            fna_present = True
            print('\n WARNING: genome sequence is already present in gff file (--fna not needed): ' + genomeID + '\n')

        # search fna file (FASTA genome sequence)
        fna_search_files = [os.path.join(fna_folder, genomeID + ending) for ending in ('.fna','.fna.gz','.fna.bz2')]
        fna_files_exist = [f for f in fna_search_files if os.path.exists(f)]
        fna_file = fna_files_exist[0] if fna_files_exist else None

        # add fna to gff file, and save in new gff folder
        gff_out = os.path.join(new_gff_folder,genomeID+'.gff')
        if fna_present: # copy gff file, if sequence already included
            print('    Keep gff file unchanged, copy gff file to\n    ' + gff_out)
            with openread(gff_file) as f_in, open(gff_out, 'w') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        elif fna_file: # add fna genome sequence, if fna file exist
            if VERBOSE: print('    Add fna sequence to gff file\n    ' + gff_out)
            with open(gff_out, 'w') as f_out:
                with openread(gff_file) as f_in:
                    for line in f_in:
                        f_out.write(line)
                f_out.write('##FASTA\n') # add '##FASTA' separator
                with openread(fna_file) as f_in: # add fna genome file
                    for line in f_in:
                        f_out.write(line)
        else:
            print('\n WARNING: Cannot find genome file: ' + os.path.join(fna_folder,genomeID+'.fna\n'))

    num_new_gff = len([os.path.join(new_gff_folder,f) for f in os.listdir(new_gff_folder) if f.endswith(('.gff'))])
    if num_new_gff==0: sys.exit('\n\n ERROR: Could not add fna genome sequence to any gff file.\n')
    print('\n[I] ' + str(num_new_gff) + ' gff files that include the .fna genome sequence are written to folder:\n    ' + new_gff_folder + '\n')
    return new_gff_folder
# ------------------------------------------------------------------------------
def read_gff_write_fna_ffn(gff_folder, VERBOSE):
    '''
    gene2loc,gene2genome,gene2description,gene2gffdata,path_fna_files,path_ffn_files = read_gff_write_fna_ffn(gff_folder)
    main parts developed by Moreno Zolfo
    a) read gff-files and extract gene location and annotation
    b) write genome fna files (for bowtie2 index database)
    c) write gene ffn files (for usearch7 clustering)
    gene2loc    = {geneID : (contig,start,stop)}
    gene2genome = {geneID : genomefilename}
    gene2description = {geneID : description}
    '''
    # read gff and prepare output folder fna ffn
    path_gff_files = [os.path.join(gff_folder,f)for f in os.listdir(gff_folder) if fnmatch(f,'*.gff')]
    if len(path_gff_files) < 1:
        print('\nERROR: Cannot find .gff files in folder\n  '+gff_folder+'\n')
        # sys.exit(2)
    # create new ffn fna folder 'fna_from_gff' 'ffn_from_gff'
    new_fna_folder = os.path.join(os.getcwd(),'fna_from_gff','') # '' to get ending '/'
    new_ffn_folder = os.path.join(os.getcwd(),'ffn_from_gff','')
    if os.path.exists(new_fna_folder) or os.path.exists(new_ffn_folder):
        print('\nERROR: Gene ffn or genome fna directory exist already:\n  '+new_fna_folder+'\n  '+new_ffn_folder+'\n  Please rename or remove the directories and try again.\n')
        sys.exit(2)
    else:
        os.makedirs(new_fna_folder)
        os.makedirs(new_ffn_folder)

    # loop over all gff files
    gene2loc = defaultdict(tuple)
    gene2genome      = {} # {geneID : genomefilename}
    gene2description = {} # {geneID : description}
    gene2gffdata = defaultdict(dict)
    for gff_file in path_gff_files:
        genomeID = os.path.splitext(os.path.basename(gff_file))[0] # get genome filename without extension (allow dots in genome-name)
        print('[I] ' + genomeID)
        if VERBOSE: print('    Read gff file:\n    ' + gff_file)
        fna_genome_seq  = []
        ffn_gene_seq    = []
        gff_version = 0
        t1 = ''
        for contig in GFF.parse(open(gff_file)):
            if type(contig.seq)==UnknownSeq: # type check for Bio.Seq.UnknownSeq
              print('\n\n Cannot find genome sequence in gff file: \n ' + gff_file)
              print(' Used option --fna to provided genome sequences.')
              print('\n ERROR: Missing genome sequence.\n\n')
              sys.exit(2)
            if not gff_version:
                gff_version = int(contig.annotations['gff-version'][0])
                if VERBOSE: print('    gff-version: ' + str(gff_version))
            # create pure contig-record without contig.features
            # add genomeID(filename) to contigID to avoid duplicated contig-names
            fna_contigID = contig.id
            if not fna_contigID.startswith(genomeID):
                fna_contigID = genomeID + ':' + fna_contigID
            contig_record = SeqRecord(contig.seq, id=fna_contigID, name='', description='')
            fna_genome_seq.append(contig_record)
            for t in contig.features:
                if not t1: t1=t # keep first record
                if 'locus_tag' in t.qualifiers:
                    # get gene sequence and final geneID -> genomeID:geneID:start-stop)
                    locus_tag = t.qualifiers['locus_tag'][0]
                    g_start = int(t.location.start) # start is python 0-based (-1 compared to gff-file)
                    g_stop  = int(t.location.end)
                    gene_seq = contig.seq[g_start:g_stop] # start is still 0-based
                    g_start = g_start + 1 # correct 0-based to 1-based (as in gff), to use in geneID and pangenome file
                    if t.strand == -1:
                        gene_seq = gene_seq.reverse_complement()
                        s_start,s_stop = str(g_stop),str(g_start)
                    elif t.strand == 1:
                        s_start,s_stop = str(g_start),str(g_stop)
                    else:
                        print('Error: no location-strand information, gene: ' + t.qualifiers['locus_tag'][0])
                        sys.exit(2)
                    ffn_geneID = fna_contigID + ':' + s_start + '-' + s_stop # contigID already contains genomefilename prefix
                    gene2loc[ffn_geneID]    = (fna_contigID, min(g_start,g_stop), max(g_start,g_stop)) # to always have start < stop
                    gene2genome[ffn_geneID] = genomeID

                    # get gff annotation data
                    if not locus_tag.startswith(genomeID): locus_tag = genomeID + ':' + locus_tag # REF_wMel_A:gene_00258
                    gene2gffdata[ffn_geneID]['locus_tag']=locus_tag
                    if 'gene' in t.qualifiers:
                        gene2gffdata[ffn_geneID]['gene'] = t.qualifiers['gene'][0]
                    if t.sub_features and 'protein_id' in t.sub_features[0].qualifiers: # protein_in in sub_features (NCBI & Tin's gff)
                        gene2gffdata[ffn_geneID]['protein_id'] = t.sub_features[0].qualifiers['protein_id'][0]
                    if 'product' in t.qualifiers: # wolbachia prokka gff
                        gene2gffdata[ffn_geneID]['product'] = t.qualifiers['product'][0] # NCBI and Tin's gff
                    elif t.sub_features and 'product' in t.sub_features[0].qualifiers:
                        gene2gffdata[ffn_geneID]['product'] = t.sub_features[0].qualifiers['product'][0]
                    if 'eC_number' in t.qualifiers: # Enzyme Commission number # wolbachia prokka gff
                        gene2gffdata[ffn_geneID]['eC_number'] = t.qualifiers['eC_number'][0]
                    elif t.sub_features and 'eC_number' in t.sub_features[0].qualifiers: # Tin's gff (not NCBI)
                        gene2gffdata[ffn_geneID]['eC_number'] = t.sub_features[0].qualifiers['eC_number'][0]
                    # use 'product' as gene sequence description
                    gene2description[ffn_geneID] = gene2gffdata[ffn_geneID].get('product', '')
                    # gene_annotations_string = '##'.join(filter(None, gene_annotations ))
                    # gene_annotations_string = gene_annotations_string.replace(' ','_')

                    # create gene sequence record
                    ffn_description = ' '.join([locus_tag, gene2description[ffn_geneID]])
                    gene_record = SeqRecord(gene_seq, id=ffn_geneID, name='', description=ffn_description)
                    ffn_gene_seq.append(gene_record)


        # write gene.ffn and genome.fna files (genomefilename added to sequence ID's)
        fna_out = os.path.join(new_fna_folder,genomeID+'.fna')
        ffn_out = os.path.join(new_ffn_folder,genomeID+'.ffn')
        numGenes = len(ffn_gene_seq)
        if VERBOSE: print('    Number of genes (locus tags):  ' + str(numGenes))
        if VERBOSE: print('    Write genome .fna and gene .ffn files:\n    ' + fna_out + '\n    '  + ffn_out)
        SeqIO.write(fna_genome_seq, fna_out,'fasta')
        SeqIO.write(ffn_gene_seq  , ffn_out,'fasta')

    if numGenes==0:
        print('\n\nExample of detected gff feature ("locus_tag" missing):\n',t1.qualifiers, '\n') # print last gff entry
        print('\n  ERROR: Could not find any gene features in gff file (no "locus_tag"s)')
        print('  ' + gff_file)
        print('\n  Please use gene feature .gff files (not CDS.gff). Try Prokka for predicting gene location .gff files from .fna genomes.\n')
        sys.exit('\n ERROR: Could not find any gene feature in gff file\n')
    # get list of all fna ffn files
    path_fna_files = [os.path.join(new_fna_folder,f) for f in os.listdir(new_fna_folder) if fnmatch(f,'*.fna')]
    path_ffn_files = [os.path.join(new_ffn_folder,f) for f in os.listdir(new_ffn_folder) if fnmatch(f,'*.ffn')]
    print('\n   Extracted gene .ffn and genome .fna sequence files are in folders:\n   ' + new_ffn_folder + '\n   '+ new_fna_folder + '\n')

    # check for duplicated locus tags,
    # Give only warning, not error as converting to ffn and fna is possible.
    # For usearch approach (ffn,fna) no problem, seqIDs are unique by added filename prefix, but Roary is using orig gff (no prefix).
    all_locus_tags = [ gene2gffdata[k].get('locus_tag') for k in gene2gffdata]
    if not len(all_locus_tags) == len(set(all_locus_tags)):
        print('\nWARNING: gff files have duplicated locus tags (across all gff files)')
        print('         Please correct before using gff files in Roary gene clustering\n\n')

    return gene2loc,gene2genome,gene2description,gene2gffdata,path_fna_files,path_ffn_files
# ------------------------------------------------------------------------------
def get_gene_locations(path_genome_fna_files, path_gene_ffn_files, VERBOSE):
    '''
    Extract gene locations from gene-identifier or blast-like search.

    Get gene locations: Read all .ffn files to extract start and stop location from gene-name,
    If location cannot be extracted from gene-name, use a blast-like mapping of genes against their genomes.

    Gene-locations are returned as dictionary {geneID:(contig,start,stop)}

    Requires: Biopython module
    '''
    if VERBOSE: print('[I] Get gene locations and contigs for each gene.')
    gene2loc = defaultdict(tuple)
    for (genomefile, genefile) in zip(path_genome_fna_files,path_gene_ffn_files):
        if VERBOSE:
            print('[I] genomefile: ' + genomefile)
            print('    genefile: '   + genefile)
        try: # extract gene-location from geneIDs
            if VERBOSE:
                print('    Extract gene-location from geneIDs')
            gene_location_check_done = False
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
                locinfo=r.id.split(':')[-1] # try to get location info from after last ':'
                pos1 = int(locinfo.split('-')[0].replace('c',''))
                pos2 = int(locinfo.split('-')[-1].replace('c',''))
                # contig = r.id.split(':')[-2]
                contig = ':'.join(r.id.split(':')[:-1]) # take all before last ':' as contigID
                start, stop = min(pos1, pos2), max(pos1, pos2) # to always have start < stop
                gene2loc[r.id] = (str(contig), start, stop)
                if not gene_location_check_done: # double check first gene for correct start-stop location
                    if not check_gene_location(genomefile, genefile, r.id, gene2loc[r.id], VERBOSE):
                        raise IndexError('Location extracted from geneID is not correct')
                    gene_location_check_done = True
        except (IndexError, ValueError) as err: # alternatively, run BLAST-like python gene-genome mapping to get locations
            if VERBOSE:
                print('    Expected geneID format: ">contigID:start-stop" (1-based)')
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
def check_gene_location(genomefile,genefile,geneID,geneLocation, VERBOSE):
    '''
    check if gene-location from geneID gives correct gene sequence
    1) extract gene sequence from genome, based on geneLocation (contig, start, stop)
    2) extract gene sequence from ffn-gene-file
    3) check: both sequences have to be identical!
    '''
    if VERBOSE: print('    Check if location extracted from geneID is correct, gene: ' + geneID)
    for contig in SeqIO.parse(open(genomefile, mode='r'), 'fasta'):
        if contig.id==geneLocation[0]:
            start=geneLocation[1]
            stop =geneLocation[2]
            gene_from_fna        = contig.seq[start-1:stop]
            gene_from_fna_0based = contig.seq[start:stop] # 0-based bedtools format
    for gene in SeqIO.parse(open(genefile, mode='r'), 'fasta'):
        if gene.id==geneID:
            gene_from_ffn = gene.seq
    if gene_from_ffn == gene_from_fna:
        return True
    if VERBOSE: print('    Sequence location in genome file is not identical with gene sequence!')
    if gene_from_ffn == gene_from_fna_0based:
        if VERBOSE: print('    GeneID might follow 0-based (start-1) BED-tools format, PanPhlAn requires 1-based locations as in .gff files.')
        if VERBOSE: print('    Use panphlan to extract .ffn gene sequences using .gff gene annotations')
    return False
# ------------------------------------------------------------------------------
def get_contigs(path_genome_fna_files, VERBOSE):
    '''
    contig2genome = get_contigs(path_genome_fna_files)
    Map each contig to its genome (filename)
    Used to
     - check if contig-name from gene2loc (contig,start,stop) is valid
     - check for duplicated contig-names across all genoems
    requires Biopython (Bio module)
    '''
    # genome2contigs = defaultdict(set) # { genome : contig-set }
    contig2genome  = {} # {contig : genomefilename}
    if VERBOSE: print('[I] Get contig-IDs from all genome fna files')
    for f in path_genome_fna_files:
        # genome = os.path.basename(f).split('.')[0] # get genome filename without extension
        genome = os.path.splitext(os.path.basename(f))[0] # get genome filename without extension (allow dots in genome-name)
        # if VERBOSE: print('    Genome: ' + genome)
        for r in SeqIO.parse(open(f, mode='r'), 'fasta'):
            contigID=r.id
            if contigID not in contig2genome:
                contig2genome[contigID] = genome
                # genome2contigs[genome].add(contigID)
            else:
                print('\n\n[I] Duplicated contig-ID: \n      ' + contigID + '\n    found in genome files: \n      ' + contig2genome[contigID] + '\n      ' + genome)
                print('\n    ERROR: Duplicated contig-IDs\n\n')
                sys.exit(2)
    return contig2genome
#-------------------------------------------------------------------------------
def check_for_valid_contigIDs(gene2loc,contig2genome):
    '''
    Check if contigID from gene2loc (contig,start,stop) is valid
    '''
    for gene in gene2loc:
        gene_contigID=gene2loc[gene][0]
        if gene_contigID not in contig2genome:
            print('\n\n[I] Gene-ID:   ' + gene)
            print('    Contig-ID: '     + gene_contigID)
            print('\n    ERROR: Contig-ID extracted from geneID is not valid (not present in any genome fna file)\n')
            sys.exit(2)
    return True
# ------------------------------------------------------------------------------
def gene2genome_mapping(path_gene_ffn_files, VERBOSE):
    '''
    Map each gene to its genome (genome-filename)

    requires Biopython (Bio module)
    '''
    gene2genome     = {} # {geneID : genomefilename}
    gene2description = {} # {geneID : description}
    if VERBOSE: print('[I] Get gene-genome list')
    for f in path_gene_ffn_files:
        # genome = os.path.basename(f).split('.')[0] # get genome filename without extension
        genome = os.path.splitext(os.path.basename(f))[0] # get genome filename without extension (allow dots in genome-name)
        if VERBOSE: print('    Genome: ' + genome)
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
# --- Roary --------------------------------------------------------------------
def read_roary_gene_clustering(roary_folder, VERBOSE):
    '''
    Read roary clustering result
    1) Search for file: 'gene_presence_absence.csv'
    2) Read file and create dict mappings:
         gene2family       = {geneID   : familyID}
         family2annotation = {familyID : function}
    '''
    # check if roary cluster file exist
    roary_cluster_csv = os.path.join(roary_folder,'gene_presence_absence.csv')
    if not os.path.exists(roary_cluster_csv):
        print('\nERROR: Cannot find Roary gene_presence_absence.csv file:\n  '+roary_cluster_csv+'\n')
        sys.exit(2)
    # read roary file
    if VERBOSE: print('    Read Roary gene cluster "gene_presence_absence.csv" file:\n    '+roary_cluster_csv+'\n')
    gene2family       = {}
    family2annotation = {}
    roary_genomeIDs = []
    NotSampleColumns = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']
    with open(roary_cluster_csv, mode='r') as f:
        headerline=f.readline().strip().strip('"').split('","')
        if 'Gene' not in headerline:
            print('\nERROR: Cannot find column "Gene" in Roary gene_presence_absence.csv file:\n  '+roary_cluster_csv+'\n')
            sys.exit(2)
        if 'Annotation' not in headerline:
            print('\nERROR: Cannot find column "Annotation" in Roary gene_presence_absence.csv file:\n  '+roary_cluster_csv)
            sys.exit(2)
        for line in f:
            # time.sleep(0.1) # <<<<<<<<<<<<<<<<<
            l=line.strip().strip('"').strip('",').split('","')
            genefamily=''
            for (i,h) in zip(l,headerline):
                # print(i + ' - ' + h)
                if h=='Gene': genefamily = i
                if h=='Annotation': family2annotation[genefamily] = i
                if h not in NotSampleColumns:
                    genomeID = h
                    geneID = i
                    if not genomeID in roary_genomeIDs: roary_genomeIDs.append(genomeID)
                    if geneID: # if not empty (gene-family present in reference genome)
                        # add filename prefix in front of gene-name
                        if not geneID.startswith(genomeID): geneID = genomeID + ':' + geneID
                        # print(genomeID + ' - ' + geneID)
                        gene2family[geneID]=genefamily
    return gene2family, family2annotation, roary_genomeIDs
# ------------------------------------------------------------------------------
def convert_roary_geneIDs(roary_gene2family,gene2gffdata, VERBOSE):
    '''
    Covert Roary geneIDs (gff) to PanPhlAn geneIDs (contig based, former NCBI format)
    example: REF_wMel_A:gene_00793 (Roary) --> REF_wMel_A:NC_002978.6:153-1535 (PanPhlAn)
    - remove Roary geneIDs, not present in gff-input (gene2gffdata) and vice versa
    '''
    # get mapping from Roary locus-tags to PanPhlAn contig-based geneIDs, 'REF_wRi_A:gene_00050': 'REF_wRi_A:NC_012416.1:59109-59450'
    locustag2gene = dict( (gene2gffdata[k].get('locus_tag', 'NaN'), k) for k in gene2gffdata.keys() )
    # create new gene2family dict, using panphlan contig-based geneIDs
    gene2family={}
    missingIDs_in_roary = []
    for g in gene2gffdata.keys():
        locustag = gene2gffdata[g].get('locus_tag', 'NaN')
        if locustag in roary_gene2family.keys(): # 'REF_wAna_A:gene_00962': 'dnaE1'
            gene2family[g] = roary_gene2family[locustag]
        else:
            missingIDs_in_roary.append(locustag)
    maxN=min(5,len(missingIDs_in_roary))
    if maxN>0: # print top 5 missing genes
        if VERBOSE: print('    ' + str(len(missingIDs_in_roary)),'genes are not present in Roary clustering: ')
        for ltag in missingIDs_in_roary[0:maxN]:
            if VERBOSE: print('      '+ ltag + '\t' + locustag2gene[ltag])
        if VERBOSE: print('      o o o')
    num_roray_gff_hits = len(gene2family)
    if VERBOSE: print('    Number of total geneIDs (present in Roary and gff): ' + str(num_roray_gff_hits))
    if num_roray_gff_hits < 1:
        sys.exit('\n\n ERROR: Could not match Roary geneIDs with gff-file geneIDs.\n')
    return gene2family, locustag2gene
#-------------------------------------------------------------------------------
def read_roary_centroids(roary_folder, clade, output_path, locustag2gene, family2annotation, VERBOSE):
    '''
    1) read Roary centroid-sequence file pan_genome_reference.fa
    2) rename geneIDs into PanPhlAn contig based geneIDs
    3) save panphlan_species_centroid file (only genomes present in gff input)
    4) return mapping from genefamily-clusterID to gentroid-geneIDs
    '''
    # check if roary centroid sequence file exist
    panphlan_centroids_ffn   = os.path.join(os.path.abspath(output_path),'panphlan_' + clade + '_centroids.ffn')
    roary_centroids_ffn      = os.path.join(roary_folder,'pan_genome_reference.fa')
    if not os.path.exists(roary_centroids_ffn):
        print('\nERROR: Cannot find Roary centroid sequence file: "pan_genome_reference.fa"')
        print('       ' + roary_centroids_ffn)
        print('       Please run Roary with option "-e"')
        sys.exit('\nERROR: Cannot fine Roary pan_genome_reference.fa\n\n')
    # read roary centroids
    if VERBOSE: print('\n    Convert Roary centroid sequence file into panphlan centroid file\n    Roary:    '+roary_centroids_ffn+'\n    PanPhlAn: '+panphlan_centroids_ffn+'\n')
    family2centroidGeneID = {}

    roary_locustags_missing_in_gff = []
    with open(panphlan_centroids_ffn, 'w') as f_out:
        for seq_record in SeqIO.parse(open(roary_centroids_ffn, mode='r'), 'fasta'):
            seq_record.description=' '.join(seq_record.description.split()[1:]) # remove .id from .description record (remove all before first space)
            locustag         = seq_record.id
            genefamily       = seq_record.description
            roary_annotation = family2annotation[genefamily]
            if locustag in locustag2gene.keys():
                geneID = locustag2gene[locustag]
                family2centroidGeneID[genefamily] = geneID
                # write new fasta file
                # clade + ':' + genefamID + ':' + seq.id
                seq_record.id = clade + ':' + genefamily + ':' + geneID
                seq_record.description = locustag + ' ' + roary_annotation
                r=SeqIO.write(seq_record, f_out, 'fasta')
                if r!=1: print('Error while writing sequence:  ' + seq_record.id)
            else: # add all centroids, even when nof in gff, but other genes of the group may exist in gff
                seq_record.id = clade + ':' + genefamily + ':' + 'Centroid_gene_not_in_gff_files'
                seq_record.description = locustag + ' ' + roary_annotation
                r=SeqIO.write(seq_record, f_out, 'fasta')
                if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                roary_locustags_missing_in_gff.append(locustag)
    if VERBOSE: print('    ' + str(len(roary_locustags_missing_in_gff)),' Roary genes are not present in gff files: ')
    maxN=min(5,len(roary_locustags_missing_in_gff))
    for ltag in roary_locustags_missing_in_gff[0:maxN]:
        if VERBOSE: print('      '+ ltag)
    if VERBOSE: print('      o o o')
    return family2centroidGeneID
# ------------------------------------------------------------------------------
def reject_genomes_not_in_roary(path_genome_fna_files, roary_genomeIDs):
    '''
    Remove genomes not present in Roary clustering, to avoid having unconsidered genomes in bowtie2 index database
    '''
    print('\n    Check presence of .fna input genomes in Roary cluster result')
    # Check GFF in Roary
    for f in path_genome_fna_files[:]:  # need a copy [:], as we remove items from list
        genome = os.path.splitext(os.path.basename(f))[0]
        if genome not in roary_genomeIDs:
            print('    [W] "' + genome + '" excluded from bowtie2 index (not present in Roary)')
            path_genome_fna_files.remove(f)
    # check Roary in GFF
    genomefiles = [os.path.splitext(os.path.basename(f))[0] for f in path_genome_fna_files]
    for g in roary_genomeIDs:
        if g not in genomefiles:
            print('    [W] "' + g + '" missing in .fna genome files')
    if len(genomefiles)==0: sys.exit('\n Error: Roary cluster results does not match with input fna genomes\n')
    return path_genome_fna_files
# ------------------------------------------------------------------------------
def write_annotations_gff(clade, output_path, family2centroidGeneID, gene2gffdata, VERBOSE):
    '''
    Write gene-family annotation file, based on gff metadata
    '''
    annotation_file = output_path + 'panphlan_' + clade + '_annotations.csv'
    if VERBOSE: print('    Write gene-family annotations based on gff metadata\n    ' + annotation_file)
    with open(annotation_file, mode='w') as ocsv:
        genefam_list = sorted(family2centroidGeneID.keys())
        ocsv.write('Gene_family' +'\t'+ 'Centroid_gene_ID' +'\t'+ 'GFF_locus_tag' +'\t'+ 'Product' +'\t'+ 'EC_number' +'\t'+ 'Protein_id' + '\n')
        for gfam in genefam_list:
            geneID   = family2centroidGeneID[gfam]
            # 'protein_id' contains no value in Prokka gff's, but in NCBI gff's
            # if gene2gffdata[geneID].get('protein_id', ''): print(gene2gffdata[geneID].get('protein_id', ''))
            ocsv.write(gfam +'\t'+ geneID +'\t'+ gene2gffdata[geneID].get('locus_tag', '') +'\t'+ gene2gffdata[geneID].get('product', '') +'\t'+ gene2gffdata[geneID].get('eC_number', '') +'\t'+ gene2gffdata[geneID].get('protein_id', '') + '\n')
    return True
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
    family2centroidGeneID = {}
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
                family2centroidGeneID[genefamID]=seq.id # for annotation file
                seq.id = clade + ':' + genefamID + ':' + seq.id
                seq.name=''
                r = SeqIO.write(seq, f, 'fasta')
                if r!=1:
                    sys.exit('[E] Error while writing centroid sequence:  ' + seq.id)
        os.remove(centroids_orig_ffn)
    return family2centroidGeneID
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
def usearch_sortbylength(path_gene_ffn_files, tmp_path, TIME, VERBOSE):
    '''
    Merge all the gene-sequence FFN files into a unique one, then sort by length
    '''

    tmp_ffn        = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.ffn',        dir=tmp_path)
    tmp_sorted_ffn = tempfile.NamedTemporaryFile(delete=False, prefix='panphlan_', suffix='.sorted.ffn', dir=tmp_path)

    try:
        with tmp_ffn:
            # cat genefiles.ffn > merged_file.ffn
            cat_cmd = ['cat'] # same as for bowtie2, but based on gene.ffn files
            for f in path_gene_ffn_files:
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
def usearch_clustering(path_gene_ffn_files, identity_threshold_perc, clade, output_path, tmp_path, KEEP_UC, gene2description, TIME, VERBOSE):
    '''
    Note: If KEEP_UC, then <clusters>.uc is a file written in the output directory.
    Otherwise, <clusters>.uc is a temp file (in /tmp), deleted at the end of the computation
    '''
    # Merge and sort genes by length
    TIME, tmp_sorted_ffn = usearch_sortbylength(path_gene_ffn_files, tmp_path, TIME, VERBOSE)
    # usearch7 clustering
    tmp_uc, TIME = run_usearch(tmp_sorted_ffn.name, identity_threshold_perc / 100.0, clade, output_path, tmp_path, KEEP_UC, TIME, VERBOSE)
    # Convert usearch7 result
    merged_txt = output_path + 'usearch7_' + clade + '_genefamily_cluster.txt'
    TIME = convert_usearch_result(tmp_uc.name, merged_txt, TIME, VERBOSE)
    # get dictionary gene2family
    gene2family = usearch_get_gene2family_dict(merged_txt, VERBOSE)
    # Add prefix clade:genefamID: to geneIDs in centroid.ffn sequence file
    #  to do: add also function (gene description)
    family2centroidGeneID = usearch_centroids_add_geneID_prefix(clade, gene2family, output_path, gene2description)
    # clean up tmp files
    if not KEEP_UC:
        os.unlink(tmp_uc.name)
    if VERBOSE: print('[I] Remove usearch7 tmp results')
    os.remove(merged_txt)

    return gene2family, family2centroidGeneID, TIME
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
        print('\n[E] Missing "usearch7". Please install.\n')
        print('    usearch7 is required to merge gene sequences into gene-family cluster\n')
        print('    1) Download version 7 from http://drive5.com/usearch/')
        print('    2) Rename into "usearch7" and place it in any tools directory present in PATH')
        print('       For example, run the following commands:\n')
        print('     mv usearch7.0.1090_i86osx32 usearch7 # rename')
        print('     chmod a+x usearch7                   # set executable permissions')
        print('     sudo mv usearch7 /usr/local/bin/     # move to system tools directory')
        print('\n')
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
def add_filename_to_seqIDs(path_gene_fxx_files, tmp_path, Fxx, VERBOSE):
    '''
    To get unique geneIDs (or contigIDs) across all genomes:
    1) Copy gene ffn (or genome fna) files to TMP
    2) add filename as prefix to geneIDs "Filename:originalGeneID", or to contigIDs "Filename:originalContigID"
       add prefix only if not exist already (prefix does not need to follow by ':', "Filename_contig01" counts also as prefix exist)

    requires Biopython (Bio module)
    '''

    if VERBOSE and Fxx=='ffn': print('[I] To get unique geneIDs   across genomes: add filename as prefix to geneIDs')
    if VERBOSE and Fxx=='ffa': print('[I] To get unique contigIDs across genomes: add filename as prefix to contigIDs')
    # create new folder 'fxx_uniqueGeneIDs' in TMP
    new_fxx_folder = os.path.join(tmp_path, Fxx+'_uniqueSeqIDs','') # '' to get ending '/'
    if os.path.exists(new_fxx_folder):
        # print('\n\n ERROR: directory exist already: ' + new_fxx_folder)
        sys.exit('\n\n ERROR: directory exist already: ' + new_fxx_folder + '\n Please rename or remove: ' + tmp_path + '\n\n')
    else:
        os.makedirs(new_fxx_folder)

    # create new list of fxx files: new_path_gene_fxx_files
    new_path_gene_fxx_files = [os.path.join(new_fxx_folder,os.path.basename(f)) for f in path_gene_fxx_files]

    for (fxx_in,fxx_out) in zip(path_gene_fxx_files,new_path_gene_fxx_files):
        filename = os.path.splitext(os.path.basename(fxx_in))[0]
        with open(fxx_out, 'w') as f_out:
            for seq in SeqIO.parse(open(fxx_in), 'fasta'):
                if seq.id == seq.name:
                    seq.name=''
                if seq.id == seq.description.split()[0]:
                    seq.description=' '.join(seq.description.split()[1:])
                if not seq.id.startswith(filename):
                    seq.id = filename + ':' + seq.id
                r = SeqIO.write(seq, f_out, 'fasta')
                if r!=1:
                    sys.exit('[E] Error while writing sequence to '+Fxx+'-file:\n    ' + fxx_out)

    return new_path_gene_fxx_files
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
    path_genome_fna_files = sorted([os.path.join(fna_folder,f) for f in genomefiles])
    path_gene_ffn_files   = sorted([os.path.join(ffn_folder,f) for f in genefiles])

    # check usearch7 cluster depends on file order? (bvulgatus14) Yes, result depends on gene-file order
    # path_genome_fna_files = [ path_genome_fna_files[i] for i in [2,1,0,3]]  # Original file order
    # path_gene_ffn_files   = [ path_gene_ffn_files[i]   for i in [2,1,0,3]]

    print('\nExpected runtime: ' + str(len(genomefiles)*20) + ' minutes (start time: ' + time.strftime("%b %d %Y %H:%M") + ')\n')

    return path_genome_fna_files, path_gene_ffn_files
# ------------------------------------------------------------------------------
def clean_up(path_gene_ffn_files, tmp_path, VERBOSE):
    '''
    Remove files not needed anymore
    1) copy of sequence ffn and fna files having prefix to geneIDs "Filename:originalSeqID"
    2) TMP folder
    '''
    if VERBOSE: print('[I] Remove copies of sequence files (ffn or fna) from ' + tmp_path)
    for f in path_gene_ffn_files:
        if tmp_path in f: # make sure we deleting in the TMP directory
            os.remove(f)

    # remove empty tmp directories, if exist
    for rmdir in [os.path.join(tmp_path,'ffn_uniqueSeqIDs'), os.path.join(tmp_path,'fna_uniqueSeqIDs'), tmp_path]:
        if os.path.isdir(rmdir) and os.listdir(rmdir) == []:
            if rmdir == tmp_path:
                if VERBOSE: print('[I] Remove '+ tmp_path +' directory')
            os.rmdir(rmdir)
    # os.rmdir(os.path.join(tmp_path,'ffn_uniqueGeneIDs'))
    # os.rmdir(tmp_path)
    return True
# ------------------------------------------------------------------------------
def check_args():
    '''
    Check input arguments
    '''
    parser = PanPhlAnGenParser()
    args = vars(parser.parse_args())
    VERBOSE = args['verbose']
    if not VERBOSE: print('\nUse option --verbose to display progress information.\n')

    if VERBOSE:
        print('\nPanPhlAn pangenome generation version '+__version__)
        print('Python version: ' + sys.version.split()[0])
        print('System: ' + sys.platform)
        print(' '.join(sys.argv))
        print('')

    # valid options:
    # panphlan_pangenome_generation.py --i_gff  ncbi_download_gff/--fna ncbi_download_fna/ # add genome fna to gff (no bowtie, no usearch)
    # panphlan_pangenome_generation.py --i_gff  prokka_gff/   # only convert gff to fna and ffn (filename prefix added to seqID) (no bowtie, no usearch) (genome seq needs to be included in gff)
    # panphlan_pangenome_generation.py    ....   -c species # generate database

    # not valid option combinations
    if args['i_ffn'] and not args['i_fna']: sys.exit('\n Error: Genome sequence .fna files required, add option:  --i_fna genomes/\n')
    if args['i_ffn'] and args['i_fna'] and args['i_gff']: sys.exit('\n Error: Too many input options, use only (--fna & --ffn), (--gff & --fna), or --gff alone.\n')
    if not args['i_ffn'] and not args['i_fna'] and not args['i_gff']: sys.exit('\n Error: Please provide input files, either (--fna & --ffn), (--gff & --fna), or --gff alone.\n')

    if args['i_ffn'] and not args['clade']: sys.exit('\n Error: Species pangenome database name required, for example add option:  --clade ecoli17\n')

    # Roary
    if args['roary_dir'] and args['i_ffn']: sys.exit('\n Error: Please use --gff input (not --ffn). Required are the same gff files as used for Roary\n')

    if args['roary_dir']:
        roary_dir = args['roary_dir']
        roary_dir = os.path.abspath(roary_dir)
        if not os.path.exists(roary_dir):
            print('\n ' + roary_dir)
            sys.exit('\n Error (--roary_dir): Roary directory does not exist. Please check the above path.\n')
        if not os.path.isdir(roary_dir):
            print('\n ' + roary_dir)
            sys.exit('\n Error (--roary_dir): Please provide the directory of Roary output (not individual files).\n')
        roary_dir = os.path.join(roary_dir,'')
        args['roary_dir']=roary_dir

    # Check: GFF_FOLDER --------------------------------------------------------
    if args['i_gff']:
        ipath = args['i_gff']
        if not os.path.exists(ipath):
            show_error_message('Input folder -i_gff does not exist.')
            sys.exit(INEXISTENCE_ERROR_CODE)
        ipath = os.path.abspath(ipath)
        ipath = os.path.join(ipath,'')
        args['i_gff'] = ipath
        if VERBOSE: print('[I] Input gene GFF folder: ' + args['i_gff'])
    # Check: FFN_FOLDER --------------------------------------------------------
    if args['i_ffn']:
        ipath = args['i_ffn']
        if not os.path.exists(ipath):
            show_error_message('Input folder -i_ffn does not exist.')
            sys.exit(INEXISTENCE_ERROR_CODE)
        ipath = os.path.abspath(ipath)
        ipath = os.path.join(ipath,'')
        args['i_ffn'] = ipath
        if VERBOSE: print('[I] Input gene FFN folder: ' + args['i_ffn'])
    # Check: FNA_FOLDER --------------------------------------------------------
    if  args['i_fna']:
        ipath = args['i_fna']
        if not os.path.exists(ipath):
            show_error_message('Input folder -i_fna does not exist.')
            sys.exit(INEXISTENCE_ERROR_CODE)
        ipath = os.path.abspath(ipath)
        ipath = os.path.join(ipath,'')
        args['i_fna'] = ipath
        if VERBOSE: print('[I] Input genome FNA folder: ' + args['i_fna'])


    if args['clade']: # don't need to check if only converting gff to fna & ffn
        # Check: CLADE -------------------------------------------------------------
        args['clade']=args['clade'].replace('panphlan_','') # remove panphlan_ prefix (added later only for bowtie2)
        args['clade']=args['clade'].replace('_','-') # convert underscore '_' to dash '-'  (underscore is used as separator in _map)
        if VERBOSE: print('[I] Species database name: ' + args['clade'])

        # Check: IDENTITY_PERCENATGE -----------------------------------------------
        if not args['roary_dir']: # threshold not used for Roary inport
            identity_threshold_perc = args['th']
            if identity_threshold_perc < 0.0 or identity_threshold_perc > 100.0:
                args['th'] = 95.0
                if VERBOSE: print('[I] Invalid value for identity threshold percentage. Default value (95.0 %) has been set.')
            else:
                if VERBOSE: print('[I] Identity threshold percentage: ' + str(args['th']) + ' %.')

        # Check: OUTPUT_FOLDER -----------------------------------------------------
        opath = os.path.join(args['output'],'')
        if not os.path.exists(os.path.dirname(opath)):
            os.makedirs(opath)
        args['output'] = opath
        if VERBOSE: print('[I] Output folder: ' + args['output'])

        # Check: TEMP_FOLDER ------------------------------------------------------
        tmp_path = os.path.join(args['tmp'],'')
        if not os.path.exists(os.path.dirname(tmp_path)):
            os.makedirs(tmp_path)
        args['tmp'] = tmp_path
        if VERBOSE: print('[I] Temporary folder: ' + args['tmp'])

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
    if VERBOSE: print('\nSTEP 1. Checking required software installations ...')
    if args['clade']:
        bowtie2   = check_bowtie2(VERBOSE, PLATFORM)  # for generating .bt2 index files
    if args['clade'] and not args['roary_dir']:
        usearch7  = check_usearch7(VERBOSE, PLATFORM) # for generating usearch7 gene-family cluster


    if VERBOSE: print('\nSTEP 2. Prepare input gene and genome files ...')

    if args['i_gff'] and args['i_fna']: # if not 'clade': skip bowtie2 & usearch, only add genome to gff files
        print('[I] Add genome fna sequences to gff files')
        args['i_gff'] = gff_add_genome_seq(args['i_gff'], args['i_fna'], VERBOSE)

    if (args['i_gff'] and not args['i_fna'])  or  (args['i_gff'] and args['i_fna'] and args['clade']):
        # extract fna & ffn from gff always if 'clade'
        # if not 'clade': only extract fna & ffn from gff, but not if 'fna' added right before
        print('[I] Read gff files, extract gene location and annotation, write genome fna and gene ffn sequences.')
        gene2loc,gene2genome,gene2description,gene2gffdata,path_genome_fna_files,path_gene_ffn_files = read_gff_write_fna_ffn(args['i_gff'], VERBOSE)

    if args['i_fna'] and not args['i_gff'] and args['clade']: # classic input: fna & ffn files (usearch approach)
        path_genome_fna_files, path_gene_ffn_files = check_genomes(args['i_ffn'], args['i_fna'], VERBOSE)
        path_gene_ffn_files   = add_filename_to_seqIDs(path_gene_ffn_files,   args['tmp'], 'ffn', VERBOSE)
        path_genome_fna_files = add_filename_to_seqIDs(path_genome_fna_files, args['tmp'], 'fna', VERBOSE)
        gene2genome, gene2description = gene2genome_mapping(path_gene_ffn_files, VERBOSE)
        gene2loc = get_gene_locations(path_genome_fna_files, path_gene_ffn_files, VERBOSE)

    # sys.exit('\n --- Stop for testing/debugging --- \n\n') # for testing

    # clustering, write pangenome, bowtie2 (otherwise only convert gff to ffn fna)
    if args['clade']:
        # check if contigIDs from geneIDs are present in contigs from fna (bowtie2 index based on fna equals pangenome contig-IDs?)
        contig2genome = get_contigs(path_genome_fna_files, VERBOSE)
        check_for_valid_contigIDs(gene2loc,contig2genome)

        if VERBOSE: print('\nSTEP 3. Get gene-family clusters (usearch7 or Roary) ...')
        if args['roary_dir']: # later:  if roary_dir or run_roary
            roary_gene2family, roary_family2annotation, roary_genomeIDs = read_roary_gene_clustering(args['roary_dir'], VERBOSE)
            gene2family, locustag2gene = convert_roary_geneIDs(roary_gene2family,gene2gffdata, VERBOSE)
            family2centroidGeneID      = read_roary_centroids(args['roary_dir'], args['clade'], args['output'], locustag2gene, roary_family2annotation, VERBOSE) # copy pan_genome_reference.fa and get centroidIDs
            path_genome_fna_files      = reject_genomes_not_in_roary(path_genome_fna_files, roary_genomeIDs) # to run bowtie2 only with genomes present in roary clustering
        else: # Run usearch7 to get gene families cluster
            gene2family, family2centroidGeneID, TIME = usearch_clustering(path_gene_ffn_files,args['th'],args['clade'],
                                           args['output'],args['tmp'],KEEP_UC,gene2description,TIME,VERBOSE)

        # Write the pangenome database file: panphlan_clade_pangenome.csv
        if VERBOSE: print('\nSTEP 4. Write pangenome file ...')
        write_pangenome_file(gene2loc, gene2family, gene2genome, args['output'], args['clade'], VERBOSE)
        # Write gene-family annotation file
        if args['i_gff']:
            write_annotations_gff(args['clade'], args['output'], family2centroidGeneID, gene2gffdata, VERBOSE)



        # Get bowtie2 index files
        if VERBOSE: print('\nSTEP 5. Get bowtie2 index database ...')
        TIME = create_bt2_indexes(path_genome_fna_files, args['clade'], args['output'], args['tmp'], TIME, VERBOSE)


        clean_up(path_genome_fna_files, args['tmp'], VERBOSE)
        clean_up(path_gene_ffn_files,   args['tmp'], VERBOSE)
        end_program(time.time() - TOTAL_TIME)

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
