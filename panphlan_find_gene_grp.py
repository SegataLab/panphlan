#!/usr/bin/env python

"""
Analyzing PanPhlAn presence absence matrices
    - plot the matrix as heatmap using seaborn
    - based on location on the contigs, assess if a group of gene can be mobile element
"""

import os, subprocess, sys, time, bz2
import re
import random
from sklearn.cluster import OPTICS
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist, jaccard
from scipy import stats
import argparse as ap

author__ = 'Leonard Dubois and Nicola Segata (contact on https://forum.biobakery.org/)'
__version__ = '3.0'
__date__ = '20 April 2020'

OPTICS_MIN_PTS = 5

# ------------------------------------------------------------------------------
"""
Reads and parses the command line arguments of the script.
:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-i','--i_matrix', type=str, default = None,
                    help='Path to presence/absence matrix')
    p.add_argument('-o', '--output', type = str, default = None,
                    help='Path to ouput file with genes groups')
    p.add_argument('-c', '--cut_top', type = float, default = 0.01,
                    help='Amount of rows to remove based on top and bottom row sums. Default 0.01')
    p.add_argument('--out_plot', type = str, default = None,
                    help='Path to heatmap plot output.')
    p.add_argument('-p', '--pangenome', type = str, default = None,
                    help='Path to pangenome file.')
                    
    p.add_argument('--optics_xi', type = float, default = 0.01,
                    help='Xi parameter for OPTICS clustering')
    p.add_argument('-a', '--analysis', action = "store_true",
                    help='Should gene groups be analyzed?')
    p.add_argument('--empirical', type = int, default = 500,
                    help='How many ramdom sample in empirical pvalue generation ? Default 500')
    
    p.add_argument('-v', '--verbose', action='store_true',
                    help='Show progress information')
    return p.parse_args()

# ------------------------------------------------------------------------------
#   READ AND PROCESS PANPHLAN MATRIX
# ------------------------------------------------------------------------------

def read_and_filter_matrix(filepath, threshold_sums, verbose):    
    panphlan_matrix = pd.read_csv(filepath, sep = '\t', header = 0, index_col = 0)
    if verbose:
        print(' [I] Reading PanPhlAn presence/absence matrix from : ' + str(filepath))
        print('     Matrix with ' + str(panphlan_matrix.shape[0]) + ' genes families (rows) and')
        print('                 ' + str(panphlan_matrix.shape[1]) + ' samples (columns)')
    row_sums = panphlan_matrix.sum(axis = 1)
    thres_bottom = round(panphlan_matrix.shape[1] * threshold_sums)
    thres_top = panphlan_matrix.shape[1] - thres_bottom
    keep = (row_sums > thres_bottom) & (row_sums < thres_top)
    panphlan_matrix = panphlan_matrix[keep]
    if verbose:
        print(' [I] After removing genes too often present/absent. Cut-off at ' + str(threshold_sums * 100) + ' % of the samples :')
        print('     Matrix with ' + str(panphlan_matrix.shape[0]) + ' genes families (rows) and')
        print('                 ' + str(panphlan_matrix.shape[1]) + ' samples (columns)')
        
    return panphlan_matrix


def compute_dist(panphlan_matrix, verbose):
    if verbose:
        print(' [I] Comuting Jaccard distance between gene families...')
    a = pdist(panphlan_matrix, 'jaccard')
    dist_matrix = pd.DataFrame(squareform(a),
                                index = panphlan_matrix.index,
                                columns = panphlan_matrix.index)
    if verbose: print(' Done')
    return dist_matrix

# ------------------------------------------------------------------------------
#   CLUSTERING AND DEFINING GROUPS
# ------------------------------------------------------------------------------

def process_OPTICS(dist_matrix, xi_value, verbose):
    if verbose:
        print(' [I] Performing OPTICS clustering...')
    clustering = OPTICS(min_samples = OPTICS_MIN_PTS, 
                        cluster_method = "xi", xi = xi_value,
                        metric="precomputed").fit(dist_matrix)
    
    a = pd.DataFrame(clustering.labels_)
    # a[0].value_counts()
    optics_res = dict(zip(dist_matrix.index, clustering.labels_) )
    # dbscan_res = {"UniRef90_XXX" : cluster_ID, "UniRef90_YYY" : cluster_ID , ...}
    if verbose:
        print(' Done')
    return optics_res


def write_clusters(dbscan_res, out_file, operon_pval = None, subspec_pval = None):
    OUT = open(out_file, mode='w')
    #OUT.write("clust_ID\tsize\toperon_pval\tsubspec_pval\tUniRef_ID\n")
    OUT.write("clust_ID\tsize\toperon_pval\tUniRef_ID\n")
    for i in range(-1, max(dbscan_res.values()) + 1 ):
        OUT.write(str(i) + "\t")
        ids = [k for k,v in dbscan_res.items() if v == i]
        OUT.write(str(len(ids)) + "\t")
        if operon_pval:
            OUT.write(str(operon_pval[i]) + "\t")
        else:
            OUT.write("NA\t")
        # if subspec_pval:
        #     OUT.write(str(subspec_pval[i]) + "\t")
        # else:
        #     OUT.write("NA\t")
        OUT.write(";".join(ids))
        OUT.write("\n")
    OUT.close()

# ------------------------------------------------------------------------------
#   PLOTTING HEATMAP WITH LABELS
# ------------------------------------------------------------------------------

def plot_heatmap(panphlan_matrix, dbscan_res, out_path ):
    from matplotlib import pyplot as plt
    import seaborn as sns
    sns.set(color_codes = True, font_scale = 0.1)

    clust_names = set(dbscan_res.values())
    clust_to_col = dict(zip(clust_names, sns.color_palette("hls", len(clust_names) )))
    # # identify cluster of core genome : biggest clust that is not -1 (noise points)
    # table_count = {a : list(dbscan_res.values()).count(a) for a in dbscan_res.values()}
    # table_count = sorted(table_count, key = table_count.get, reverse = True)
    # if table_count[0] == -1 :
    #     core_gene_clust = table_count[1]
    # else:
    #     core_gene_clust = table_count[0]
    # # changing colors for core genes and rare genes 
    #clust_to_col.update({core_gene_clust : (0.2,0.2,0.2)})
    clust_to_col.update({-1 : (1,1,1)})

    lbl_to_col = {lbl : clust_to_col[clust] for lbl, clust in dbscan_res.items()  }
    lbl_to_col = pd.Series(lbl_to_col)

    plot_width = round(panphlan_matrix.shape[0] / 50)
    plot_height = round(panphlan_matrix.shape[1] / 50) 
    
    plt.figure()
    p = sns.clustermap(data = panphlan_matrix,
                    metric = 'jaccard',
                    row_colors = lbl_to_col,
                    cmap = sns.color_palette(["#f7f7f7", "#ca0020"])
                    )
    p.cax.set_visible(False)
    plt.savefig(out_path, dpi = 500)

# ------------------------------------------------------------------------------
#   ANALYZING GENES GROUPS
# ------------------------------------------------------------------------------

def get_span_ratio(gene_list, df):
    """Compute the ratio beetween sum of gene length and actual span of gene set on genome
    In case of multiple genome and multiple values return the highest one
    """
    df = df[df['UniRef'].isin(gene_list)]
    df = df.assign(length = lambda df : df["stop"] - df["start"])
    sum_length = df["length"].sum()
    tot_span = df["stop"].max() - df["start"].min()
    return tot_span / sum_length


def assessment_operon(clust_res, args):
    """For each group of genes detected by dbscan clustering, compute its
    spanning ratio (function above) and randomized spanning ratio of groups
    with the same size. An empirical pvalue is thus computed.
    """
    if args.verbose:
        print(' [I] Computing empirical pvalue for span of gene families clusters...')
        print('     Generating ' + str(args.empirical) + ' empirical samples to compute the pvalue (arg --empirical, default 500)')
        
    clust_is_operon = dict()
    table_count = {a : list(clust_res.values()).count(a) for a in clust_res.values()}
    table_count = sorted(table_count, key = table_count.get, reverse = True)

    pangenome_df = pd.read_csv(args.pangenome, sep = '\t', header = None)
    pangenome_df.columns = ["UniRef", "name", "genome", "contig", "start", "stop"]

    for cluster in table_count:
        if args.verbose : print("Analysing cluster : {} ".format(cluster))
        if cluster == table_count[0]:
            clust_is_operon[cluster] = "NA"
            continue
        genes = [k for k,v in clust_res.items() if v == cluster]
        cluster_info = pangenome_df[pangenome_df['UniRef'].isin(genes)]
        contigs = cluster_info.contig.unique()
        if args.verbose : print("Cluster span across {} contigs ".format(len(contigs)))
        result_contig = []
        for c in contigs:
            if args.verbose : print("Analysing contig : {} ".format(c))
            contig_info = cluster_info[cluster_info['contig'] == c]
            gene_of_contig = contig_info['UniRef']
            if len(gene_of_contig) < OPTICS_MIN_PTS:
                if args.verbose : print("Contig discarded. Not enough genes on it")
                result_contig.append(None)
                continue
            # assess genes span ratio
            cluster_span_ratio = get_span_ratio(gene_of_contig, contig_info)
            random_ratio = list()
            for i in range(args.empirical):
                contig_all_genes_info = pangenome_df[pangenome_df['contig'] == c]
                random_set = random.sample(range(0, contig_all_genes_info.shape[0]), k=len(gene_of_contig))
                random_set = list(contig_all_genes_info.iloc[random_set]["UniRef"])
                random_ratio.append(get_span_ratio(random_set, contig_info))
            # create sample of span values
            pvalue = sum([a > cluster_span_ratio for a in random_ratio]) / len(random_ratio)
            result_contig.append(pvalue)
            
        result_contig = [x for x in result_contig if not x is None]    
        
        if len(result_contig) >= 1:
            clust_is_operon.update({cluster : min(result_contig)})
        else:
            clust_is_operon.update({cluster : "NA"})
    
    return clust_is_operon


def assessment_subspecies_gene(dbscan_res, panphlan_matrix):
    """WIP 
    Check if group of genes is driving the phylogenetic dendrogramm
    WIP
    """
    clust_subspec = dict()
    table_count = {a : list(dbscan_res.values()).count(a) for a in dbscan_res.values()}
    table_count = sorted(table_count, key = table_count.get, reverse = True)
    for cluster in table_count:
        if cluster in table_count[0:2]:
            clust_subspec[cluster] = "NA"
            continue
        genes = [k for k,v in dbscan_res.items() if v == cluster]
        slice = panphlan_matrix.loc[genes]
        colsums = slice.sum()
        # DO NOT WORK
        present_grp = [a == len(genes) for a in colsums]
        absent_grp = [not a for a in present_grp]
        #present_grp = slice[slice.columns[present_grp]]
        #absent_grp = slice[slice.columns[absent_grp]]
        present_grp = panphlan_matrix[slice.columns[present_grp]]
        absent_grp = panphlan_matrix[slice.columns[absent_grp]]
        print(present_grp)

        pres_dist = pdist(present_grp.T, 'jaccard')
        absent_dist = pdist(absent_grp.T, 'jaccard')
        x= stats.mannwhitneyu(pres_dist, absent_dist)
        clust_subspec[cluster] = x.pvalue
        # pvalue always 0  does not work :(

    return clust_subspec

# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------
def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python 3, please update Python')
    args = read_params()

    panphlan_matrix = read_and_filter_matrix(args.i_matrix, args.cut_top, args.verbose)
    dist_matrix = compute_dist(panphlan_matrix, args.verbose)
    
    optics_res = process_OPTICS(dist_matrix, args.optics_xi, args.verbose)

    if args.analysis and args.output:

        operon_pval = assessment_operon(optics_res, args)
        #subspec_pval = assessment_subspecies_gene(dbscan_res, panphlan_matrix)
        #write_clusters(dbscan_res, args.output, operon_pval = operon_pval, subspec_pval = subspec_pval)
        write_clusters(optics_res, args.output, operon_pval = operon_pval)
    elif args.output:
        write_clusters(optics_res, args.output)

    if args.out_plot:
        plot_heatmap(panphlan_matrix, optics_res, args.out_plot)


if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
