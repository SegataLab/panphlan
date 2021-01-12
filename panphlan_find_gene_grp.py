#!/usr/bin/env python

"""
Analyzing PanPhlAn presence absence matrices
    - plot the matrix as heatmap using seaborn
    - based on location on the contigs, assess if a group of gene can be mobile element
"""

import os, subprocess, sys, time, bz2
import re
import random
import itertools
from sklearn.cluster import OPTICS
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist, jaccard
from scipy import stats
import argparse as ap
import multiprocessing
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests


author__ = 'Leonard Dubois and Nicola Segata (contact on https://forum.biobakery.org/)'
__version__ = '3.0.2'
__date__ = '4 Dec 2020'

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
    p.add_argument('--cut_core_thres', type = float, default = 0.9,
                    help='Remove gene families present in [cut_core_thres] or more of the sample. Defaul is 90 percent (extended core)')
    p.add_argument('--cut_cloud_thres', type = float, default = 0.1,
                    help='Remove gene families present in [cut_core_thres] or less of the sample. Defaul is 5 percent')
    p.add_argument('--out_plot', type = str, default = None,
                    help='Path to heatmap plot output.')
    p.add_argument('--plot_clust_pa', type = str, default = None,
                    help='Output : Path to heatmap of clusters presence-absence plot')                
    p.add_argument('-p', '--pangenome', type = str, default = None,
                    help='Path to pangenome file.')

    p.add_argument('--optics_xi', type = float, default = 0.05,
                    help='Xi parameter for OPTICS clustering')
    p.add_argument('--n_jobs', type = int, default = 4,
                    help='How many cores to use for OPTICS clustering.  Default 4')

    p.add_argument('--close_analysis', action='store_true',
                    help='Compute analysis of genes proximity in genomes')
    p.add_argument('--ssp_analysis', action='store_true',
                    help='Compute analysis of genes cluster links to intra species phylogeny')
    p.add_argument('--empirical', type = int, default = 1000,
                    help='How many ramdom sample in empirical pvalue generation ? Default 1000')

    p.add_argument('-v', '--verbose', action='store_true',
                    help='Show progress information')
    return p.parse_args()

# ------------------------------------------------------------------------------
#   READ AND PROCESS PANPHLAN MATRIX
# ------------------------------------------------------------------------------

def read_and_filter_matrix(filepath, core_thres, cloud_thres, verbose):
    panphlan_matrix = pd.read_csv(filepath, sep = '\t', header = 0, index_col = 0)
    if verbose:
        print(' [I] Reading PanPhlAn presence/absence matrix from : ' + str(filepath))
        print('     Matrix with ' + str(panphlan_matrix.shape[0]) + ' genes families (rows) and')
        print('                 ' + str(panphlan_matrix.shape[1]) + ' samples (columns)')

    row_sums = panphlan_matrix.sum(axis = 1)
    thres_bottom = round(panphlan_matrix.shape[1] * core_thres)
    keep = row_sums < thres_bottom
    panphlan_matrix = panphlan_matrix[keep]
    if verbose:
        print(' [I] After removing genes families too often present. Cut-off at ' + str(core_thres * 100) + ' % of the samples :')
        print('     Matrix with ' + str(panphlan_matrix.shape[0]) + ' genes families (rows) and')
        print('                 ' + str(panphlan_matrix.shape[1]) + ' samples (columns)')
        
    thres_top = round(panphlan_matrix.shape[1] * cloud_thres)
    keep = row_sums > thres_top
    panphlan_matrix = panphlan_matrix[keep]
    if verbose:
        print(' [I] After removing genes too often absent. Cut-off at ' + str(cloud_thres * 100) + ' % of the samples :')
        print('     Matrix with ' + str(panphlan_matrix.shape[0]) + ' genes families (rows) and')
        print('                 ' + str(panphlan_matrix.shape[1]) + ' samples (columns)')

    return panphlan_matrix


def compute_dist(panphlan_matrix, verbose):
    if verbose:
        print(' [I] Computing Jaccard distance between gene families...')
    a = pdist(panphlan_matrix, 'jaccard')
    dist_matrix = pd.DataFrame(squareform(a),
                                index = panphlan_matrix.index,
                                columns = panphlan_matrix.index)
    if verbose: print(' Done')
    return dist_matrix

# ------------------------------------------------------------------------------
#   CLUSTERING AND DEFINING GROUPS
# ------------------------------------------------------------------------------

def process_OPTICS(dist_matrix, xi_value, n_jobs, verbose):
    if verbose:
        print(' [I] Performing OPTICS clustering...')
    clustering = OPTICS(min_samples = OPTICS_MIN_PTS,
                        cluster_method = "xi", xi = xi_value,
                        metric="precomputed", n_jobs = n_jobs).fit(dist_matrix)

    optics_res_df = pd.DataFrame(clustering.labels_)
    # optics_res_df[0].value_counts()
    optics_res_dict = dict(zip(dist_matrix.index, clustering.labels_) )
    # optics_res = {"UniRef90_XXX" : cluster_ID, "UniRef90_YYY" : cluster_ID , ...}
    if verbose:
        print(' Done')
    return optics_res_dict, optics_res_df


def write_clusters(dbscan_res, out_file, operon_pval = None, subspec_pval = None):
    OUT = open(out_file, mode='w')
    OUT.write("clust_ID\tsize\toperon_pval\tsubspec_pval\tUniRef_ID\n")
    for i in range(0, max(dbscan_res.values()) + 1 ):
        OUT.write(str(i) + "\t")
        ids = [k for k,v in dbscan_res.items() if v == i]
        OUT.write(str(len(ids)) + "\t")
        if operon_pval:
            OUT.write(str(operon_pval[i]) + "\t")
        else:
            OUT.write("NA\t")
        if subspec_pval:
            OUT.write(str(subspec_pval["clust_" + str(i)]) + "\t")
        else:
            OUT.write("NA\t")
        OUT.write(";".join(ids))
        OUT.write("\n")
    OUT.close()

# ------------------------------------------------------------------------------
#   PLOTTING HEATMAP WITH LABELS
# ------------------------------------------------------------------------------

def plot_heatmap(panphlan_matrix, out_path, clust_res = None  ):
    import seaborn as sns
    from matplotlib import pyplot as plt


    sys.setrecursionlimit(10**7)

    if clust_res is None:
        plt.figure(figsize=(20,20))
        p = sns.clustermap(data = panphlan_matrix, metric = 'jaccard',
                        cmap = sns.color_palette(["#f7f7f7", "#ca0020"]),
                        xticklabels = [], yticklabels = [],
                        dendrogram_ratio=0.1)
        p.cax.set_visible(False)
        plt.axis('off')
        p.savefig(out_path, dpi = 300)
    else :

        sns.set(color_codes = True, font_scale = 0.1)

        clust_names = set(clust_res.values())
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

        lbl_to_col = {lbl : clust_to_col[clust] for lbl, clust in clust_res.items()  }
        lbl_to_col = pd.Series(lbl_to_col)

        #plot_width = round(panphlan_matrix.shape[0] / 50)
        #plot_height = round(panphlan_matrix.shape[1] / 50)
        #plt.figure()
        p = sns.clustermap(data = panphlan_matrix,
                        metric = 'jaccard',
                        row_colors = lbl_to_col,
                        cmap = sns.color_palette(["#f7f7f7", "#ca0020"]),
                        xticklabels = [], yticklabels = [],
                        dendrogram_ratio=0.1)
        p.cax.set_visible(False)
        p.savefig(out_path, dpi = 300)


def plot_hm_genes_clusters(clust_PA_matrix, out_path, ssp_pval, args  ):
    import seaborn as sns
    from matplotlib import pyplot as plt

    clust_PA_matrix = clust_PA_matrix.astype(int)
    print(' [I] Exporting HM of clusters')
    sys.setrecursionlimit(10**7)

    sns.set(color_codes = True, font_scale = 0.1)

    clust_names = ["clust_" + str(x) for x in set(ssp_pval.keys())]
    annot_col_map = {True : "#f7f7f7", False : "#1a9641"}
    lbl_to_col = dict()
    for lbl, pval in ssp_pval.items():
        col = annot_col_map[pval <= 0.01]
        lbl_to_col.update({lbl : col})
    lbl_to_col = pd.Series(lbl_to_col)

    #plot_width = round(panphlan_matrix.shape[0] / 50)
    #plot_height = round(panphlan_matrix.shape[1] / 50)
    #plt.figure()
    p = sns.clustermap(data = clust_PA_matrix,
                    metric = 'jaccard',
                    row_colors = lbl_to_col,
                    cmap = sns.color_palette(["#f7f7f7", "#2c7bb6"]),
                    xticklabels = [], yticklabels = [],
                    dendrogram_ratio=0.1)
    p.cax.set_visible(False)
    p.savefig(out_path, dpi = 300)
# ------------------------------------------------------------------------------
#   ANALYZING GENES GROUPS - Proximity
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
        print('     Generating ' + str(args.empirical) + ' empirical samples to compute the pvalue (arg --empirical, default 1000)')

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
    
# ------------------------------------------------------------------------------
#   ANALYZING GENES GROUPS - Phylogeny
# ------------------------------------------------------------------------------


def create_cluster_PA_matrix(clust_res, panphlan_matrix):
    """
    Collapse the presence-absence signature of genes families to genes clusters level
    Gene cluster is considered present if more than half of the genes are present
    """
    clust_res.index = panphlan_matrix.index
    clust_IDs = set(clust_res[0].values)
    clust_IDs.remove(-1)
    clust_PA_matrix = pd.DataFrame(index = ["clust_" + str(x) for x in clust_IDs], 
                                    columns = panphlan_matrix.columns )
    for clust in clust_IDs :
        uniref_in_clust = clust_res[clust_res[0] == clust].index
        pres_abs_submat = panphlan_matrix.loc[uniref_in_clust]
        clust_PA = pres_abs_submat.apply(lambda x : int(sum(x) > 0.5 * len(uniref_in_clust))).T
        clust_PA_matrix.loc["clust_" + str(clust)] = clust_PA
    return clust_PA_matrix


def assessment_cluster_phylo_link(cluster, clust_PA_matrix, samples_dist, args):
    if args.verbose:
        print(' [I] Analyzing genes OPTICS cluster {}'.format(cluster))
    grp_present = clust_PA_matrix.columns[clust_PA_matrix.loc[str(cluster)] == 1]
    dist_grp_present = samples_dist.loc[samples_dist['sampleA'].isin(grp_present) & samples_dist['sampleB'].isin(grp_present)]['dist']
    dist_inter_grp = samples_dist.loc[~samples_dist['sampleA'].isin(grp_present) ^ ~samples_dist['sampleB'].isin(grp_present)]['dist']
    #mean_base = dist_inter_grp.mean() - dist_grp_present.mean() 
    mean_base = np.std(dist_inter_grp) - np.std(dist_grp_present)
    mean_base = abs(mean_base)
    
    empirical_means = []
    for i in range(1, args.empirical):
        rand_samples = random.sample(set(clust_PA_matrix.columns), k = len(grp_present))
        dist_rand_samples = samples_dist.loc[samples_dist['sampleA'].isin(rand_samples) & samples_dist['sampleB'].isin(rand_samples)]['dist']
        dist_not_rand_samples = samples_dist.loc[samples_dist['sampleA'].isin(rand_samples) ^ samples_dist['sampleB'].isin(rand_samples)]['dist']
        # here to avoid computing the mean over way too many values and having a narrow range of empirical value, a subset is extracted
        dist_rand_samples = random.sample(list(dist_rand_samples), k = 300)
        dist_not_rand_samples = random.sample(list(dist_not_rand_samples), k = 300)
        # empirical_means.append(dist_rand_samples.mean())
        empirical_means.append(abs(np.mean(dist_not_rand_samples) - np.mean(dist_rand_samples)))
        
    clust_pval = sum([mean_base < x for x in empirical_means]) / args.empirical
    print("Cluster {} has pval {}".format(cluster, clust_pval))
    print("Mean Base {} ".format(mean_base))
    
    return([cluster, clust_pval])


# TEST WITH MEDOID
def assessment_cluster_phylo_link_medoid(cluster, clust_PA_matrix, sample_dist_matrix, args):
    if args.verbose:
        print(' [I] Analyzing genes OPTICS cluster {}'.format(cluster))
    grp_present = clust_PA_matrix.columns[clust_PA_matrix.loc[str(cluster)] == 1]
    
    tmp = sample_dist_matrix[grp_present]
    tmp = tmp.T[grp_present]
    medoid_index = np.argmin(tmp.sum(axis=0))
    mean_dist_to_medoid = tmp.iloc[medoid_index].mean()
    empirical_means = []
    for i in range(1, args.empirical):
        rand_samples = random.sample(set(clust_PA_matrix.columns), k = len(grp_present ))
        tmp = sample_dist_matrix[rand_samples]
        tmp = tmp.T[rand_samples]
        medoid_index_emp = np.argmin(tmp.sum(axis=0))
        empirical_means.append(tmp.iloc[medoid_index_emp].mean())
        
    clust_pval = sum([mean_dist_to_medoid > x for x in empirical_means]) / args.empirical
    print("Cluster {} has pval {}".format(cluster, clust_pval))
    print("Mean Base {} ".format(mean_dist_to_medoid))
    return([cluster, clust_pval])


def assessment_subspecies_gene(clust_res, panphlan_matrix, clust_PA_matrix, args):
    """WIP
    Check if group of genes is linked to the phylogenetic dendrogramm
    Check if distance between samples is significantly different between
    the samples which have the gene cluster (or not) and the others.
    Use empirical pvalue (random samples)
    WIP
    """
    clust_res.index = panphlan_matrix.index
    clust_IDs = set(clust_res[0].values)
    clust_IDs.remove(-1)
    
    # then ompute pairwise distances of samples based on PanPhlAn presence_absence
    # Shown to be correlated with StrainPhlAn distances
    samples_dist = pd.DataFrame(itertools.combinations(panphlan_matrix.columns, 2),
                                columns=['sampleA','sampleB'])
    samples_dist['dist'] = pdist(panphlan_matrix.T, 'jaccard')
    sample_dist_matrix = pd.DataFrame(squareform(samples_dist['dist'] ),
                                index = panphlan_matrix.columns,
                                columns = panphlan_matrix.columns)
    
    clust_phylo_test = dict()
    
    par_res = Parallel(n_jobs=args.n_jobs)(
                delayed(assessment_cluster_phylo_link)(i, clust_PA_matrix, samples_dist, args) for i in clust_PA_matrix.index)
    
    for x,y in par_res :
        clust_phylo_test.update({x : y})    
    
    # stick with bonferroni ? or BH ?
    clust_phylo_test_df = pd.DataFrame.from_dict(clust_phylo_test, orient='index', columns = ['raw_pval'])
    clust_phylo_test_df['pval_bonf'] = multipletests(clust_phylo_test_df['raw_pval'], method = 'bonferroni')[1]
    clust_phylo_test_bonf = dict()
    for index, row in clust_phylo_test_df.iterrows():
        clust_phylo_test_bonf.update({index : row['pval_bonf']})

    
    # return clust_phylo_test_bonf
    return clust_phylo_test

# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------
def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python 3, please update Python')
    args = read_params()

    panphlan_matrix = read_and_filter_matrix(args.i_matrix, args.cut_core_thres, args.cut_cloud_thres, args.verbose)

    if args.output:
        dist_matrix = compute_dist(panphlan_matrix, args.verbose)

        optics_res_dict, optics_res_df = process_OPTICS(dist_matrix, args.optics_xi, args.n_jobs, args.verbose)
        if args.close_analysis:
            operon_pval = assessment_operon(optics_res_dict, args)
        else:
            operon_pval = None
        if args.ssp_analysis:
            clust_PA_matrix = create_cluster_PA_matrix(optics_res_df, panphlan_matrix)
            ssp_pval = assessment_subspecies_gene(optics_res_df, panphlan_matrix, clust_PA_matrix, args)
            # with open('ssp_pval.pkl', 'wb') as out_pkl:
            #     pickle.dump(ssp_pval, out_pkl)
    
            # with open('ssp_pval.pkl', 'rb') as in_pkl:
            #     ssp_pval = pickle.load(in_pkl)     
            # del ssp_pval[-1]
            # with open('PA_clust.pkl', 'wb') as out_pkl:
            #     pickle.dump(clust_PA_matrix, out_pkl)
            
        else :
            ssp_pval = None
            
        # print(clust_PA_matrix)    
        write_clusters(optics_res_dict, args.output, operon_pval = operon_pval, subspec_pval = ssp_pval)
        if args.out_plot:
            plot_heatmap(panphlan_matrix, args.out_plot,  optics_res_dict)
        if args.plot_clust_pa:
            plot_hm_genes_clusters(clust_PA_matrix, args.plot_clust_pa, ssp_pval, args )
    elif args.out_plot:
        plot_heatmap(panphlan_matrix, args.out_plot)



if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
