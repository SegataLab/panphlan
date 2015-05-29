# README #

----

## PanPhlAn

#### Pangenome-based Phylogenomic Analysis 

PanPhlAn is a strain-level metagenomic profiling tool
for identifying the gene composition and *in-vivo* transcriptional activity of individual strains
in metagenomic samples. PanPhlAn’s ability for strain-tracking and functional analysis of unknown
pathogens makes it an efficient tool for culture-free infectious outbreak epidemiology and
microbial population studies.

PanPhlAn is written in Python and covers the three main tasks:

* `panphlan_pangenome_generation.py`, to create the pangenome database of a bacterial species
* `panphlan_map.py`, to profile each metagenomic sample by mapping it against the species specific database
* `panphlan_profile.py`, to merge and process the mapping results for getting the final gene presence/absence and transcriptional matrices

PanPhlAn runs under Ubuntu/Linux and requires the following software tools to be installed on your system:

* Bowtie2 
* Samtools
* Usearch 7 (only if you want to generate your own species database)
* Python version 2.6, 2.7 or 3.x (including the Biopython module) 

For more information, see our [wiki](https://bitbucket.org/CibioCM/panphlan/wiki).

### User support ###

For questions and to keep up to date with information you can join our [email-based group and discussion forum](https://groups.google.com/forum/#!forum/panphlan-users) 

Or, you can directly write your question to [panphlan-users@googlegroups.com](mailto:panphlan-users@googlegroups.com)

The user group is managed by

* Matthias Scholz, algorithm design
* Thomas Tolio, programmer
* Nicola Segata, principal investigator

[PanPhlAn](http://cibiocm.bitbucket.org/tools/panphlan.html) is a project of the [Computational Metagenomics Lab at CIBIO](http://cibiocm.bitbucket.org/), University of Trento, Italy

----

## PanPhlAn modules usages

Pangenome generation:

```
./panphlan_pangenome_generation.py
    --i_ffn FFN_FOLDER/ --i_fna FNA_FOLDER/ --clade CLADE --output OUTPUT_FOLDER/
    [--tmp TEMP_FILES_FOLDER/] [--uc] [--th UCLUST_PERC_IDENTITY] [--verbose]
```
[Pre-processed pangenome databases for download](https://bitbucket.org/CibioCM/panphlan/wiki/Pangenome%20databases) are provided for more than 400 species.


Gene abundances mapping:

```
./panphlan_map.py
    --input INPUT_FILE --clade CLADE [--output OUTPUT_FILE]
    [--th_mismatches NUMBER_OF_MISMATCHES] [--out_bam OUTPUT_BAM_FILE]
    [--nproc NUMBER_OF_PROCESSORS] [--mGB MEMORY-GIGABYTES]
    [--readLength READS_LENGTH] [--tmp TEMP_FILES_FOLDER/] [--verbose]
```

Metagenomic/-transcriptomic profiling:

```
./panphlan_profile.py
    --i_dna GENE_ABUNDANCES_FOLDER/ --clade CLADE
    [--i_rna TRANSCRIPT_ABUNDANCES_FOLDER/]
    [--sample_pairs DNA_RNA_SAMPLE_PAIRS_FILE]
    [--th_present PRESENCE_THRESHOLD [--th_zero ABSENCE_THRESHOLD
                                      --th_multicopy MULTICOPY_THRESHOLD]]
    [--min_coverage MINIMUM_MEDIAN_COVERAGE]
    [--left_max MAXIMUM_LEFT_THRESHOLD] [--right_min MINIMUM_RIGHT_THRESHOLD]
    [--rna_max_zeros MAXIMUM_ZERO_COVERED_GENEFAMILIES_PERCENTAGE]
    [--strain_similarity_perc STRAIN_SIMILARITY_PERCENTAGE]
    [--np NON_PRESENCE_TOKEN] [--nan NOT_A_NUMBER_TOKEN]
    [--o_dna GENEFAMILY_PRESENCE_MATRIX] [--o_cov GENEFAMILY_COVERAGE_MATRIX]
    [--o_idx GENEFAMILY_DNAINDEX_MATRIX] [--o_rna GENEFAMILY_RNA_EXPRESSION_MATRIX]
    [--o_covplot PLOT_FILE] [--o_covplot_normed NORMED_PLOT_FILE]
    [--add_strains] [--interactive] [--verbose]
```


To find out more about how to use PanPhlAn, read our [User Tutorial](https://bitbucket.org/CibioCM/panphlan/wiki/Tutorial).

Or, to learn more about the algorithm behind PanPhlAn, read our [Technical Guide](https://bitbucket.org/CibioCM/panphlan/wiki/Home).

----