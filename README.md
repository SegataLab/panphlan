# README #

----

## General info

PanPhlAn is a complete framework for strain-level metagenomic profiling. PanPhlAn is written in Python and works with versions 2.6, 2.7 and 3.x.

PanPhlAn is composed of three different modules (executable programs), which would define a computational pipeline:

* `panphlan_pangenome_generation.py`, to create the pangenome of a bacterial species starting from the DNA sequences of many genomes;
* `panphlan_map.py`, to map each gene to its abundance in a given input metagenomic sample;
* `panphlan_profile.py`, to execute a metagenomic/-transcriptomic profiling of the species in many metagenomic samples.

You can find the instruction for downloading and installing PanPhlAn in the [Wiki](https://bitbucket.org/CibioCM/panphlan/wiki/Home) page of the project, but briefly what you need to do is the following:

* install Usearch 7 in your system
* install Bowtie2 in your system
* install Samtools in your system
* install Python v. 2.6, 2.7 or 3.x (and Biopython module) in your system
* download PanPhlAn from this repository and install it in your system

### Mailing list ###

For questions and to keep up to date with information you can join our mailing list: [panphlan-users@googlegroups.com](mailto:panphlan-users@googlegroups.com)

The list is managed by

* Nicola Segata, principal investigator
* Thomas Tolio, programmer
* Matthias Scholz, researcher

----

## PanPhlAn modules usages

Here follows the usage of the three PanPhlAn modules.

Pangenome generation:

```
./panphlan_pangenome_generation.py
    --i_ffn FFN_FOLDER/ --i_fna FNA_FOLDER/ --clade CLADE --output OUTPUT_FOLDER/
    [--tmp TEMP_FILES_FOLDER/] [--uc] [--th UCLUST_PERC_IDENTITY] [--verbose]
```

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


Learn to use PanPhlAn reading this useful [tutorial](https://bitbucket.org/CibioCM/panphlan/wiki/Tutorial) we wrote for the beginners!

Or rather [know more](https://bitbucket.org/CibioCM/panphlan/wiki/Home) about PanPhlAn pipeline reading this technical guide.

----

[Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)