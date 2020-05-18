# README #

## PanPhlAn 3 - strain detection and characterization 

#### Pangenome-based Phylogenomic Analysis

PanPhlAn is a strain-level metagenomic profiling tool for identifying
the gene composition of individual strains in metagenomic samples.
PanPhlAn’s ability for strain-tracking and functional analysis of unknown
pathogens makes it an efficient tool for culture-free microbial population studies.

PanPhlAn is written in Python and covers the 4 main tasks:

* `panphlan_download_pangenome.py`, to download pangenome files (fasta, BowTie2 indexes and general information) for over 3,000 species
* `panphlan_map.py`, to profile each metagenomic sample by mapping it against the species of interest
* `panphlan_profile.py`, to merge and process the mapping results in order to get the final gene presence/absence matrix
* `panphlan_find_gene_grp.py`, organise OPTICS clustering to find some group of gene with similar profile and assess if they could be mobile elements in the genome. Also plot the presence/absence matrix as Heatmap. 

PanPhlAn runs under Ubuntu/Linux and requires the following software tools to be installed on your system:

* Bowtie2
* Samtools
* Python 3

And the following Python libraries:

* numpy
* pandas
* scipy
* sklearn (only if using `panphlan_find_gene_grp.py`)  
If visualizations are made, one also needs :
* matplotlib
* seaborn

For any help see the wiki or the [bioBakery forum](https://forum.biobakery.org/)

----

[PanPhlAn] is a project of the [Computational Metagenomics Lab at CIBIO](http://segatalab.cibio.unitn.it/), University of Trento, Italy
