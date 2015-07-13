# README #

----

## PanPhlAn - strain detection and characterization

#### Pangenome-based Phylogenomic Analysis 

PanPhlAn is a strain-level metagenomic profiling tool
for identifying the gene composition and *in-vivo* transcriptional activity of individual strains
in metagenomic samples. PanPhlAn’s ability for strain-tracking and functional analysis of unknown
pathogens makes it an efficient tool for culture-free infectious outbreak epidemiology and
microbial population studies.

PanPhlAn is written in Python and covers the three main tasks:

* `panphlan_pangenome_generation.py`, to create the pangenome database of a bacterial species [→wiki](https://bitbucket.org/CibioCM/panphlan/wiki/panphlan_pangenome_generation)
* `panphlan_map.py`, to profile each metagenomic sample by mapping it against the species specific database
* `panphlan_profile.py`, to merge and process the mapping results for getting the final gene presence/absence and transcriptional matrices



PanPhlAn runs under Ubuntu/Linux and requires the following software tools to be installed on your system:

* Bowtie2 
* Samtools
* Usearch 7 (only for generating your own species database)
* Python version 2.7 or 3.x (including the Biopython module) 

For more information, see our wiki [→ Download and Installation](https://bitbucket.org/CibioCM/panphlan/wiki/Download_and_Installation).

### User support ###

[→ User Tutorial](https://bitbucket.org/CibioCM/panphlan/wiki/Home)

For questions and to keep up to date with information you can join our [email-based group and discussion forum](https://groups.google.com/forum/#!forum/panphlan-users) 

Or, you can directly write your question to [panphlan-users@googlegroups.com](mailto:panphlan-users@googlegroups.com)

The user group is managed by

* [Matthias Scholz](http://www.matthias-scholz.de/), algorithm design
* [Thomas Tolio](https://www.linkedin.com/in/thomastolio), programmer
* [Nicola Segata](http://segatalab.cibio.unitn.it/), principal investigator

[PanPhlAn](http://segatalab.cibio.unitn.it/tools/panphlan/) is a project of the [Computational Metagenomics Lab at CIBIO](http://segatalab.cibio.unitn.it/), University of Trento, Italy



# Example of *E. coli* strain profiling #

## Characterization of the German 2011 *E. coli* outbreak strain ##
![PanPhlAn heatmap](http://segatalab.cibio.unitn.it/images/panphlan_heatmap_ecoli_outbreak.png "PanPhlAn heatmap")

PanPhlAn profiling of the German outbreak metagenomes using a reference database in which the target outbreak genome is missing. **(a)** Hierarchical clustering. The heatmap displays presence/absence gene-family profiles of 110 reference strains (bright colored columns) and of 12 metagenomically detected strains (darker columns). Most outbreak samples cluster together due to almost identical profiles (right), with four samples (left) showing different profiles due to the presence of additional dominant *E. coli* strains overlying the target outbreak strain. **(b)** Functional analysis of outbreak-specific gene-families (Fisher exact test) confirmed that the outbreak strain is a combination of a EAEC pathogen (pAA plasmid) with acquired Shiga toxin and antibiotic resistance genes, complemented with a set of enriched virulence-related functions and pathway modules. 

### Reference ###
Matthias Scholz\*, Doyle V. Ward\*, Thomas Tolio, Moreno Zolfo, Francesco Asnicar, Duy Tin Truong, Edoardo Pasolli, Adrian Tett, Ardythe L. Morrow, and Nicola Segata (\* Equal contribution)  
**Strain-level microbial epidemiology and population genomics from shotgun metagenomics**  
*Under review* 

----