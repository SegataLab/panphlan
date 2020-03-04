## PanPhlAn - strain detection and characterization

#### Pangenome-based Phylogenomic Analysis

PanPhlAn is a strain-level metagenomic profiling tool
for identifying the gene composition and *in-vivo* transcriptional activity of individual strains
in metagenomic samples. PanPhlAn’s ability for strain-tracking and functional analysis of unknown
pathogens makes it an efficient tool for culture-free infectious outbreak epidemiology and
microbial population studies.

PanPhlAn is written in Python and covers the three main tasks:

* `panphlan_pangenome_generation.py`, to create the pangenome database of a bacterial species
  [→ read more](https://bitbucket.org/CibioCM/panphlan/wiki/panphlan_pangenome_generation)
* `panphlan_map.py`, to profile each metagenomic sample by mapping it against the species specific database
  [→ read more](https://bitbucket.org/CibioCM/panphlan/wiki/panphlan_map)
* `panphlan_profile.py`, to merge and process the mapping results for getting the final gene presence/absence and transcriptional matrices
  [→ read more](https://bitbucket.org/CibioCM/panphlan/wiki/panphlan_profile)



PanPhlAn runs under Ubuntu/Linux and requires the following software tools to be installed on your system:

* Bowtie2
* Samtools
* Python 3.x (including the Biopython module)

For more information, see our wiki [→ Download and Installation](https://bitbucket.org/CibioCM/panphlan/wiki/Download_and_Installation).

## Contact & User support ##

[→ User Tutorial](https://bitbucket.org/CibioCM/panphlan/wiki)

The PanPhlAn software team: [Matthias Scholz](http://www.matthias-scholz.de/) (algorithm design), [Thomas Tolio](https://www.linkedin.com/in/thomastolio) (programmer), Leonard Dubois and [Nicola Segata](http://segatalab.cibio.unitn.it/) (principal investigator).  
For help, use the [bioBakery help forum](https://forum.biobakery.org/).

----

[PanPhlAn](http://segatalab.cibio.unitn.it/tools/panphlan/) is a project of the [Computational Metagenomics Lab at CIBIO](http://segatalab.cibio.unitn.it/), University of Trento, Italy
