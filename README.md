
# PanPhlAn 3 - strain detection and characterization

### Pangenome-based Phylogenomic Analysis

PanPhlAn is a strain-level metagenomic profiling tool for identifying
the gene composition of individual strains in metagenomic samples.
PanPhlAn’s ability for strain-tracking and functional analysis of unknown
pathogens makes it an efficient tool for culture-free microbial population studies.

PanPhlAn is written in Python and covers the 3 main tasks:

* `panphlan_download_pangenome.py`, to download pangenome files (fasta file of contigs, BowTie2 indexes and general information) for over 3,000 species
* `panphlan_map.py`, to profile each metagenomic sample by mapping it against the species of interest
* `panphlan_profile.py`, to merge and process the mapping results in order to get the final gene presence/absence matrix

For custom pangenome generation (advanced) see the [PanPhlAn exporter](https://github.com/SegataLab/PanPhlAn_pangenome_exporter)

---
PanPhlAn runs under Ubuntu/Linux and requires the following software tools to be installed on your system:

* Bowtie2 (bowtie2 >=2.3.0)
* Samtools (samtools >=1.9)
* Python 3 (python >=3.7)

And the following Python libraries:

* numpy
* pandas
* scipy

---

If you use PanPhlAn, please cite:

> [**Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3**](https://www.biorxiv.org/content/10.1101/2020.11.19.388223v1) *Francesco Beghini, Lauren J McIver, Aitor Blanco-Miguez, Leonard Dubois, Francesco Asnicar, Sagun Maharjan, Ana Mailyan, Andrew Maltez Thomas, Paolo Manghi, Mireia Valles-Colomer, George Weingart, Yancong Zhang, Moreno Zolfo, Curtis Huttenhower, Eric A Franzosa, Nicola Segata*. bioRxiv preprint (2020)

---

### Help

For any help see :
* The [PanPhlan wiki](https://github.com/SegataLab/panphlan/wiki/Home_3_0)
* The PanPhlAn step-by-step [tutorial](https://github.com/SegataLab/panphlan/wiki/Tutorial-3_0)
* The [bioBakery help forum](https://forum.biobakery.org/) for overall discussions on various metagenomics softwares. (Purely technical issues should better be raised on GitHub than on the forum.)




----

PanPhlAn is a project of the [Computational Metagenomics Lab at CIBIO](http://segatalab.cibio.unitn.it/), University of Trento, Italy
