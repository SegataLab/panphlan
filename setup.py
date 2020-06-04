import setuptools
from setuptools.command.install import install
from io import open
import os

install_requires = ["biopython"]
setuptools.setup(
    name='panphlan',
    version='1.3',
    author='Matthias Scholz',
    author_email='matthias.scholz.de@gmail.com',
    url='http://github.com/SegataLab/panphlan/',
    py_modules=['panphlan_map', 'panphlan_profile', 'panphlan_new_pangenome_generation', 'panphlan_pangenome_generation', 'utils'],
    entry_points = { "console_scripts" : [ "panphlan_map = panphlan_map:main",
                                           "panphlan_profile = panphlan_profile:main",
                                           "panphlan_new_pangenome_generation = panphlan_new_pangenome_generation:main",
                                           "panphlan_pangenome_generation = panphlan_pangenome_generation:main"
                                         ] },
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    description='PanPhlAn is a strain-level metagenomic profiling tool for identifying the gene composition and *in-vivo* transcriptional activity of individual strains in metagenomic samples. PanPhlAn’s ability for strain-tracking and functional analysis of unknown pathogens makes it an efficient tool for culture-free infectious outbreak epidemiology and microbial population studies.',
    install_requires=install_requires
)
