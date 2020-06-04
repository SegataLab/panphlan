import setuptools
from setuptools.command.install import install
from io import open
import os

install_requires = ["numpy", "pandas", "scipy", "sklearn", "matplotlib", "seaborn"]
setuptools.setup(
    name='panphlan',
    version='3.0',
    author='Leonard Dubois',
    author_email='leonard.dubois@unitn.it',
    url='http://github.com/SegataLab/panphlan/',
    py_modules=['panphlan_map', 'panphlan_profiling', 'panphlan_download_pangenome', 'panphlan_find_gene_grp', 'misc'],
    entry_points = { "console_scripts" : [ "panphlan_map = panphlan_map:main",
                                           "panphlan_profiling = panphlan_profiling:main",
                                           "panphlan_download_pangenome = panphlan_download_pangenome:main",
                                           "panphlan_find_gene_grp = panphlan_find_gene_grp:main"
                                         ] },
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    description='PanPhlAn is a strain-level metagenomic profiling tool for identifying the gene composition and *in-vivo* transcriptional activity of individual strains in metagenomic samples. PanPhlAn’s ability for strain-tracking and functional analysis of unknown pathogens makes it an efficient tool for culture-free infectious outbreak epidemiology and microbial population studies.',
    install_requires=install_requires
)