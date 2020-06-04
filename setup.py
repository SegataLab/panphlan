import setuptools
from setuptools.command.install import install
from io import open
import os

install_requires = ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "seaborn"]
setuptools.setup(
    name='panphlan',
    version='3.0',
    author='Leonard Dubois',
    author_email='leonard.dubois@unitn.it',
    url='http://github.com/SegataLab/panphlan/',
    packages = setuptools.find_packages(),
    package_dir = {'panphlan' : '' },
    scripts=['panphlan_map.py', 'panphlan_profiling.py', 'panphlan_download_pangenome.py', 'panphlan_find_gene_grp.py', 'misc.py'],
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    description='PanPhlAn is a strain-level metagenomic profiling tool for identifying the gene composition and *in-vivo* transcriptional activity of individual strains in metagenomic samples. PanPhlAn’s ability for strain-tracking and functional analysis of unknown pathogens makes it an efficient tool for culture-free infectious outbreak epidemiology and microbial population studies.',
    install_requires=install_requires
)