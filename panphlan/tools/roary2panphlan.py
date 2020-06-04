#!/usr/bin/env python 

__author__  = 'Moreno Zolfo (panphlan-users@googlegroups.com)'
__version__ = '0.1'
__date__    = '1 September 2017'


class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	OKGREEN2 = '\033[42m\033[30m'
	RED = '\033[1;91m'
	CYAN = '\033[0;37m'

def fancy_print(mesg,label,type,reline=False,newLine=False):
	opening = "\r" if reline else ''
	ending = "\r\n" if not reline or newLine else ''

	if len(mesg) < 65:
	
		sys.stdout.write(opening+mesg.ljust(66)+(type+'[ - '+label.center(5)+' - ]'+bcolors.ENDC).ljust(14)+ending)
	else: 
		c=0
		wds = []
		lines=[]
		for word in mesg.split(' '):

				if c + len(word)+2 > 65:
					print ' '.join(wds)
					c=0
					wds=[word]
					continue
				c = c+len(word)+2
				wds.append(word)
		sys.stdout.write(opening+(' '.join(wds)).ljust(66)+(type+'[ - '+label.center(5)+' - ]'+bcolors.ENDC).ljust(14)+ending)

	sys.stdout.flush()
	


try:
	import argparse
	import subprocess
	import argparse
	import os
	import sys
	import glob

except ImportError as e:
	print "Error while importing python modules! Remember that this script requires: argparse, subptocess, argparse, os ,sys, glob"
	sys.exit(1)

try:
	from Bio import SeqIO
	from BCBio import GFF
	from Bio.SeqUtils.CheckSum import seguid
except ImportError as e:
	fancy_print("Failed in importing Biopython. Please check Biopython is installed properly on your system!",'FAIL',bcolors.FAIL)
	sys.exit(1)
 

try:
	import pandas as pd 
	import numpy as np
except ImportError as e:
	fancy_print("Failed in importing Pandas and Numpy. Please check they are installed properly on your system!",'FAIL',bcolors.FAIL)
	sys.exit(1)


parser = argparse.ArgumentParser(description='Converts Roary output file into Panphlan Pangenome')
parser.add_argument('--gff_folder', help='The folder containing the reference genomes, in GFF format. Files must end with .gff. Required', required=True)
parser.add_argument('--roary', help='Roary file (gene_presence_absence csv file. Required.', required=True)
parser.add_argument('--index', help='If set, builds the Bowtie2 index instead of a reference fasta file', action='store_true')

parser.add_argument('--out_folder', help='Output folder for pangenome files. By default is the current folder.',default='./')
parser.add_argument('--wname', help='Pangenome name (e.g ecoli)',default='bug')
args = parser.parse_args()

genomeInfo = {} 
trackGenemoes = []

if args.gff_folder: genomesList = glob.glob(args.gff_folder+'/*.gff')
elif args.gff_genome: genomesList = [args.gff_genome]

tpla = []

try:
	if not os.path.isdir(args.out_folder): os.mkdir(args.out_folder)
except OSError as e:
	fancy_print("Failed in writing in "+args.out_folder,'FAIL',bcolors.FAIL)
	sys.exit(1)

fancy_print('Reading '+str(len(genomesList))+' genomes','...',bcolors.OKBLUE)

for fil in genomesList:
	baseName='.'.join(os.path.basename(fil).split('.')[:-1])
	trackGenemoes.append(baseName)

	for contig in GFF.parse(open(fil)):
		contig.description = ''
		contig.id = baseName+'_'+contig.id
		tpla.append(contig)
		for t in contig.features:
			if 'locus_tag' in t.qualifiers:
				genomeInfo[baseName+'__'+t.qualifiers['locus_tag'][0]] = (contig.id,baseName,t.location)
				if len(t.qualifiers['locus_tag']) > 1:
					fancy_print("More than one locus_tag for this annotation! This should not happen.",'???',bcolors.WARNING)

SeqIO.write(tpla,args.out_folder+'/panphlan_'+args.wname+'.fasta','fasta')
if args.index:
	fancy_print('Generating Bowtie2 Index ('+args.out_folder+'/panphlan_'+args.wname+'.fasta)','...',bcolors.OKBLUE)
		
	try:
		subprocess.call(['bowtie2-build',args.out_folder+'/panphlan_'+args.wname+'.fasta',args.out_folder+'/panphlan_'+args.wname,'--quiet'])
		os.remove(args.out_folder+'/panphlan_'+args.wname+'.fasta')
	except OSError as e:
		fancy_print("Failed in running Bowtie2-build!",'FAIL',bcolors.FAIL)
		sys.exit(1)
 
	fancy_print('Generating Bowtie2 Index ('+args.out_folder+'/panphlan_'+args.wname+'.fasta)','DONE',bcolors.OKGREEN)

fancy_print('Generating pangenome file','...',bcolors.OKBLUE)

of = open(args.out_folder+'/panphlan_'+args.wname+'_pangenome.csv','w')

if not os.path.isfile(args.roary):
	fancy_print("Unable to access "+args.roary,'FAIL',bcolors.FAIL)
	sys.exit(1)

ian = pd.read_csv(args.roary,low_memory=False)

for k,ser in ian.iterrows():
	for l in ser.keys():
		if l in trackGenemoes and ser[l] is not np.nan:
			for j in ser[l].split('\t'):
				element = str(l)+'__'+str(j.split('___')[0])
				
				of.write(str(ser['Gene'])+'\t'+str(l)+':'+str(l)+'_'+str(j)+'__'+str(ser['Annotation'])+('##'+str(ser['Non-unique Gene name']) if ser['Non-unique Gene name'] is not np.nan else '')+'\t'+str(l)+'\t'+str(genomeInfo[element][0])+'\t'+str(genomeInfo[element][2].start)+'\t'+str(genomeInfo[element][2].end)+'\n')
of.close()

fancy_print('Task completed. Have a nice day','DONE',bcolors.OKGREEN)
