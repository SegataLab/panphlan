# RNA-seq: get transcripts that align with putatively non-­coding DNA sequences of unknown gene families
# get reads that maps against the reference genomes but are not assigned to any of the clustered  gene  families

# please adapt to following variables
SAMPLE_PATH=Samples/G38878.tar.bz2
SPECIES=ecoli14    # panphlan species database
GENE_SEQ_PATH=ecoli14/ffn  # folder of gene sequences (.ffn files of pangenome database)
NUMBER_OF_PROCESSORS=12



###############################################################################
# Step A) generate bowtie2 index of all gene sequences
# all genes from all reference genomes of the species database
cat ${GENE_SEQ_PATH}/*.ffn > genes.fna
bowtie2-build genes.fna genes_${SPECIES}



###############################################################################
# Step B) get non-coding reads of a sample

SAMPLE_ID=`basename ${SAMPLE_PATH%%.*}`;

# panphlan mapping against reference genomes, get bam (option --out_bam)
./panphlan/panphlan_map.py -c ${SPECIES} -i ${SAMPLE_PATH} -o Cov/${SAMPLE_ID} --out_bam Bam/${SAMPLE_ID} -p ${NUMBER_OF_PROCESSORS}

# convert bam to fastq
samtools bam2fq Bam/${SAMPLE_ID}.bam > ${SAMPLE_ID}.fastq

# map all genome related reads against gene-sequences
bowtie2 -x genes_${SPECIES} -U ${SAMPLE_ID}.fastq -S ${SAMPLE_ID}_mapped_and_unmapped.sam -p ${NUMBER_OF_PROCESSORS}
# convert sam to bam
samtools view -bS ${SAMPLE_ID}_mapped_and_unmapped.sam > ${SAMPLE_ID}_mapped_and_unmapped.bam

# get reads not mapped to gene-sequences (but mapped to reference genomes)
samtools view -b -f 4 -F 256 ${SAMPLE_ID}_mapped_and_unmapped.bam > ${SAMPLE_ID}_non_coding_reads.bam

echo -e '\n Result: '${SAMPLE_ID}'_non_coding_reads.bam \n'

echo 'Number of reads:'
samtools view -c ${SAMPLE_ID}_non_coding_reads.bam


