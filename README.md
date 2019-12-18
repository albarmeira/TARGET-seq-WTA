# TARGET-seq-WTA

Example scripts to perform whole transcriptome analysis of TARGET-seq datasets


1.	First, demultiplex your files using bcl2fastq (Illumina). Edit RunInfo.xlm file read1, to change read1 to an index read.

<Read Number="1" NumCycles="15" IsIndexedRead="Y" />

2.	Run bcl2fastq using a sample sheet containing barcode R1 (cell barcode ) and index read (i7 pool barcode) sequences. For example, considering you are using barcoded oligodT containing a 14 cell barcode sequence, and i7 Illumina indexes (8 bp), you should use the following read configuration R1=15 cycles, I=8 cycles and R2=70 cycles, and run the following command line:

module load bcl2fastq/2.20.0.422

bcl2fastq -o output_dir/ --sample-sheet example_sheet.csv --use-bases-mask I14N*,I8,Y70 --no-lane-splitting

3.	Then use STAR to align each *.fastq file:

#!/bin/sh
#$ -N starAlign
#$ -cwd

module load rna-star/2.4.2a

genomeDir='/path_to_genome/hg19_STAR' 

#directory where your STAR reference genome is located; use STAR --runMode genomeGenerate to generate such reference genome

mkdir -p output
mkdir -p tmp

for file1 in ../fastq/*R1*.gz
do

        outPrefix=`basename $file1 _R1_001.fastq.gz`
        mkdir output/$outPrefix
        file2name=${file1%_R1*}
        file2="${file2name}_R2_001.fastq.gz"

        STAR --runThreadN 4 \
                --genomeLoad LoadAndKeep \
                --genomeDir  $genomeDir \
                --readFilesIn $file1 $file2 \
                --readFilesCommand zcat \
                --outFileNamePrefix output/$outPrefix \
                --outTmpDir tmp/$outPrefix \
                --outReadsUnmapped Fastx \
                --outSAMtype BAM Unsorted
done

STAR --genomeDir $genomeDir --genomeLoad Remove

4.	Finally, generate a counts table using FeatureCounts:

#!/bin/sh
#$ -N SortBam
#$ -cwd

module load samtools/0.1.19

mkdir -p output_counts

for file in ../STAR/output/*.bam
do
        prefix=`basename ${file}`
        samtools sort -@ 4 -m 4G $file output_counts/$prefix.sorted
        samtools index output_counts/$prefix.sorted.bam

done

