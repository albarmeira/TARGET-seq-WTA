# TARGET-seq <img align="right" width="250" height="125" src="https://github.com/albarmeira/TARGET-seq/blob/master/target.png">

"Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing"

Rodriguez-Meira, A., Buck, G., Clark, S.-A., Povinelli, B.J., Alcolea, V., Louka, E., McGowan, S., Hamblin, A., Sousos, N., Barkas, N., et al. (2018). Unraveling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing. Molecular Cell.

https://doi.org/10.1016/j.molcel.2019.01.009

# TARGET-seq whole transcriptome analysis

Below are example scripts on how to run preprocessing steps to analyze whole transcriptome TARGET-seq data, from fastq generation, alignment (using STAR) and generation of counts tables using FeatureCounts. 

You can find the pipeline to analyze TARGET-seq single cell genotyping data (SCpipeline) in https://github.com/albarmeira/TARGET-seq/

Author: Alba Rodriguez-Meira.

1.First, demultiplex your files using bcl2fastq (Illumina). Edit RunInfo.xlm file read1, to change read1 to an index read.

```diff
+ Read Number="1" NumCycles="15" IsIndexedRead="Y"
```
2.Run bcl2fastq using a sample sheet containing barcode R1 (cell barcode ) and index read (i7 pool barcode) sequences. For example, considering you are using barcoded oligodT containing a 14 cell barcode sequence, and i7 Illumina indexes (8 bp), you should use the following read configuration R1=15 cycles, I=8 cycles and R2=70 cycles, and run the following command line:

```
module load bcl2fastq/2.20.0.422

bcl2fastq -o output_dir/ --sample-sheet example_sheet.csv --use-bases-mask I14N*,I8,Y70 --no-lane-splitting
```
You can find an example_sample_sheet.csv for bcl2fastq demultiplexing of TARGET-seq fastq files in this repository.

3.Then use STAR to align each * .fastq file:

```
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
```

4.Finally, generate a counts table using FeatureCounts:

```
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
```
