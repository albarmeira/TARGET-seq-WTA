# TARGET-seq <img align="right" width="250" height="125" src="https://github.com/albarmeira/TARGET-seq/blob/master/target.png">

"Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing"

Rodriguez-Meira, A., Buck, G., Clark, S.-A., Povinelli, B.J., Alcolea, V., Louka, E., McGowan, S., Hamblin, A., Sousos, N., Barkas, N., et al. (2018). Unraveling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing. Molecular Cell.

https://doi.org/10.1016/j.molcel.2019.01.009

# TARGET-seq whole transcriptome analysis

Below are example scripts on how to run preprocessing steps to analyze whole transcriptome TARGET-seq data, from fastq generation, alignment (using STAR) and generation of counts tables using FeatureCounts. 

You can find the pipeline to analyze TARGET-seq single cell genotyping data (SCpipeline) in https://github.com/albarmeira/TARGET-seq/

Author: Alba Rodriguez-Meira.

# 1. High throughput 3'-TARGET-seq whole transcriptome analysis

1. 1.First, demultiplex your files using bcl2fastq (Illumina). Edit RunInfo.xlm file read1, to change read1 to an index read.

```diff
+ Read Number="1" NumCycles="15" IsIndexedRead="Y"
```
Run bcl2fastq using a sample sheet containing barcode R1 (cell barcode ) and index read (i7 pool barcode) sequences. For example, considering you are using barcoded oligodT containing a 14 cell barcode sequence, and i7 Illumina indexes (8 bp), you should use the following read configuration R1=15 cycles, I=8 cycles and R2=70 cycles, and run the following command line:

```
module load bcl2fastq/2.20.0.422

bcl2fastq -o output_dir/ --sample-sheet example_sample_sheet_3TARGETseq.csv --use-bases-mask I14N*,I8,Y70 --no-lane-splitting
```
You can find an example_sample_sheet.csv for bcl2fastq demultiplexing of TARGET-seq fastq files in this repository.

1. 2. Trim poly-A tails using TrimGalore:

```
module load trim_galore/0.4.1

mkdir -p galore

for file1 in ../fastq/*R1*.gz
do
        trim_galore --stringency 3 --fastqc --adapter AAAAAAAAAAAAAAAAAAAAA --length 25 --output_dir galore $file1

done

```

1. 3.Then use STAR to align each trimmed * .fastq file:

```
#!/bin/sh
#$ -N starAlign
#$ -cwd

module load rna-star/2.4.2a

genomeDir='/path_to_genome/hg19_STAR' 

#directory where your STAR reference genome is located; use STAR --runMode genomeGenerate to generate such reference genome

mkdir -p output
mkdir -p tmp

for file1 in ../galore/*R1*.gz
do

        outPrefix=`basename $file1 _R1_001_trimmed.fq.gz` #this is the trimmed fastq file prefix

        STAR --runThreadN 4 \
                --genomeLoad LoadAndKeep \
                --genomeDir  $genomeDir \
                --readFilesIn $file1 \
                --readFilesCommand zcat \
                --outFileNamePrefix output/$outPrefix \
                --outTmpDir tmp/$outPrefix \
                --outReadsUnmapped Fastx \
                --outSAMtype BAM Unsorted
done

STAR --genomeDir $genomeDir --genomeLoad Remove
```

1. 4.Finally, generate a counts table using FeatureCounts:

```
#!/bin/sh
#$ -cwd
#$ -N featureCount

module load subread/1.4.5-p1

featureCounts --primary -T 4 -a /path_to_annotation/annotation.gtf  -o counts.txt STAR/output/*.bam
```

# 2. Full length TARGET-seq whole transcriptome analysis

2. 1. Full-length TARGET-seq dataset contain two index reads (i7/i5). First, demultiplex your files using bcl2fastq (Illumina). Run bcl2fastq using a sample sheet containing index read 1 (i7) and index read 2 (i5) sequences.

If using single-end 75 cycle reads:

```
module load bcl2fastq/2.20.0.422

bcl2fastq -o output_dir/ --sample-sheet example_sample_sheet_FLTARGETseq.csv --no-lane-splitting
```
If using paired-end 75 cycle reads:

```
module load bcl2fastq/2.20.0.422

bcl2fastq -o output_dir/ --sample-sheet example_sample_sheet_FLTARGETseq.csv --no-lane-splitting
```

You can find an example_sample_sheet_FLTARGETseq.csv for bcl2fastq demultiplexing of full-length TARGET-seq fastq files in this repository.

2. 2. Trim Nextera adaptors using TrimGalore. If using single-end reads, run:

```
module load trim_galore/0.4.1

mkdir -p galore

for file1 in ../fastq/*R1*.gz
do
        trim_galore --stringency 3 --fastqc --nextera --output_dir galore $file1

done

```

If using paired-end reads, run:

```
module load trim_galore/0.4.1

mkdir -p galore

for file1 in ../fastq/*R1*.gz
do
        file2name=${file1%_R1*}
        file2="${file2name}_R2_001.fastq.gz" #this is the prefix of the R2 fastq file
        trim_galore --stringency 3 --fastqc --nextera --output_dir galore --paired $file1 $file2

done

```

2. 3. Then use STAR to align each * .fastq file. If using single-end reads, run:

```
#!/bin/sh
#$ -N starAlign
#$ -cwd

module load rna-star/2.4.2a

genomeDir='/path_to_genome/hg19_STAR' 

#directory where your STAR reference genome is located; use STAR --runMode genomeGenerate to generate such reference genome

mkdir -p output
mkdir -p tmp

for file1 in ../galore/*R1*.gz
do

        outPrefix=`basename $file1 _R1_001_trimmed.fq.gz` #this is the trimmed fastq file prefix

        STAR --runThreadN 4 \
                --genomeLoad LoadAndKeep \
                --genomeDir  $genomeDir \
                --readFilesIn $file1 \
                --readFilesCommand zcat \
                --outFileNamePrefix output/$outPrefix \
                --outTmpDir tmp/$outPrefix \
                --outReadsUnmapped Fastx \
                --outSAMtype BAM Unsorted
done

STAR --genomeDir $genomeDir --genomeLoad Remove
```

If using paired-end reads, run:

```
#!/bin/sh
#$ -N starAlign
#$ -cwd

module load rna-star/2.4.2a

genomeDir='/path_to_genome/hg19_STAR' 

#directory where your STAR reference genome is located; use STAR --runMode genomeGenerate to generate such reference genome

mkdir -p output
mkdir -p tmp

for file1 in ../galore/*R1*.gz
do

        outPrefix=`basename $file1 _R1_001_trimmed.fq.gz`
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

2. 4. Finally, generate a counts table using FeatureCounts. If using single-end reads, run:

```
#!/bin/sh
#$ -cwd
#$ -N featureCount

module load subread/1.4.5-p1

featureCounts --primary -T 4 -a /path_to_annotation/annotation.gtf  -o counts.txt STAR/output/*.bam
```

If using paired-end reads, run:

```
#!/bin/sh
#$ -cwd
#$ -N featureCount

module load subread/1.4.5-p1

featureCounts -p --primary -T 4 -a /path_to_annotation/annotation.gtf  -o counts.txt STAR/output/*.bam
```

