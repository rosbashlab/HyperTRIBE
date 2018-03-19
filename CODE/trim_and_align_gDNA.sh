#!/bin/sh

#----- Update the following varibles ----
bowtie_indexes="/home/analysis/genome/dm6/genome"
TRIMMOMATIC_JAR="/opt/TRIMMOMATIC/0.36/trimmomatic.jar"
#GENOME_FAI_FILE="/home/analysis/genome/dm3/Sequence/bowtie_indexes/genome.fa.fai"
#------End update variable -----

file=$1
prefix=${file%.fastq*}
trim_input=$file
trim_outfile=$prefix.trim.fastq 
avgquality="25"

# ----------- trim low quality bases and remove low quality reads ---------
# The first six nucleotide of the read is removed due to pontential error from random hexamer mispriming. Trimmomatics is used to remove low quality bases from either end of the reads. Reads with avg quality score of less than 25 are also removed. Please free to use similar software instead of Trimmomatics

java -jar $TRIMMOMATIC_JAR SE -phred33 $trim_input $trim_outfile HEADCROP:6 LEADING:25 TRAILING:25 AVGQUAL:$avgquality MINLEN:19

#--- align gDNA to genome with bowtie2 -----
input=$trim_outfile
bowtie2_out=$prefix"_unfilterrm.sam"
bowtie2 --sensitive -p 9 -x $bowtie_indexes -U $input -S $bowtie2_out

#following code is tested with samtools 1.3.1, you might have to tweak it a bit bases your installed verison of samtools (these flags can be problematic for older version of samtools: -@, -o)
samtools view -@ 4 -Sh -q 10 $bowtie2_out >$prefix".sam"
rm $bowtie2_out

samtools view -@ 4 -bSh  $prefix".sam" >  $prefix".bam"

# sort the bam file
sort_out=$prefix".sort.bam"
samtools sort -@ 6 $prefix".bam" -o $sort_out

# The next step of HyperTRIBE requires the sam file to be sorted 
#Create a SAM file from this sorted bam file
samtools view -@ 4 -h $sort_out > $prefix".sort.sam"



rm $prefix".sam"
rm $prefix".bam"
