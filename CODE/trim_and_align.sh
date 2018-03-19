#!/bin/sh

#----- Update the following varibles
star_indices="/home/analysis/genome/dm6/star_indices"
TRIMMOMATIC_JAR="/opt/TRIMMOMATIC/0.36/trimmomatic.jar"
PICARD_JAR="/opt/PICARD/2.8.2/picard.jar"
#GENOME_FAI_FILE="/home/analysis/genome/dm6/genome.fa.fai"
#---------End Update variable------------

# Use the loop below to run the code for multiple fastq files
#filelist=`ls *.fastq`#
#for file in ${filelist[@]}
#do

# the input fastq file should have .fastq prefix, for example s2_mRNA.fastq or HyperTRIBE_rep1.fastq
file=$1
prefix=${file%.fastq*}
trim_input=$file
trim_outfile=$prefix.trim.fastq 
avgquality="25"

# The first six nucleotide of the read is removed due to pontential error from random hexamer mispriming. Trimmomatics is used to remove low quality bases from either end of the reads. Reads with avg quality score of less than 25 are also removed. Please free to use similar software instead of Trimmomatics

java -jar $TRIMMOMATIC_JAR SE -phred33 $trim_input $trim_outfile HEADCROP:6 LEADING:25 TRAILING:25 AVGQUAL:$avgquality MINLEN:19

# Align library with STAR
input=$trim_outfile
STAR  --runThreadN 8 --outFilterMismatchNoverLmax 0.07 --outFileNamePrefix $prefix"_" --outFilterMatchNmin 16 --outFilterMultimapNmax 1  --genomeDir $star_indices --readFilesIn $input
# --outFilterMismatchNoverLmax 0.07: number of mismatches is <= 7% of mapped read length
# --outFilterMatchNmin 16: min numberof bases mapped genome per read
# --outFilterMultimapNmax 1: output reads that only map to one loci

output=$prefix".sam"
mv $prefix"_"Aligned.out.sam $output

#following code is tested with samtools 1.3.1, you might have to tweak it a bit bases your installed verison of samtools (these flags can be problematic for older version of samtools: -@, -o)
#remove low quality alignment
samtools view -@ 4 -Sh -q 10 $output > $prefix"_highquality.sam"
mv $prefix"_highquality.sam" $output
bam_out=$prefix".bam"
#convert sam to bam 
samtools view -@ 4 -bhS $output > $bam_out
rm $output
#sort the bam file before using picard
sort_out=$prefix".sort.bam"
samtools sort -@ 6 $bam_out -o $sort_out
rm $bam_out

#run Picard to remove duplicates
input_for_picard=$sort_out
dupremove_bam=$prefix"_nodup.bam"
java -Xmx4g -jar $PICARD_JAR MarkDuplicates INPUT=$input_for_picard OUTPUT=$dupremove_bam METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=tmp ASSUME_SORTED=true
rm $input_for_picard

# sort the output bam file from picard
sort_out=$prefix".sort.bam"
samtools sort -@ 6 $dupremove_bam -o $sort_out
rm $dupremove_bam

# The next step of HyperTRIBE requires the sam file to be sorted 
#Create a SAM file from this sorted bam file
samtools view -@ 4 -h $sort_out > $prefix".sort.sam"
samtools index $sort_out
#------------
echo "Done with STAR mapping and PCR duplicate removal with PICARD"
echo "created sam file: $prefix.sam"
#-------------
## end of do for fastq file
#done

