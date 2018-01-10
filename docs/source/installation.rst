Installation
============

Installation of HyperTRIBE's computational pipeline involves installing a set of softwares, downloading a set of annotation files, and updating the shells scripts to provide their locations.


Software Dependencies (tested version)
--------------------------------------
- Trimmomatics (v. 0.36)
- Bowtie2 (v. 2.1.0, 2.2.9)
- STAR (v. 2.5.2b)
- samtools (v. 1.3.1)
- bedtools suite (v. 2.16.2)
- Perl (5.8.8, 5.12.5, 5.22.1) 
- Perl module DBI.pm (1.631, 1.636) and DBD:mysql (4.042)
- MySQL database (MySQL, MariDB)
- Python (2.7.2, other versions should work)
- SRA Toolkit 

Operating systems used: RHEL 5.11 and RHEL 7.2. HyperTRIBE is likely to work seamlessly with the latest version of these software/packages as well.

Source Code
-----------
Download the source code from github.
::

   cd /path_from_root/desired_location
   git clone https://github.com/rosbashlab/HyperTRIBE

Install MySQL (MariaDB) database and Perl modules
-------------------------------------------------
Install MariaDB, here is a helpful `resource <http://tribe-tool.readthedocs.io/en/latest/mariadb.html>`_. Once the root password is set, log in to MySQL using the root password. 
::

    mysql -h localhost -u root -p

If MySQL is installed on a different machine, then replace “localhost” with the IP address of that machine.

Create a database which will be used by HyperTRIBE software. The scripts *load_matrix_data.pl* and *find_rnaeditsites.pl* assumes that a MySQL database called **dmseq** has already been created. If you want to create a different database for this purpose, update the scripts to change the default value of the "$database" variable on line 24 in each perl scripts.
::

    CREATE DATABASE dmseq;   

After exiting MySQL, install MariaDB-devel
::

    yum install mariadb-devel

Install the perl modules using cpan or cpanm if they are not pre-installed in your operating system.
::

    cpan DBI
    cpan DBD::mysql    


Update Perl Scripts
-------------------
Edit the Perl Scripts “load_matrix_data.pl” and “find_rnaeditsites.pl” to update the mysql- related variables that are needed to communicate with the MySQL database tables. Provide MySQL username (in this case “root”), password and name of database, which was created earlier, in the perl scripts load_matrix_data.pl (line 22-25) and find_rnaeditsites.pl (line 23-26). If the MySQL database is hosted on a different machine, change the “$host” variable from “localhost” to the IP address of that machine. 


Download Annotations and Genome Sequences
-----------------------------------------
HyperTRIBE requires transcriptome annotation in two formats and genome sequence in fasta format. These files should be updated by the user based on the organism and genome build of interest. HyperTRIBE need these two annotation files at different step of the pipeline.

**Download refseq annotation for Drosophila** (dm6) from `UCSC Genome Browser <https://genome.ucsc.edu/index.html>`_ Tools=TableBrowser; clade=Insect; genome=D.melanogaster; assembly=dm6; group:Genes and Gene Predictions; track=RefSeq Genes; table=refFlat; output format=all field from selected tables; output file: choose a filename (dm6_refFlat.txt).

**Optional: Recommended procedure for getting the refseq annotation file when the Refseq Genes track is not available as an option in Table Browser.** Here is an example of his alternate approach downloading the refseq annotation for human genome from `UCSC Genome Browser <https://genome.ucsc.edu/index.html>`_ Tools:TableBrowser; clade=Mammal; genome=human; assembly=hg38; group:All Tables; database=hg38; table=refFlat; output format=all field from selected tables; output file: choose a filename (hg38_refFlat.txt). 

**Download refseq annotation in GTF format** from Illumina’s `Igenome page <https://support.illumina.com/sequencing/sequencing_software/igenome.html>`_ . This file can also be downloaded from UCSC browser using the instruction above with one small change, choose output file format as *GTF – gene transfer format*.

These annotation files are provided as part of HyperTRIBE code distribution as an example of required annotation file.


These `annotation files  <https://github.com/rosbashlab/HyperTRIBE/tree/master/annotation>`_ for Drosophila (dm6) are automatically downloaded when  HyperTRIBE source code is cloned. Uncompress the annotation files, which creates a directory with all the annotation files.
::

    cd /path_from_root/HyperTRIBE/annotation
    #uncompress the annotaion files as needed
    gunzip genes.gtf.gz

After downloading the genome sequence from UCSC genome browser or from any other appropriate place, create bowtie2 indices for the genome and STAR indices for the transcriptome. The *genes.gtf* file is used during the creation of STAR indices. 
::

    #Create bowtie2 indices.
    cd /location_of_genome/
    bowtie2-build genome.fa genome
    
Create STAR indices
::

     STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_indices --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf
     #star indices are created at /location_of_genome/star_indices

Update Shell Scripts
--------------------
Update Shell Scripts to reflect the location of software location, annotation and genome indices. Open the shell scripts with a text editor like nano or emacs and update the following lines of code with the location of HyperTRIBE code, annotation files, Bowtie2 and STAR indices.

Edit these variables in **trim_and_align.sh**
::

    star_indices="/path_from_root/star_indices"
    TRIMMOMATIC_JAR="/path_from_root/trimmomatic.jar"
    PICARD_JAR="/path_from_root/picard.jar"

If you want to use a different trimmer or aligner, feel free to change the code

Edit these variables in **trim_and_align_gDNA.sh**
::

    #location of bowtie2 indices
    bowtie_indexes="/path_from_root/genome"
    TRIMMOMATIC_JAR="/path_from_root/trimmomatic.jar"


Edit **load_table.sh**
::

    #location of HyperTRIBE code
    HyperTRIBE_DIR="/path_from_root/HyperTRIBE/CODE"

**Congratulations!!! Now, you are ready to run HyperTRIBE.**



