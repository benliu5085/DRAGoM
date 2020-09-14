==============================================

# DRAGoM: Detection of RNA using Assembly Graph from Metagenomics data

==============================================
# Description:

DRAGoM is a tool designed to predict and assemble ncRNA from next generation sequencing data.
DRAGoM is written in C++ and has been tested on a 64-bit Linux system.
The input for DRAGoM includes:
1, raw reads, pair-end or single-end, trim if necessary;
2, CM files of interested ncRNA from Rfam(https://rfam.xfam.org)
The output of DRAGoM includes:
1, summary.csv: which summarized the read counts for each ncRNA family;
2, predicted reads for each ncRNA family in seperate file;
3, assembled ncRNA homologs for each ncRNA family in seperate file;

==============================================
# Prerequisites:

1. gcc compiler (version > 4.8.5)
2. boost-1.54.0 or newer
3. python2.7

==============================================
# Installation:

To install DRAGoM, please follow the steps below:

1. Untar the downloaded file "DRAGoM.tar.gz". This will generate the directory "DRAGoM".
    $ tar xzvf DRAGoM.tar.gz

2. Install DRAGoM by running the "install.sh" script.
    $ bash install.sh

Note:
To use other version of any third party software listed in DRAGoM/lib:
  1, replace the corresponding path to the executable in "env.config" and,
  2, delete the line starting with "sanity" in "env.config".

==============================================
# Running the program:

1.  The runDragom.py wrapper is used to run the DRAGoM.

USAGE: ./runDragom.py [options]

Input  options:
-1   <filename>    : fastq file with forward paired-end reads
-2   <filename>    : fastq file with reverse paired-end reads
                   only use -1 <filename> means that input reads are interleaved paired-end reads
-s   <filename>    : fastq file with single-end reads
                   if both pair-end or single-end input are provided, the single-end input will be ignored
-p   <filename>    : parameter file
-m   <int>         : [optional] maximum extension length for anchors [default: 100]
-d                 : [optional] anchor masking flag [use this to disable anchor masking]
-k                 : [optional] keep and gzip all intermediate files [use this to keep files]
-h                 : print help message


Note:
1,The parameter file specifies the parameters used for running all programs including in DRAGoM;
2,The parameter file must contain the path to input CM files;
3,You must provide at least one cm file to run the program;
4,Don't change the parameters if you don't understand what it means.

==============================================
# Example:
An example of simulated pair-end reads file is provided in the ~/example/ directory.
You can use the example to test the installation using:
$ cd ~/example/
$ python2 ../runDragom.py -1 example.read1.fq -2 example.read2.fq -p example.config

Note:
1, you might need to unzip the example file using:
$ tar -xvf example.tar.gz

==============================================
# Output:

The final output of the example should contain 5 files.

1. summary.csv, this file summarized the read counts for each ncRNA family:
(Note this is a csv file)
acc,predicted reads,name
RF00002,11,5_8S_rRNA
RF00010,21,RNaseP_bact_a

2, RF00010.assembled_rna.fa, this file listed all homologs of RF00010:
E.g.
>24,39,1,100,7,67,100_0_166_5.8e-28_RF00010
TGGATGATAGATGGAGGAGAGGAAAGTCCGGGCTCCACAGGGCAGGGTGCCAGATAACGTGTGGGGGGTGAAAGCCCACGACCAGTGCAACAGAGAGCAA
ACCGCCGATGGCTTACTTAATGTAGGATCAGGTAAGGGTGAACGGGTGCGGTAAGAGCGCACCGCA
>69,1,1,100_27_404_6.7e-98_RF00010
GGAGTTGACTAGACATTCGCTGCTTTATCATTAATCCTTTGGATGATAGATGGAGGAGAGGAAAGTCCGGGCTCCACAGGGCAGGGTGCCAGATAACGTC
TGGGGGGTGAAAGCCCACGACCAGTGCAACAGAGAGCAAACCGCCGATGGCTTACTTAATGTAGGATCAGGTAAGGGTGAAAGGGTGCGGTAAGAGCGCA
CCGCACGGCTGGTAACAGTTCGTGGCACGGTAAACTCCACTCGGAGCAAGGCCAAATAGGGGTTCATTAGGTACGGCCCGTATCGAACCCGGGTAGGCTG
CTTGAGCCAGTGCGTAAGTGCTGGCCTAGATGAATGATTGTCCACGACAGAACCCGGCTTATCGGTCAACTTCACAA

3. RF00010.prediction.csv, this file listed all predicted reads for RF00010:
(Note this is a csv file)
E.g.
Read_Header,E-value
test_203_341_0:0:0_1:0:0_16/1,6.7e-98
test_114_297_2:0:0_0:0:0_9/1,6.7e-98
test_114_297_2:0:0_0:0:0_9/2,6.7e-98

4, RF00002.assembled_rna.fa, this file listed all homologs of RF00002:
E.g.
>4,33,1,100,11,53,100,43,92,100,49,105,100_123_204_6e-11_RF00002
GACTCTCGGCAACGGATATCTTGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAAGATCC
>12,58,1,100,4,16,100,46,46,100,9,85,100,57,117,100_130_216_1.2e-16_RF00002
AACTTTCAACAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGGGAATTGCAGAATTCAGTG

5. RF00002.prediction.csv, this file listed all predicted reads for RF00002:
(Note this is a csv file)
E.g.
Read_Header,E-value
test_543_727_0:0:0_2:0:0_2/2,1.2e-16
test_512_722_1:0:0_2:0:0_b/1,1.2e-16
test_853_1013_1:0:0_3:0:0_0/2,6e-11

==============================================
