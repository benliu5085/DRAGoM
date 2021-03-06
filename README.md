==============================================

# DRAGoM: Detection of RNA using Assembly Graph from Metagenomics data

==============================================
# Description:

DRAGoM is a tool designed to predict and assemble ncRNA from next generation sequencing data.
DRAGoM is written in C++ and has been tested on a 64-bit Linux system, Red hat (7) and Ubuntu (20.04 LTS).
DRAGoM is freely available under the Creative Commons BY-NC-ND 4.0 License Agreement (https://creativecommons.org/licenses/by-nc-nd/4.0/).
If you have any troubles when using DRAGoM, feel free to contact Ben Liu (ben_0522@ku.edu) or Cuncong Zhong (cczhong@ku.edu).

The input for DRAGoM includes:
1. raw reads, pair-end or single-end, trim if necessary;
2. CM files of interested ncRNA from Rfam(https://rfam.xfam.org)

The output of DRAGoM includes:
1. summary.csv: which summarized the read counts for each ncRNA family;
2. predicted reads for each ncRNA family in seperate file;
3. assembled ncRNA homologs for each ncRNA family in seperate file;

==============================================
# Prerequisites:

1. g++ compiler ( >= 4.8.5)
2. boost ( >= 1.53)
3. python 2.7
4. make ( >= 3.82)
5. cmake ( >= 3.12)

Note:
We assume the boost library is installed using package manager to default directory (/usr/include/boost). If you installed boost library at other directory, please update the ~/DRAGoM/src/CMakeLists.txt "set(Boost_INCLUDE_DIR /usr/include)" into "set(Boost_INCLUDE_DIR WHERE/BOOST/IS)", WHERE/BOOST/IS is the directory where you installed boost library.

# Third-party software:
1. sga (https://github.com/jts/sga)
2. SPAdes (https://cab.spbu.ru/software/spades/)
3. bwa (http://bio-bwa.sourceforge.net)
4. infernal (http://eddylab.org/infernal/)
5. cd-hit (http://weizhongli-lab.org/cd-hit/)
6. samtools (http://www.htslib.org)

Note:
1. You must install the prerequisites to successfully run DRAGoM.
2. We provided one set of pre-built binary of all third-party softwares with DRAGoM.
   If you prefer to use another version of any third-party software, remember to update the env.config file.
3. We have tested the softwares on Red hat ( > 7 ) and Ubuntu ( > 18.04 ), if you are going to use DRAGoM on
   other OS, you might need to install the third-party softwares on your own.

==============================================
# Installation:

To install DRAGoM, please follow the steps below:

1. Untar the downloaded file "DRAGoM.tar.gz". This will generate the directory "DRAGoM".
    $ tar xzvf DRAGoM.tar.gz

2. Install DRAGoM by running the "install.sh" script.
    $ bash install.sh

Note:
1. An 'env.config' file will be created automatically after runninng 'install.sh', which contains
   the absolute directory to all executable of all third-party softwares. If you want to use another
   version, update the corresponding directory to the executable you want DRAGoM to use.
2. The 'cmsearch' and 'cmpress' are included in the software infernal.

==============================================
# Running the program:

1.  The runDragom.py wrapper is used to run the DRAGoM.

USAGE: python2 runDragom.py [options]

Input  options:
-1   <filename>    : fastq file with forward paired-end reads
-2   <filename>    : fastq file with reverse paired-end reads
                   only use -1 <filename> means that input reads are interleaved paired-end reads
-s   <filename>    : fastq file with single-end reads
                   if both pair-end or single-end input are provided, the single-end input will be ignored
-p   <filename>    : parameter file
-m   <int>         : [optional] maximum extension length for anchors [default: 100]
-d                 : [optional] anchor masking flag [use this to disable anchor masking]
-h                 : print help message


Note:
1. The parameter file specifies the parameters used for running all programs including in DRAGoM;
2. The parameter file must contain the path to input CM files;
3. You must provide at least one cm file to run the program;
4. Don't change the parameters if you don't understand what it means.

==============================================
# Example:
An example of simulated pair-end reads file is provided in the ~/example/ directory.
You can use the example to test the installation using:
$ cd ~/example/
$ python2 ../runDragom.py -1 example.read1.fq -2 example.read2.fq -p example.config

Note:
you might need to unzip the example file using:
$ tar -xvf example.tar.gz

==============================================
# Output:

The final output of the example should contain 5 files.

1. summary.csv, this file summarized the read counts for each ncRNA family:
(Note this is a csv file)
acc,predicted reads,name
RF00002,11,5_8S_rRNA
RF00010,21,RNaseP_bact_a

2. RF00010.assembled_rna.fa, this file listed all homologs of RF00010:
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

4. RF00002.assembled_rna.fa, this file listed all homologs of RF00002:
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
