# -- Script d'installation des outils et téléchargement des fichiers
# TP 2 - Exome Sequencing Pipeline
# Soundous Bella Baci

#!/bin/bash

conda init bash
conda activate

# Creating a directory that will contain all the necessary files for this pipeline
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|   Step 1/4 : Creating a working directory to host all files for this session  |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Check if file exist, if not create them
if [ ! -d TP_variant_SBB ];then
  mkdir TP_variant_SBB
fi

echo -e "\e[1;96mThe directory for this session is called : ./TP_variant_SBB  \e[0m\n"

cd TP_variant_SBB

# Data retrieval
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                             Step 2/4 : Extracting data                        |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Downloading the fastq files from an online repository
echo -e "\e[1;96mDownloading the fastq files :\e[0m\n"
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt
echo -e "\e[1;33mDon't mind the 'tar: Ignoring unknown extended header keyword' messages, it's because the files were created on macOS and are now being extracted on Linux.\e[0m \n"
# uncompressing the file
tar -zxvf patient7.tar.gz
gunzip ./patient7.exome/TCRBOA7*

# Check if the fastq files have been downloaded
if test -f ./patient7.exome/TCRBOA7-N-WEX-chr16_r1F.fastq && test -f ./patient7.exome/TCRBOA7-N-WEX-chr16_r2F.fastq && test -f ./patient7.exome/TCRBOA7-T-WEX-chr16_r1F.fastq && test -f ./patient7.exome/TCRBOA7-T-WEX-chr16_r2F.fastq; then
  echo -e "\n\e[1;32mThe fastq files has been downloaded and unzipped.\e[0m\n"
else
  echo -e "\n\e[1;31mSomething went wrong with the download or the unzipping of the fastq files. Please check what went wrong.\e[0m\n"
  exit
fi

# Check if file exist, if not create them
if [ ! -d index ];then
  mkdir index
  echo -e "\e[1;96mCreating directory for the index : ./index \e[0m\n"
fi

cd index

# Chromosome 16 retrieval from an online source
echo -e "\e[1;96mDownloading the human chromosome 16 in fasta format :\e[0m\n"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz
# uncompressing the file
gunzip chr16.fa.gz

# Check if the chromosome file have been downloaded
if test -f ./chr16.fa ; then
  echo -e "\e[1;32mThe chromosome 16 has been downloaded and unzipped.\e[0m\n"
else
  echo -e "\n\e[1;31mSomething went wrong with the download or the unzipping of the chromosome file. Please check what went wrong.\e[0m\n"
  exit
fi

cd ..

# Check if file exist, if not create them
if [ ! -d annotation ];then
  mkdir annotation
  echo -e "\e[1;96mCreating directory for the annotation : ./annotation \e[0m\n"
fi

cd annotation

# Download the annotation file from the ebi website
echo -e "\e[1;96mDownloading the annotattion file :\e[0m\n"
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz
# uncompressing the file
gunzip gencode.v24lift37.basic.annotation.gtf.gz

# Check if the annotation was done
if test -f ./gencode.v24lift37.basic.annotation.gtf ; then
  echo -e "\e[1;32mThe annotation file has been downloaded and unzipped.\e[0m\n"
else
  echo -e "\n\e[1;31mSomething went wrong with the download or the unzipping of the annotation file. Please check what went wrong.\e[0m\n"
  exit
fi

cd ..

# Installation of the programs needed for the pipeline

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                         Step 3/4 : Installing programs                        |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Installing fastQC
echo -e "\e[1;31mInstalling fastQC : \e[0m\n "
conda install -c bioconda fastqc
echo -e "\e[1;32mFastQC is installed on your machine!\e[0m\n"

# Installing trimmomatic
echo -e "\e[1;31mInstalling trimmomatic : \e[0m\n "
conda install -c bioconda trimmomatic
echo -e "\e[1;32mTrimmomatic is installed on your machine!\e[0m\n"

# Installing bwa
echo -e "\e[1;31mInstalling bwa : \e[0m\n"
conda install -c bioconda bwa
echo -e "\e[1;32mBwa is installed on your machine!\e[0m\n"

# Installing samtools
echo -e "\e[1;31mInstalling samtools : \e[0m\n"
sudo apt install samtools
echo -e "\e[1;32mSamtools is installed on your machine!\e[0m\n"

# Installing varscan
echo -e "\e[1;31mInstalling varscan :\e[0m\n "
conda install -c bioconda varscan
echo -e "\e[1;32mVarscan is installed on your machine!\e[0m\n"

# Installing bedtools
echo -e "\e[1;31mInstalling bedtools :\e[0m\n "
conda install -c bioconda bedtools
echo -e "\e[1;32mBedtools is installed on your machine!\e[0m\n"

# Creating an index of the partial genome (chromosome 16)

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                         Step 4/4 : Creating a BWA index                       |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Indexing
bwa index -a bwtsw ./index/chr16.fa

# Check if the index was done
if test -f ./index/chr16.fa.amb && test -f ./index/chr16.fa.ann && test -f ./index/chr16.fa.bwt && test -f ./index/chr16.fa.pac && test -f ./index/chr16.fa.sa; then
  echo -e "\n\e[1;32mThe index has been created in the directory : /TP_variant_SBB/index \e[0m\n"
else
  echo -e "\n\e[1;31mSomething went wrong with the indexing. Please check what went wrong.\e[0m\n"
  exit
fi

echo -e "\e[1;46;30m-----------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;46;30m|                                                                                 |\e[0m"
echo -e "\e[1;46;30m|You have set up everything you need, you can now run the exome sequencing script.|\e[0m"
echo -e "\e[1;46;30m|                                                                                 |\e[0m"
echo -e "\e[1;46;30m-----------------------------------------------------------------------------------\e[0m"
