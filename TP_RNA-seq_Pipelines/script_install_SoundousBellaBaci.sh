# -- Script d'installation des outils et téléchargement des fichiers
# TP1 - RNAseq
# Soundous Bella Baci

#!/bin/bash

conda init bash
conda activate

# Creating a directory that will contain all the necessary files for this pipeline
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|   Step 1/3 : Creating a working directory to host all files for this session  |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Check if file exist, if not create them
if [ ! -d TP_rnaseq_SBB ];then
  mkdir TP_rnaseq_SBB
fi

echo -e "\e[1;96mThe directory for this session is called : ./TP_rnaseq_SBB  \e[0m\n"

cd TP_rnaseq_SBB

# Installation of the programs needed for the pipeline

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                         Step 2/3 : Installing programs                        |\e[0m"
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

# Installing STAR
echo -e "\e[1;31mInstalling Star : \e[0m\n"
conda install -c bioconda star
echo -e "\e[1;32mStar is installed on your machine!\e[0m\n"

# Installing samtools
echo -e "\e[1;31mInstalling samtools : \e[0m\n"
sudo apt install samtools
echo -e "\e[1;32mSamtools is installed on your machine!\e[0m\n"

# Installing featureCounts
echo -e "\e[1;31mInstalling featureCounts :\e[0m\n "
sudo apt-get update -y
sudo apt-get install -y subread
echo -e "\e[1;32mFeatureCounts is installed on your machine!\e[0m\n"


# Data retrieval
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                             Step 3/3 : Extracting data                        |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

## -- Récupération des fichiers fastq

# Check if file exist, if not create them
if [ ! -d fastq_files ];then
  mkdir fastq_files
  echo -e "\e[1;96mCreating directory for the fastq files: ./fastq_files \e[0m\n"
fi

cd fastq_files

echo -e "\e[1;96mDownloading the fastq files :\e[0m\n"
wget  http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz
tar -zxvf TPrnaseq.tar.gz

# Check if the fastq files have been downloaded
sample_list=(Day_0_1_chr18.sampled Day_0_2_chr18.sampled Day_0_3_chr18.sampled Day_7_1_chr18.sampled Day_7_2_chr18.sampled Day_7_3_chr18.sampled)

for i in ${sample_list[*]}
do
  if ! test -f ./$i'.R1.fastq' && ! test -f ./$i'.R2.fastq'; then
    echo -e "\n\e[1;31mSomething went wrong with the download or the unzipping of the fastq files. Please check what went wrong.\e[0m\n"
    exit
  fi
done

echo -e "\n\e[1;32mThe fastq files has been downloaded and unzipped.\e[0m\n"
cd ..


# -- Recupération du génome
if [ ! -d genome ];then
  mkdir genome
  echo -e "\e[1;96mCreating directory for the genome (chromosome 18) : ./genome \e[0m\n"
fi

cd genome

echo -e "\e[1;96mDownloading the chromosome 18 file :\e[0m\n"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz
gunzip chr18.fa.gz

if test -f ./chr18.fa ; then
  echo -e "\e[1;32mThe chromosome 18 has been downloaded and unzipped.\e[0m\n"
else
  echo -e "\n\e[1;31mSomething went wrong with the download or the unzipping of the chromosome file. Please check what went wrong.\e[0m\n"
  exit
fi

cd ..

# -- Recupération de l'annotation

if [ ! -d annotation ];then
  mkdir annotation
  echo -e "\e[1;96mCreating directory for the annotation file : ./annotation \e[0m\n"
fi

cd annotation

echo -e "\e[1;96mDownloading the annotation file :\e[0m\n"
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz
gunzip gencode.v24lift37.basic.annotation.gtf.gz

if test -f ./gencode.v24lift37.basic.annotation.gtf ; then
  echo -e "\e[1;32mThe annotation file has been downloaded and unzipped.\e[0m\n"
else
  echo -e "\n\e[1;31mSomething went wrong with the download or the unzipping of the annotation file. Please check what went wrong.\e[0m\n"
  exit
fi

cd ..
