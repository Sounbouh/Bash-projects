# -- Script du pipeline du séquençage d'exome
# TP 2 - Exome Sequencing Pipeline
# Soundous Bella Baci

#!/bin/bash

cd TP_variant_SBB
echo -e "\e[1;96mThe directory for this session is called : ./TP_variant_SBB  \e[0m\n"


# -- 1/ Data source

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                           Step 1/7 : Data source                              |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Information about the data used in the pipeline
echo -e "\e[1;33mIn this study, we are working on data from patient ID TCRBOA7 from the Open Access pancreatic cancer dataset from Baylor College of Medicine, Houston.\e[0m"
echo -e "\e[1;33mSequences mapping to Chromosome 16 have been extracted from the original file to simulate the actual output from a sequencing platform.\e[0m"
echo -e "\e[1;33mN stands for normal sample and T stands for tumoral sample.\e[0m\n"

# -- 2/ Downloading fastq files

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                     Step 2/7 : Inspecting the fastq files                     |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Samples lists of the fastq files names
sample_list=(TCRBOA7-N-WEX-chr16_r1F TCRBOA7-N-WEX-chr16_r2F TCRBOA7-T-WEX-chr16_r1F TCRBOA7-T-WEX-chr16_r2F)
sample_list_2=(TCRBOA7-N-WEX-chr16 TCRBOA7-T-WEX-chr16)

# Checking file sizes and estimating read numbers for each sample :
for i in ${sample_list[*]}
do
  echo -e "\e[1;32mThe file $i contains $(wc -l < ./patient7.exome/$i'.fastq') reads.\e[0m\n"
done

# -- 3/ Quality control + trimming
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                   Step 3/7 : Quality control and trimming                     |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

echo -e "\e[1;96mCreating fastQC files for all the samples : \e[0m\n"

for i in ${sample_list[*]}
do
  # Check if file exist, if not create them
  if [ ! -d fastQC_files ];then
    mkdir fastQC_files
    echo -e "\e[1;96mCreating a directory : ./TP_variant_SBB/fastQC_files \e[0m\n"
  fi

  # Quality control : generating fastQC files for each fastq files
  fastqc -o fastQC_files ./patient7.exome/$i'.fastq'

  # Check if the files were generated, if not stop the script
  if test -f ./fastQC_files/$i'_fastqc.html'; then
    echo -e "\n\e[1;32mFastQC files were created for the file $i. You can find them in the directory : fastQC_files\e[0m\n"
  else
    echo -e "\n\e[1;31mSomething went wrong with the fastQC command. Please check what went wrong.\e[0m\n"
    exit
  fi
done

echo -e "\e[1;96mEliminating low quality sequences at read extermities : \e[0m\n"

for i in ${sample_list_2[*]}
do
  # Check if file exist, if not create them
  if [ ! -d trimmed_fastq ];then
    mkdir trimmed_fastq
    echo -e "\e[1;96mCreating a directory : ./TP_variant_SBB/trimmed_fastq \e[0m\n"
  fi

  # Trimming the reads extremities if they have a low quality
  trimmomatic PE -threads 8 ./patient7.exome/$i'_r1F.fastq' ./patient7.exome/$i'_r2F.fastq' -baseout ./trimmed_fastq/$i'.fastq' LEADING:20 TRAILING:20 MINLEN:50

  # Check if the files were generated, if not stop the script
  if test -f ./trimmed_fastq/$i'_1P.fastq' && test -f ./trimmed_fastq/$i'_2P.fastq' && test -f ./trimmed_fastq/$i'_1U.fastq' && test -f ./trimmed_fastq/$i'_2U.fastq'; then
    echo -e "\n\e[1;32mTrimmed fastq created for the file $i. You can find them in the directory : trimmed_fastq\e[0m\n"
  else
    echo -e "\n\e[1;31mSomething went wrong with the trimmomatic command. Please check what went wrong.\e[0m\n"
    exit
  fi
done

# -- 5/ Mapping with BWA

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                            Step 4/7 : Mapping with BWA                        |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

for i in ${sample_list_2[*]}
do
  # Check if file exist, if not create them
  if [ ! -d mapping ];then
    mkdir mapping
    echo -e "\e[1;96mCreating a directory for the mapping : ./TP_variant_SBB/mapping \e[0m\n"
  fi

  # Mapping the fastq files to the reference chromosome 16
  echo -e "\e[1;96mMapping for $i :\e[0m\n"
  bwa mem -M -t 8 -A 2 -E 1 ./index/chr16.fa ./trimmed_fastq/$i'_1P.fastq' ./trimmed_fastq/$i'_2P.fastq' > ./mapping/$i'.sam'

  # Check if the files were generated, if not stop the script
  if test -f ./mapping/$i'.sam'; then
    echo -e "\n\e[1;32mMapping done for $i. You can find it in the directory : mapping\e[0m\n"
  else
    echo -e "\n\e[1;31mSomething went wrong with the mapping. Please check what went wrong.\e[0m\n"
    exit
  fi
done


# -- 6/ Processing SAM files

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                          Step 5/7 : Processing SAM files                      |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Processing the mapped sam file with samtools
for i in ${sample_list_2[*]}
do
  echo -e "\e[1;96mProcessing for $i :\e[0m\n"

  # Converting sam file to bam file
  echo -e "\e[1;96mConverting Sam file to Bam file :\e[0m"
  samtools view -@ 8 -S -b ./mapping/$i'.sam' > ./mapping/$i'.bam'

  # Check if the files were generated, if not stop the script
  if test -f ./mapping/$i'.bam'; then
    echo -e "\e[1;32mDone \e[0m\n"
  else
    echo -e "\n\e[1;31mThe file has not been created. Please check what went wrong.\e[0m\n"
    exit
  fi

  # Display the statistics of the mapped file
  echo -e "\e[1;96mNumber of alignments for each FLAG type :\e[0m\n"
  samtools flagstat ./mapping/$i'.bam'
  echo -e "\n"

  # Sorting the bam file
  echo -e "\e[1;96mSorting Bam file :\e[0m"
  samtools sort -@ 8 ./mapping/$i'.bam' > ./mapping/$i'_sorted.bam'

  # Check if the files were generated, if not stop the script
  if test -f ./mapping/$i'_sorted.bam'; then
    echo -e "\e[1;32mDone \e[0m\n"
  else
    echo -e "\n\e[1;31mThe file has not been created. Please check what went wrong.\e[0m\n"
    exit
  fi

  # Indexing the bam file
  echo -e "\e[1;96mIndexing bam file :\e[0m"
  samtools index -@ 8 ./mapping/$i'_sorted.bam'

  # Check if the files were generated, if not stop the script
  if test -f ./mapping/$i'_sorted.bam.bai'; then
    echo -e "\e[1;32mDone \e[0m\n"
  else
    echo -e "\n\e[1;31mThe file has not been created. Please check what went wrong.\e[0m\n"
    exit
  fi

  # Creating a pileup file
  echo -e "\e[1;96mConverting bam file to pileup format :\e[0m"
  samtools mpileup -B -A -f ./index/chr16.fa  ./mapping/$i'_sorted.bam' > ./mapping/$i'_pileupfile'

  # Check if the files were generated, if not stop the script
  if test -f ./mapping/$i'_pileupfile'; then
    echo -e "\e[1;32mDone \e[0m\n"
  else
    echo -e "\n\e[1;31mThe file has not been created. Please check what went wrong.\e[0m\n"
    exit
  fi
done

# -- 9/ Calling somatic variants with Varscan

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                 Step 6/7 : Calling somatic variants with Varscan              |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# Check if file exist, if not create them
if [ ! -d varscan_output ];then
  mkdir varscan_output
  echo -e "\e[1;96mCreating a directory : ./TP_variant_SBB/varscan_output \e[0m\n"
fi

echo -e "\e[1;96mCreating variant lists from Normal and Tumor Pileup files :\e[0m\n"

# Create variant list
varscan somatic ./mapping/TCRBOA7-N-WEX-chr16_pileupfile ./mapping/TCRBOA7-T-WEX-chr16_pileupfile ./varscan_output/variant

# Check if the files were generated, if not stop the script
if test -f ./varscan_output/variant.indel && test -f ./varscan_output/variant.snp; then
  echo -e "\e[1;32mBoth the output files were created.\e[0m\n"
else
  echo -e "\n\e[1;31mThe files have not been created. Please check what went wrong.\e[0m\n"
  exit
fi

# Checking that the output files are correct by seeing if they contain somatic variants
if grep somatic ./varscan_output/variant.indel; then
    echo -e "\e[1;32mThe file variant.indel is correct.\e[0m"
else
    echo -e "\e[1;31mThe file variant.indel is incorrect.\e[0m\n"
fi

if grep somatic ./varscan_output/variant.snp; then
    echo -e "\e[1;32mThe file variant.snp is correct.\e[0m"
else
    echo -e "\e[1;31mThe file variant.snp is incorrect.\n"
fi

# -- 10/ “Basic” VCF Annotation

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                          Step 7/7 : Basic VCF Annotation                      |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

echo -e "\e[1;96mExtracting all the somatic mutations from the variant list :\e[0m\n"

# Extracting only the somatic variants from the variant list
grep -i 'somatic' ./varscan_output/variant.snp ./varscan_output/variant.indel > ./varscan_output/somatic_variants.vcf
# Creating a file containing only the column of interest
awk '{OFS="\t"; if (!/^#/){split($1,a,":"); print a[2],$2-1,$2,$4"/"$5,"+"}}' ./varscan_output/somatic_variants.vcf > ./varscan_output/somatic_variants_temp.bed
# Removing the two "header" lines
awk '!($1=="chrom")' ./varscan_output/somatic_variants_temp.bed > ./varscan_output/somatic_variants.bed

# Check if the files were generated, if not stop the script
if test -f ./varscan_output/somatic_variants.bed ; then
  echo -e "\e[1;32mYou can find the list of somatic variants in the file : ./varscan_output/somatic_variants.bed\e[0m\n"
else
  echo -e "\n\e[1;31mThe file has not been created. Please check what went wrong.\e[0m\n"
  exit
fi

echo -e "\e[1;32mExtracting the gene annotations for the somatic mutations :\e[0m\n"

# Extracting the gene annotations for the somatic mutations
bedtools intersect -a ./annotation/gencode.v24lift37.basic.annotation.gtf -b ./varscan_output/somatic_variants.bed > ./annotation/intersect-file

# Check if the files were generated, if not stop the script
if test -f ./annotation/intersect-file ; then
  echo -e "\e[1;32mDone \e[0m\n"
else
  echo -e "\n\e[1;31mThe file has not been created. Please check what went wrong.\e[0m\n"
  exit
fi

echo -e "\e[1;96mCreating a list of the genes containing a somatic mutation :\e[0m\n"

# Creating a file containing the genes with a somatic mutation
grep '\sgene\s' ./annotation/intersect-file | awk '{print " " $1 " " $4 " " $5 " " $16}' > ./annotation/variants_genes

# Check if the files were generated, if not stop the script
if test -f ./annotation/intersect-file ; then
  echo -e "\e[1;32mYou can find the list of genes in the file : ./annotation/variants_genes\e[0m\n"
else
  echo -e "\n\e[1;31mThe file has not been created. Please check what went wrong.\e[0m\n"
  exit
fi

echo -e "\e[1;46;30m-----------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;46;30m|                                                                                 |\e[0m"
echo -e "\e[1;46;30m|                             End of the pipeline                                 |\e[0m"
echo -e "\e[1;46;30m|                                                                                 |\e[0m"
echo -e "\e[1;46;30m-----------------------------------------------------------------------------------\e[0m"
