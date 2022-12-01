# -- Script du pipeline d'analyse RNAseq
# TP 1 - RNA-seq Pipelines
# Soundous Bella Baci

#!/bin/bash

cd TP_rnaseq_SBB
echo -e "\e[1;96mThe directory for this session is called : ./TP_rnaseq_SBB  \e[0m\n"


echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                  Step 1/ : Quality control and trimming                       |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"


sample_list=(Day_0_1_chr18.sampled Day_0_2_chr18.sampled Day_0_3_chr18.sampled Day_7_1_chr18.sampled Day_7_2_chr18.sampled Day_7_3_chr18.sampled)

echo -e "\e[1;96mCreating fastQC files for all the samples : \e[0m\n"

for i in ${sample_list[*]}
do
  # Check if file exist, if not create them
  if [ ! -d output_fastQC_nontrimmed ];then
    mkdir output_fastQC_nontrimmed
    echo -e "\e[1;96mCreating a directory : ./output_fastQC_nontrimmed \e[0m\n"
  fi

  # -- Etape de contrôle qualité : génération de fichiers fastQC pour chaque fichier fastq
  fastqc -o output_fastQC_nontrimmed ./fastq_files/$i'.R1.fastq' ./fastq_files/$i'.R2.fastq'

  # Check if the files were generated, if not stop the script
  if test -f ./output_fastQC_nontrimmed/$i'.R1_fastqc.html' && test -f ./output_fastQC_nontrimmed/$i'.R2_fastqc.html'; then
    echo -e "\n\e[1;32mFastQC files were created for the file $i. You can find them in the directory : output_fastQC_nontrimmed\e[0m\n"
  else
    echo -e "\n\e[1;31mSomething went wrong with the fastQC command. Please check what went wrong.\e[0m\n"
    exit
  fi
done

echo -e "\e[1;96mEliminating low quality sequences at read extermities : \e[0m\n"

for i in ${sample_list[*]}
do
  # Check if file exist, if not create them
  if [ ! -d trimmed_fastq ];then
    mkdir trimmed_fastq
    echo -e "\e[1;96mCreating a directory : ./trimmed_fastq \e[0m\n"
  fi

  ## -- Elimination des séquences de basses qualité à l'extrémité finale des reads
  trimmomatic PE ./fastq_files/$i'.R1.fastq' ./fastq_files/$i'.R2.fastq' -baseout ./trimmed_fastq/$i'_trimmed.fastq'  LEADING:20 TRAILING:20 MINLEN:50

  # Check if the files were generated, if not stop the script
  if test -f ./trimmed_fastq/$i'_trimmed_1P.fastq' && test -f ./trimmed_fastq/$i'_trimmed_2P.fastq' && test -f ./trimmed_fastq/$i'_trimmed_1U.fastq' && test -f ./trimmed_fastq/$i'_trimmed_2U.fastq'; then
    echo -e "\n\e[1;32mTrimmed fastq created for the file $i. You can find them in the directory : trimmed_fastq\e[0m\n"
  else
    echo -e "\n\e[1;31mSomething went wrong with the trimmomatic command. Please check what went wrong.\e[0m\n"
    exit
  fi
done

echo -e "\e[1;96mCreating fastQC files for all the paired trimmed samples : \e[0m\n"

for i in ${sample_list[*]}
do
  if [ ! -d output_fastQC_trimmed ];then #
    mkdir output_fastQC_trimmed
    echo -e "\e[1;96mCreating a directory : ./output_fastQC_trimmed \e[0m\n"
  fi

  # -- Deuxième étape de contrôle qualité : génération de fichiers fastQC pour les fichiers fastq trimmés
  fastqc -o output_fastQC_trimmed ./trimmed_fastq/$i'_trimmed_1P.fastq' ./trimmed_fastq/$i'_trimmed_2P.fastq'

  # Check if the files were generated, if not stop the script
  if test -f ./output_fastQC_trimmed/$i'_trimmed_1P_fastqc.html' && test -f ./output_fastQC_trimmed/$i'_trimmed_2P_fastqc.html' ; then
    echo -e "\n\e[1;32mFastQC files were created for the file $i. You can find them in the directory : output_fastQC_trimmed\e[0m\n"
  else
    echo -e "\n\e[1;31mSomething went wrong with the fastQC command. Please check what went wrong.\e[0m\n"
    exit
  fi
done

echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m|                              Step 2/ : Mapping                                |\e[0m"
echo -e "\e[1;47;30m|                                                                               |\e[0m"
echo -e "\e[1;47;30m---------------------------------------------------------------------------------\e[0m"

# -- Indexation du génome avec STAR seulement si le fichier n'existe pas déjà

echo -e "\e[1;96mIndexing the genome  :\e[0m\n"
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ./indexgenome --genomeFastaFiles ./genome/chr18.fa --sjdbGTFfile ./annotation/gencode.v24lift37.basic.annotation.gtf --genomeSAindexNbases 12
echo -e "\n\e[1;32mYou can find the indexed genome in the directory : ./indexgenome\e[0m\n"

for i in ${sample_list[*]}
do
  STAR --runThreadN 8 --outFilterMultimapNmax 1 --genomeDir ./mapping --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./star/$i'_mapping' --readFilesIn ./trimmed_fastq/$i'_trimmed_1P.fastq' ./trimmed_fastq/$i'_trimmed_2P.fastq'
done

echo -e "\nTerminé"
echo -e "\e[1;96mMapping  :\e[0m\n"
for i in ${sample_list[*]}
do
  if [ ! -d mapping ];then
    mkdir mapping
  fi
  # -- Mapping avec STAR
  STAR --runThreadN 8 --outFilterMultimapNmax 1 --genomeDir mapping --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./mapping/$i'_mapping' --readFilesIn ./trimmed_fastq/$i'_trimmed_1P.fastq' ./trimmed_fastq/$i'_trimmed_2P.fastq'
  # -- Tri et indexation du fichier BAM
  samtools index -b ./mapping/$i'_mappingAligned.sortedByCoord.out.bam' #output .bai
done


# -- Création des dossiers pour organiser les résultats

if [ ! -d result_mapping ];then
  mkdir result_mapping
fi
# ---------------------------------


for i in ${sample_list[*]}
do
  # -- Mapping avec STAR
  STAR --runThreadN 8 --outFilterMultimapNmax 1 --genomeDir indexgenome --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./result_mapping/$i'_res_mapping' --readFilesIn ./trimmed_fastq/$i'_trimmed_1P.fastq' ./trimmed_fastq/$i'_trimmed_2P.fastq'
  # -- Tri et indexation du fichier BAM
  samtools index -b ./result_mapping/$i'_res_mappingAligned.sortedByCoord.out.bam' #output .bai
done

# -- Compatage des reads pour tous les résultats de l'alignement
# -- pas de boucle pour avoir un seul fichier de sortie
featureCounts -p -t exon -g gene_id -a gencode.v24lift37.basic.annotation.gtf -o ./output_comptage.txt ./result_mapping/*_res_mappingAligned.sortedByCoord.out.bam

# -- Création d'une table d'équivalence entre les noms ENCODE et les noms HUGO
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' gencode.v24lift37.basic.annotation.gtf | sort | uniq > encode-to-hugo.tab

# -- Tri du fichier de sortie du compate et de la table d'équivalence
sort output_comptage.txt > temp1
sort encode-to-hugo.tab > temp2

# -- Création d'un fichier trié sur le chromosome 18 à partir des deux fichiers triés précédents
join temp1 temp2 |grep "chr18" > temp3

# -- Création de la table de comptage en récupérant les colonnes d'intérêt dans le fichier précédent
awk '{print $13 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12}' temp3 > table_count
