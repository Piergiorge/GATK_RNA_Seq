#!/bin/bash

# Definir as variáveis de ambiente
export PATH=/usr/local/bin:$PATH
export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
export GATK=/home/user/gatk-4.2.2.0/gatk
export PICARD=/home/user/picard.jar
export SAMTOOLS=/home/user/samtools-1.13/samtools
export BWA=/home/user/bwa-0.7.17/bwa
export FASTQC=/home/user/FastQC/fastqc
export TRIMMOMATIC=/home/user/Trimmomatic-0.39/trimmomatic-0.39.jar
export MULTIQC=/home/user/miniconda3/bin/multiqc

# Definir os arquivos de referência
REF=/home/user/reference/GRCh38.fa # Genoma de referência humano GRCh38
GTF=/home/user/reference/gencode.v38.annotation.gtf # Arquivo de anotação gênica humano GENCODE v38
DBSNP=/home/user/reference/dbsnp_153.b38.vcf.gz # Banco de dados de variantes humano dbSNP 153

# Criar um diretório para armazenar os resultados
mkdir -p results

# Loop para processar cada amostra
for SAMPLE in $(cat samples.txt)
do
  # Criar um diretório para cada amostra
  mkdir -p results/$SAMPLE
  
  # Fazer o controle de qualidade das reads
  $FASTQC -o results/$SAMPLE data/$SAMPLE\_R1.fastq.gz data/$SAMPLE\_R2.fastq.gz
  
  # Cortar as reads de baixa qualidade e remover os adaptadores
  java -jar $TRIMMOMATIC PE -threads 4 -phred33 data/$SAMPLE\_R1.fastq.gz data/$SAMPLE\_R2.fastq.gz \
  results/$SAMPLE/$SAMPLE\_R1.trimmed.fastq.gz results/$SAMPLE/$SAMPLE\_R1.unpaired.fastq.gz \
  results/$SAMPLE/$SAMPLE\_R2.trimmed.fastq.gz results/$SAMPLE/$SAMPLE\_R2.unpaired.fastq.gz \
  ILLUMINACLIP:/home/user/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  
  # Alinhar as reads ao genoma de referência usando BWA
  $BWA mem -t 4 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" $REF results/$SAMPLE/$SAMPLE\_R1.trimmed.fastq.gz results/$SAMPLE/$SAMPLE\_R2.trimmed.fastq.gz | \
  $SAMTOOLS view -bS - > results/$SAMPLE/$SAMPLE.bam
  
  # Ordenar e indexar o arquivo BAM
  java -jar $PICARD SortSam I=results/$SAMPLE/$SAMPLE.bam O=results/$SAMPLE/$SAMPLE.sorted.bam SORT_ORDER=coordinate
  $SAMTOOLS index results/$SAMPLE/$SAMPLE.sorted.bam
  
  # Marcar e remover os duplicados
  java -jar $PICARD MarkDuplicates I=results/$SAMPLE/$SAMPLE.sorted.bam O=results/$SAMPLE/$SAMPLE.dedup.bam M=results/$SAMPLE/$SAMPLE.metrics.txt REMOVE_DUPLICATES=true
  $SAMTOOLS index results/$SAMPLE/$SAMPLE.dedup.bam
  
  # Realinhar os indels
  $GATK --java-options "-Xmx4g" IndelRealigner -R $REF -I results/$SAMPLE/$SAMPLE.dedup.bam -targetIntervals results/$SAMPLE/$SAMPLE.intervals -o results/$SAMPLE/$SAMPLE.realign.bam
  
  # Recalibrar a qualidade das bases
  $GATK --java-options "-Xmx4g" BaseRecalibrator -R $REF -I results/$SAMPLE/$SAMPLE.realign.bam -knownSites $DBSNP -o results/$SAMPLE/$SAMPLE.recal.table
  $GATK --java-options "-Xmx4g" PrintReads -R $REF -I results/$SAMPLE/$SAMPLE.realign.bam -BQSR results/$SAMPLE/$SAMPLE.recal.table -o results/$SAMPLE/$SAMPLE.recal.bam
  
  # Chamar os variantes usando o HaplotypeCaller
  $GATK --java-options "-Xmx4g" HaplotypeCaller -R $REF -I results/$SAMPLE/$SAMPLE.recal.bam -O results/$SAMPLE/$SAMPLE.vcf -ERC GVCF --dbsnp $DBSNP
  
done

# Gerar um relatório de qualidade das amostras usando MultiQC
$MULTIQC -o results/multiqc_report results
