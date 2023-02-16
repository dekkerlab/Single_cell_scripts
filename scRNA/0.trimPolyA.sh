
#!/bin/bash

mkdir zUMI_reads_150

for condition in Fix-c-1127-G1 Fix-c-1127-proM; #add your conditions here ("Fix-c-1127-G1"as an example):
do
  #trim PolyA
  for i in {1..48}; #here can be changed according to how many cells in this condition
  do
    java -jar /home/huox/software/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 \
    ${condition}-${i}_R1_001.fastq.gz ${condition}-${i}_R2_001.fastq.gz \
    ${condition}_${i}_R1_001.paired.fastq.gz ${condition}_${i}_R1_001.unpaired.fastq.gz \
    ${condition}_${i}_R2_001.paired.fastq.gz ${condition}_${i}_R2_001.unpaired.fastq.gz \
    ILLUMINACLIP:/home/huox/software/trimmomatic/Trimmomatic-0.36/adapters/PolyT.fa:2:30:10:8:true; \
    echo 'first trim';
    seqkit seq -m 150 ${condition}_${i}_R1_001.paired.fastq.gz > ${condition}_${i}_R1_001.paired.larger150.fastq ; \
    seqkit seq -m 150 ${condition}_${i}_R2_001.paired.fastq.gz > ${condition}_${i}_R2_001.paired.larger150.fastq ;\
    echo 'seqkit';
    gzip ${condition}_${i}_R1_001.paired.larger150.fastq; \
    gzip ${condition}_${i}_R2_001.paired.larger150.fastq; \
    echo 'gzip';
    java -jar /home/huox/software/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 \
    ${condition}_${i}_R1_001.paired.larger150.fastq.gz ${condition}_${i}_R2_001.paired.larger150.fastq.gz \
    ${condition}_${i}_R1_001.paired2.larger150.fastq.gz ${condition}_${i}_R1_001.unpaired3.fastq.gz \
    ${condition}_${i}_R2_001.paired2.larger150.fastq.gz ${condition}_${i}_R2_001.unpaired3.fastq.gz \
    ILLUMINACLIP:/home/huox/software/trimmomatic/Trimmomatic-0.36/adapters/PolyT.fa:2:30:10:8:true; \
    echo 'second trim';
  done 



  #organize files
  mv *.paired2.larger150.fastq.gz ./zUMI_reads_150


  #remove useless files
  rm *.unpaired3.fastq.gz
  rm *.unpaired.fastq.gz 
  rm *.paired.larger150.fastq.gz
  rm *.paired.fastq.gz

  #rename file into zUMI recognized format
  cd ./zUMI_reads_150
  for i in {1..48}; #here can be changed according how many cells
  do
    mv ${condition}_${i}_R2_001.paired2.larger150.fastq.gz ${condition}_${i}_R2_001.fastq.gz; \
    mv ${condition}_${i}_R1_001.paired2.larger150.fastq.gz ${condition}_${i}_R1_001.fastq.gz; \
   done

done

Rscript /home/huox/software/zUMI/zUMIs/misc/merge_demultiplexed_fastq.R --dir ./


  
