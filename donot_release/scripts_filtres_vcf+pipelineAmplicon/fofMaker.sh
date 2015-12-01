#!/bin/bash

#bash fofMaker.sh  <pathToIncaFile(fq)> <pathToFofFile>

#Mandatory :
#575C-a1_R1_Q30.fastq.gz
#575C-a1_R2_Q30.fastq.gz
#575C-a2_R1_Q30.fastq.gz
#/575C-a2_R2_Q30.fastq.gz
pathToIncaFile=$1
pathToFofFile=$2

for fichier in $(ls $pathToIncaFile|cut -f1 -d"-"|sort|uniq);do
        nameFile=$(echo $fichier|cut -f1 -d"-")
        touch fof$nameFile\-a1-R1_R2_Q30.txt
        touch fof$nameFile\-a2-R1_R2_Q30.txt
        touch fof$nameFile\.txt
        echo $pathToIncaFile$nameFile\-a1_R1_Q30.fastq.gz >> $pathToFofFile\fof$nameFile\-a1-R1_R2_Q30.txt
        echo $pathToIncaFile$nameFile\-a1_R2_Q30.fastq.gz >> $pathToFofFile\fof$nameFile\-a1-R1_R2_Q30.txt
        echo $pathToIncaFile$nameFile\-a2_R1_Q30.fastq.gz >>$pathToFofFile\fof$nameFile\-a2-R1_R2_Q30.txt
        echo $pathToIncaFile$nameFile\-a2_R2_Q30.fastq.gz >> $pathToFofFile\fof$nameFile\-a2-R1_R2_Q30.txt
        echo $pathToFofFile\fof$nameFile\-a1-R1_R2_Q30.txt >> $pathToFofFile\fof$nameFile\.txt
        echo $pathToFofFile\fof$nameFile\-a2-R1_R2_Q30.txt >> $pathToFofFile\fof$nameFile\.txt
        echo "done with $nameFile"
done        
