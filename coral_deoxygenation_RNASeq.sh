#!/bin/bash

#########################################
#########  Trimming adaptors  ##########
#########################################
for i in *_R1_001.fastq.gz ; do java -jar /home/cardena/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 16 $i ${i//_R1_/_R2_} -baseout ${i//_R1_001.fastq.gz/.fastq} ILLUMINACLIP:/home/cardena/software/Trimmomatic-0.38/adapters/Truseq_V.3_edited.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50; done
########################################################
#########  Mapping To A. millepora gene set   ##########
########################################################
bowtie2-build Amil.all.maker.transcripts.fasta  Amil_DB ## available from https://przeworskilab.com/data/
for i in *_1P.fastq ; do bowtie2 -p 14 --very-sensitive-local -x Amil_DB -1 $i -2 ${i//1P.fastq/2P.fastq} -S $(basename $i | sed  's/_1P.fastq/.sam/') ; done

###############################################
#########  Quantifying transcripts   ##########
###############################################
for i in *R*0-*.sam ; do /home/cardena/software/express-1.5.1-linux_x86_64/express -o $(basename $i | sed 's/M_19_[0][34][0-9][0-9]_//' | sed 's/_D[0-9][0-9][0-9]-D[0-9][0-9][0-9]_L00[0-8]//' | sed 's/.\///'| sed 's/.sam/_express/') Trinity_A.tenuis.fasta $i ; done
#count table
for f in $(ls ./*/results.xprs) ; do cat $f | cut -f2,8 > ./counts/"$(dirname "$f" | sed 's/M_19_[0][34][0-9][0-9]_//' | sed 's/_D[0-9][0-9][0-9]-D[0-9][0-9][0-9]_L00[0-8]//' | sed 's/.\///'| sed 's/_express/_counts/')" ; done
python /home/cardena/scripts/MergingMultipleTables.py . .
#fpkm table
for f in $(ls ./*/results.xprs) ; do cat $f | cut -f2,11 > ./fpkm/"$(dirname "$f" | sed 's/M_19_[0][34][0-9][0-9]_//' | sed 's/_D[0-9][0-9][0-9]-D[0-9][0-9][0-9]_L00[0-8]//' | sed 's/.\///'| sed 's/_express/_fpkm/')" ; done
python /home/cardena/scripts/MergingMultipleTables.py . .
