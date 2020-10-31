#!/bin/bash

#$ -l h_rt=336:0:0
#$ -cwd
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -N cirseq

set -e

WORKDIR=$1
REFDIR_IDX=$2
THREAD=$3
SCRIPTDIR=$4
SPLICESITE=$5
FASTQS="${@:6}"

#check if working directory exists
if [ ! -d "$WORKDIR" ];then mkdir $WORKDIR; fi
if [ ! -d "$SCRIPTDIR" ];then 
	echo "Directory for scripts does not exist!"
	exit 0
fi

# generate consensus
python ${SCRIPTDIR}/ConsensusGeneration.py $WORKDIR $FASTQS

# align consensus
#bowtie2 -q --phred33 --no-unal --no-hd --no-sq --local -x ${REFIDX} -U ${WORKDIR}/4_rearranged.fastq.gz | gzip -c > ${WORKDIR}/6_alignment.sam.gz
hisat2 -p $THREAD -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/1_consensus.fastq.gz | samtools view -@ $THREAD -b -o ${WORKDIR}/2_alignment.bam

# preprocess 1
python ${SCRIPTDIR}/preprocessing_1.py ${WORKDIR}

# align again
#bowtie2 -q --phred33 --no-unal --no-hd --no-sq --local -x ${REFIDX} -U ${WORKDIR}/4_rearranged.fastq.gz | gzip -c > ${WORKDIR}/6_alignment.sam.gz
hisat2 -p $THREAD --no-unal -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/4_rearranged.fastq.gz | samtools view -@ $THREAD -b -o ${WORKDIR}/6_alignment.bam

# preprocess 2
python ${SCRIPTDIR}/preprocessing_2.py ${WORKDIR} 

#bowtie2 -q --phred33 --no-unal --no-hd --no-sq --local -x ${REFIDX} -U ${WORKDIR}/5_rotated.fastq.gz | gzip -c > ${WORKDIR}/9_alignment.sam.gz
hisat2 -p $THREAD --no-unal -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/5_rotated.fastq.gz | samtools view -@ $THREAD -b -o ${WORKDIR}/9_alignment.bam

#bowtie2 -q --phred33 --no-unal --no-hd --no-sq --local -x ${REFIDX} -U ${WORKDIR}/8_rotated.fastq.gz | gzip -c > ${WORKDIR}/10_alignment.sam.gz
hisat2 -p $THREAD --no-unal -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/8_rotated.fastq.gz | samtools view -@ $THREAD -b -o ${WORKDIR}/10_alignment.bam

# preprocess 3
python ${SCRIPTDIR}/preprocessing_3.py ${WORKDIR} 

samtools merge -@ $THREAD ${WORKDIR}/consensus_align.bam ${WORKDIR}/3_alignment.bam ${WORKDIR}/7_alignment.bam ${WORKDIR}/11_alignment.bam
#rm -f ${WORKDIR}/3_alignment.bam ${WORKDIR}/7_alignment.bam ${WORKDIR}/11_alignment.bam ${WORKDIR}/2_alignment.bam ${WORKDIR}/4_rearranged.fastq.gz ${WORKDIR}/5_rotated.fastq.gz ${WORKDIR}/6_alignment.bam ${WORKDIR}/8_rotated.fastq.gz ${WORKDIR}/9_alignment.bam ${WORKDIR}/10_alignment.bam

#Scripts below are specific for analysing poliovirus genome. Don't apply it to human or other genome!!!
#python ${SCRIPTDIR}/QualityFilter.py ${WORKDIR} ${REFFILE}
#Rscript ${SCRIPTDIR}/ParameterPlots.R ${WORKDIR}
