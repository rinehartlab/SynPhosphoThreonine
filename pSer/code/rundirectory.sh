#!/bin/bash
#SBATCH -J km_fe
#SBATCH -p pi_gerstein
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR
#SBATCH --mail-user=paul.muir@yale.edu
#SBATCH --mail-type=ALL


rm -r ../temp/*
rm -r ../output/*

code_dir="$(pwd)"
SAMPLES=$(ls ../input/)


#### PREPROCESSING ####


for folder in $SAMPLES
do
#    cd $code_dir
    SAMPLENAME=${folder#*_}
#    mkdir ../temp/$folder
#    cd ../temp/$folder
#    cat ../../input/$folder/*R1*.gz > ${SAMPLENAME}_r1.fastq.gz
#    cat ../../input/$folder/*R2*.gz > ${SAMPLENAME}_r2.fastq.gz

    cd $code_dir
    mkdir -p ../temp/trimmed_reads/trimmed_${SAMPLENAME}

    for qual in 30 25 20 # Add back 20 and 25 later if Karl wants
    do
        cd $code_dir
        mkdir ../temp/trimmed_reads/trimmed_${SAMPLENAME}/qual_${qual}
        cd ../temp/trimmed_reads/trimmed_${SAMPLENAME}/qual_${qual}
        /usr/bin/TrimmomaticPE -threads 8 -phred33 ../../../../input/$folder/${SAMPLENAME}_r1.fastq ../../../../input/$folder/${SAMPLENAME}_r2.fastq ${SAMPLENAME}_r1_trim_paired.fastq ${SAMPLENAME}_r1_trim_unpaired.fastq ${SAMPLENAME}_r2_trim_paired.fastq ${SAMPLENAME}_r2_trim_unpaired.fastq HEADCROP:3 SLIDINGWINDOW:2:${qual}
    done
done


#### MERGING ####


cd $code_dir
SAMPLES=$(ls ../temp/trimmed_reads/)

for folder in $SAMPLES
do
    cd $code_dir
    SAMPLENAME=${folder#*_}
    mkdir -p  ../temp/merged_reads/merged_${SAMPLENAME}
    for qual in 30 25 20 # Add back 20 and 25 later if Karl wants
    do
      	cd $code_dir
        mkdir ../temp/merged_reads/merged_${SAMPLENAME}/qual_${qual}
        cd ../temp/merged_reads/merged_${SAMPLENAME}/qual_${qual}
#        gunzip ${SAMPLENAME}_r1_trim_paired.fastq.gz
#        gunzip ${SAMPLENAME}_r2_trim_paired.fastq.gz

        /mnt/c/Users/jmoen/Documents/Ubuntu/bbmap/bbmerge.sh in1=../../../../temp/trimmed_reads/trimmed_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_r1_trim_paired.fastq in2=../../../../temp/trimmed_reads/trimmed_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_r2_trim_paired.fastq out=${SAMPLENAME}_bbmerge_trim.fastq outu=unmerged_${SAMPLENAME}_bbmerge_trim.fastq outinsert=insert_length.txt ihist=len_hist.txt strict=t
    done
done


#### ALIGN READS ####
cd $code_dir
SAMPELS=$(ls ../temp/merged_reads/)

for folder in $SAMPLES
do
    cd $code_dir
    SAMPLENAME=${folder#*_}
    mkdir -p  ../output/aligned_${SAMPLENAME}
    for qual in 30 25 20 # Add back 20 and 25 later if Karl wants
    do
	cd $code_dir
	mkdir -p ../output/aligned_${SAMPLENAME}/qual_${qual}
	cd ../output/aligned_${SAMPLENAME}/qual_${qual}
        bwa mem -t 8 -M ../../../code/library_fasta/seq_list_ref.fasta ../../../temp/merged_reads/merged_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbmerge_trim.fastq > aligned.sam
        samtools view -S -bh aligned.sam > aligned.bam
		samtools sort -o aligned.sam > aligned_sorted.bam
        samtools index aligned_sorted.bam aligned_sorted.bai
        /mnt/c/Users/jmoen/Documents/Ubuntu/bbmap/pileup.sh in=aligned.sam out=PileupStats.txt secondary=false
        samtools idxstats aligned_sorted.bam | cut -f 1,3 > AlignmentStats.txt
	done
done




#### EXTRACT READS ####


cd $code_dir
SAMPLES=$(ls ../temp/merged_reads/)

for folder in $SAMPLES
do
    cd $code_dir
    SAMPLENAME=${folder#*_}
    mkdir -p ../temp/extracted_reads/extracted_${SAMPLENAME}
    for qual in 30 25 20 # Can add back 20 and 25 if Karl wants later
    do
      	cd $code_dir
        mkdir ../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}
        cd ../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}
        mkdir -p ../../../${SAMPLENAME}/qual_${qual}
        cp ../../../merged_reads/merged_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbmerge_trim.fastq ../../../${SAMPLENAME}/qual_${qual}/
#	gunzip ../../../${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbmerge_trim.fastq
        cat ../../../${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbmerge_trim.fastq | awk '{if(NR%4==2) print $1}' > ${SAMPLENAME}_bbseq_unrev.txt
    done
done


#### CLUSTER SORT  ####


cd $code_dir
SAMPLES=$(ls ../temp/extracted_reads)

for folder in $SAMPLES
do
    cd $code_dir
    SAMPLENAME=${folder#*_}
    mkdir ../output/clustered_${SAMPLENAME}
    for qual in 30 25 20 # Add back 20 and 25 if needed later
    do
      	cd $code_dir
        mkdir ../output/clustered_${SAMPLENAME}/qual_${qual}
        mkdir ../output/clustered_${SAMPLENAME}/qual_${qual}/wbarcode
        mkdir ../output/clustered_${SAMPLENAME}/qual_${qual}/wobarcode

        cd ../output/clustered_${SAMPLENAME}/qual_${qual}/wbarcode
        echo ${SAMPLENAME}_bbseq_unrev.txt
        python2 ../../../../code/demultiplex_and_rc_v6_kyle.py _${qual}_${SAMPLENAME} ../../../../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -p -b
        sort bbseq_${qual}_${SAMPLENAME}.txt -T ../../../../temp | uniq -c | sort -r -n > bbseq_clust_${qual}_${SAMPLENAME}.txt
        sort -k3,3n -k1,1nr bbseq_clust_${qual}_${SAMPLENAME}.txt -T ../../../../temp > bbseq_clust_sorted_${qual}_${SAMPLENAME}.txt
        python2 ../../../../code/demultiplex_and_rc_v6_kyle.py _${qual}_${SAMPLENAME} ../../../../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -r
        sort -k3,3n -k1,1nr clustabund_${qual}_${SAMPLENAME}.txt -T ../../../../temp > sorted_clust_w_abund_${qual}_${SAMPLENAME}.txt
        python2 ../../../../code/demultiplex_and_rc_v6_kyle.py _${qual}_${SAMPLENAME} ../../../../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -d
        sort  -k1,2nr debug_data_${qual}_${SAMPLENAME}.txt -T ../../../../temp | uniq -c | sort -r -n > sorted_debug_data_${qual}_${SAMPLENAME}.txt

        cd $code_dir
        cd ../output/clustered_${SAMPLENAME}/qual_${qual}/wobarcode
        python2 ../../../../code/demultiplex_and_rc_v6_kyle.py _${qual}_${SAMPLENAME} ../../../../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -p
        sort bbseq_${qual}_${SAMPLENAME}.txt -T ../../../../temp | uniq -c | sort -r -n > bbseq_clust_${qual}_${SAMPLENAME}.txt
        sort -k3,3n -k1,1nr bbseq_clust_${qual}_${SAMPLENAME}.txt -T ../../../../temp > bbseq_clust_sorted_${qual}_${SAMPLENAME}.txt
        python2 ../../../../code/demultiplex_and_rc_v6_kyle.py _${qual}_${SAMPLENAME} ../../../../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -r
        sort -k3,3n -k1,1nr clustabund_${qual}_${SAMPLENAME}.txt -T ../../../../temp > sorted_clust_w_abund_${qual}_${SAMPLENAME}.txt
        python2 ../../../../code/demultiplex_and_rc_v6_kyle.py _${qual}_${SAMPLENAME} ../../../../temp/extracted_reads/extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -d
        sort  -k1,2nr debug_data_${qual}_${SAMPLENAME}.txt -T ../../../../temp | uniq -c | sort -r -n > sorted_debug_data_${qual}_${SAMPLENAME}.txt
    done
done
