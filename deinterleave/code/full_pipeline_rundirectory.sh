#!/bin/bash
#SBATCH -J nk_fp
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
SAMPLES=$(ls ../input/ | grep CA)

echo $SAMPLES

for folder in $SAMPLES
do
    cd $code_dir
    SAMPLENAME=${folder#*_}
    mkdir ../temp/$folder
    cd ../temp/$folder
    ../../code/deinterleave_fastq.sh < ../../input/$folder/*.fastq ${SAMPLENAME}_r1.fastq ${SAMPLENAME}_r2.fastq

    cd $code_dir
    mkdir ../output/trimmed_${SAMPLENAME}
    mkdir ../output/merged_${SAMPLENAME}
    mkdir ../output/aligned_${SAMPLENAME}
    mkdir ../output/extracted_${SAMPLENAME}
    mkdir ../output/clustered_${SAMPLENAME}

    for qual in 20 25 30
    do
        cd $code_dir
        mkdir ../output/trimmed_${SAMPLENAME}/qual_${qual}
	mkdir ../output/merged_${SAMPLENAME}/qual_${qual}
	mkdir ../output/aligned_${SAMPLENAME}/qual_${qual}
   	mkdir ../output/extracted_${SAMPLENAME}/qual_${qual}
        mkdir ../output/clustered_${SAMPLENAME}/qual_${qual}	

	cd ../output/trimmed_${SAMPLENAME}/qual_${qual}
        /usr/bin/TrimmomaticPE -threads 8 -phred33 ../../../temp/$folder/${SAMPLENAME}_r1.fastq ../../../temp/$folder/${SAMPLENAME}_r2.fastq ${SAMPLENAME}_r1_trim_paired.fastq.gz ${SAMPLENAME}_r1_trim_unpaired.fastq.gz ${SAMPLENAME}_r2_trim_paired.fastq.gz ${SAMPLENAME}_r2_trim_unpaired.fastq.gz SLIDINGWINDOW:2:${qual}

	cd $code_dir
	cd ../output/merged_${SAMPLENAME}/qual_${qual}	
	/mnt/c/Users/jmoen/Documents/Ubuntu/bbmap/bbmerge.sh in1=../../trimmed_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_r1_trim_paired.fastq.gz in2=../../trimmed_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_r2_trim_paired.fastq.gz out=${SAMPLENAME}_bbmerge_trim.fastq.gz outu=unmerged_${SAMPLENAME}_bbmerge_trim.fastq.gz outinsert=insert_length.txt ihist=len_hist.txt strict=t

	cd $code_dir
	cd ../output/aligned_${SAMPLENAME}/qual_${qual}

	bwa mem -t 8 -M ../../../input/library_fasta/seq_list_ref.fasta ../../merged_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbmerge_trim.fastq.gz > aligned.sam
	samtools view -bS aligned.sam | samtools sort -o aligned_sorted.bam
	samtools index aligned_sorted.bam aligned_sorted.bai
	pileup.sh in=aligned.sam out=PileupStats.txt secondary=false
	samtools idxstats aligned_sorted.bam | cut -f 1,3 > AlignmentStats.txt

	cd $code_dir
	cd ../output/extracted_${SAMPLENAME}/qual_${qual}
        cp ../../merged_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbmerge_trim.fastq.gz .
        gunzip ${SAMPLENAME}_bbmerge_trim.fastq.gz
        cat ${SAMPLENAME}_bbmerge_trim.fastq | awk '{if(NR%4==2) print $1}' > ${SAMPLENAME}_bbseq_unrev.txt

        cd $code_dir
        cd ../output/clustered_${SAMPLENAME}/qual_${qual}/
        python2 ../../../code/demultiplex_and_rc_v5.py _${qual}_${SAMPLENAME} ../../extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -j
        sort bbseq_${qual}_${SAMPLENAME}.txt -T ../../../temp | uniq -c | sort -r -n > bbseq_clust_${qual}_${SAMPLENAME}.txt
#        sort -k3,3n -k1,1nr bbseq_clust_${qual}_${SAMPLENAME}.txt -T ../../../../temp > bbseq_clust_sorted_${qual}_${SAMPLENAME}.txt
#        python2 ../../../../code/demultiplex_and_rc_v4.py _${qual}_${SAMPLENAME} ../../../extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -r
#        sort -k3,3n -k1,1nr clustabund_${qual}_${SAMPLENAME}.txt -T ../../../../temp > sorted_clust_w_abund_${qual}_${SAMPLENAME}.txt
#        python2 ../../../../code/demultiplex_and_rc_v4.py _${qual}_${SAMPLENAME} ../../../extracted_${SAMPLENAME}/qual_${qual}/${SAMPLENAME}_bbseq_unrev.txt -d
#        sort  -k1,2nr debug_data_${qual}_${SAMPLENAME}.txt -T ../../../../temp | uniq -c | sort -r -n > sorted_debug_data_${qual}_${SAMPLENAME}.txt
	
    done
done
