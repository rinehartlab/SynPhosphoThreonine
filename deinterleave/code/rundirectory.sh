#!/bin/bash
#SBATCH -J km_dei
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
    mkdir ../output/$folder
    cd ../output/$folder
    ../../code/deinterleave_fastq.sh < ../../input/$folder/*.fastq ${SAMPLENAME}_r1.fastq ${SAMPLENAME}_r2.fastq
done
