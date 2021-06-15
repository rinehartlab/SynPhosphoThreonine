# SynPhosphoThreonine
Files related to the mapping of Thr and pThr phosphosites in a human representative phospho-library

An example dataset is provided for demonstration. 

Data is first seperated into files by deinterleave before be taken forward with analysis. Data is filtered for quality using Trimmomatic, which applies a sliding window filter of width 2 bp and a Phred score cutoff of 30. If the average quality score over two consecutive bases fell below 30, the read is trimmed to remove the remaining bases. Quality trimmed read pairs are then merged using BBMerge with the stringency set to “strict” (https://sourceforge.net/projects/bbmap/). Using the associated scripts, the merged reads are then sorted and assigned to the various input libraries based on barcodes added during the PCR amplification step, libraries are provided in the associated csv files. The variable sequence region for each amplicon was then extracted, and for each input library the abundance of every unique sequence was calculated. The code has been adapted for CRISPR amplicon sequencing through MGH CCIB DNA Core.
