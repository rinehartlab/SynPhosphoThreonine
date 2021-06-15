#! /usr/bin/env python

import sys

def sort_seqs():
    """ Calls bash code to sort the sequences and count the unique sequences for each label and then sort again. When the code is used on the cluster this will be done by a line of bash code in the PBS file.
    """
    import subprocess
    subprocess.call("~/scripts/karl_sort.sh", shell=True, executable='/bin/sh')



def sort_seqs2():
    """ Calls bash code to sort the sequences and count the unique sequences for each label and then sort again. When the code is used on the cluster this will be done by a line of bash code in the PBS file.
    """
    import subprocess
    subprocess.call("~/scripts/karl_sort_2.sh", shell=True, executable='/bin/sh')



def sort_seqs3():
    """ Calls bash code to sort the sequences and count the unique sequences for each label and then sort again. When the code is used on the cluster this will be done by a line of bash code in the PBS file.
    """
    import subprocess
    subprocess.call("~/scripts/karl_sort_3.sh", shell=True, executable='/bin/sh')



def reverse_complement(seq):
    seq_dict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])



def extract_variable_regions(seq):
    """Extracts the variable regions from a read in order to reduce potential for sequencing errors in the fixed regions from leading to misassignment of reads to different clusters
    """
    var_start = 'GGTACCAAG'
    var_end = 'AAGCTT'
    var_seq_start_index = seq.find(var_start)+len(var_start)
    var_seq_end_index = seq.find(var_end)-1
    var_seq = seq[var_seq_start_index:var_seq_end_index]
    return var_seq



def match_barcodes(seq):
    """Searches a sequence for start and end barcode and returns the group associated with the barcode.
    """
    barcode_dict = {1: {'start': 'AGTCTGGGTCGACT', 'stop': 'TACATGGTACGCT', 'start_degenerate': 6, 'stop_degenerate': 9},
                    2: {'start': 'TCTCTGGGTCGACT', 'stop': 'TACATGGTACGGA', 'start_degenerate': 7, 'stop_degenerate': 8},
                    3: {'start': 'AGTCTGGGTCGACT', 'stop': 'TACATGGTACGGA', 'start_degenerate': 6, 'stop_degenerate': 8},
                    4: {'start': 'TCTCTGGGTCGACT', 'stop': 'TACATGGTACGCT', 'start_degenerate': 7, 'stop_degenerate': 9},
                    5: {'start': 'GATCTGGGTCGACT', 'stop': 'TACATGGTACGTC', 'start_degenerate': 8, 'stop_degenerate': 7},
                    6: {'start': 'CTTCTGGGTCGACT', 'stop': 'TACATGGTACGAG', 'start_degenerate': 9, 'stop_degenerate': 6},
                    7: {'start': 'GATCTGGGTCGACT', 'stop': 'TACATGGTACGAG', 'start_degenerate': 8, 'stop_degenerate': 6},
                    8: {'start': 'CTTCTGGGTCGACT', 'stop': 'TACATGGTACGTC', 'start_degenerate': 9, 'stop_degenerate': 7}}
    group = None
    for key, value in barcode_dict.iteritems():
        start = seq.find(value['start'])
        end = seq.find(value['stop'])
        leading = start
        trailing = len(seq) - (end + len(value['stop']))

        if start != -1 and end != -1 and leading == value['start_degenerate'] and trailing == value['stop_degenerate']:
#        if start != -1 and end != -1:
            group = key
            break
    return group, start, end



def reverse_sequence(read):
    """ Takes in a sequence and makes sure its in the correct orientation for downstream analysis then returns the sequence. This is done in order to ensure all reads from sequencing experiments in which the sequencing primers were ligated on are in the same orientation. This uses a unique sequence found in the correct orientation version of the sequence.
    """
    read_key = 'TCTGGGTCGACTGGTG'
    read_key_short = 'GGGTCGACTGGTG'
    if read.find(read_key_short) == -1:
        read = reverse_complement(read)
    return read



def process_reads_kyle(FileName, suffix, barcode):
    """ This code takes in the assembled paired end data from Karl's sequencing experiments. It removes the the variable random bases at the ends. It sorts them into groups based on the barcodes on each end. The categorized reads are then exported as a tsv.

The reads that this process misses as well as the total number of reads missed is exported as a second file.
    """
    out = open('bbseq'+suffix+'.txt','w')
    out2 = open('missedSeqs'+suffix+'.txt','w')
# Keep track of the number of missed reads as well as the missed reads.
    missed = 0
    with open(FileName) as InFile:
# Go through each sequence and match the barcodes in oder to determine which group to which it belongs.
# The default group is set such that if the read doesn't match any group it is exported to a junk file
# Why are the start and stop indices +2 and -1 off the read.find() index?
        for line in InFile:
            read = line.rstrip('\n')
            read = reverse_sequence(read)
#            group, start, end = match_barcodes(read)

#            if group is None:
#                missed += 1
#                out2.write("%s\n" % read)
#            else:
#                trimmed_read = read[start:end]
            trimmed_read = read
            var_read = extract_variable_regions(trimmed_read)
            if barcode:
                export_read = trimmed_read + '\n'
            else:
                export_read = var_read + '\n'
            out.write(export_read)
    out.close()
    out2.write(str(missed)+'\n')
    out2.close()


def process_reads(FileName, suffix, barcode):
    """ This code takes in the assembled paired end data from Karl's sequencing experiments. It removes the the variable random bases at the ends. It sorts them into groups based on the barcodes on each end. The categorized reads are then exported as a tsv.

The reads that this process misses as well as the total number of reads missed is exported as a second file.
    """
    out = open('bbseq'+suffix+'.txt','w')
    out2 = open('missedSeqs'+suffix+'.txt','w')
# Keep track of the number of missed reads as well as the missed reads.
    missed = 0
    with open(FileName) as InFile:
# Go through each sequence and match the barcodes in oder to determine which group to which it belongs.
# The default group is set such that if the read doesn't match any group it is exported to a junk file
# Why are the start and stop indices +2 and -1 off the read.find() index?
        for line in InFile:
            read = line.rstrip('\n')
            read = reverse_sequence(read)
            group, start, end = match_barcodes(read)

            if group is None:
                missed += 1
                out2.write("%s\n" % read)
            else:
                trimmed_read = read[start:end]
                var_read = extract_variable_regions(trimmed_read)
                if barcode:
                    export_read = trimmed_read + '\t' + str(group) + '\n'
                else:
                    export_read = var_read + '\t' + str(group) + '\n'
                out.write(export_read)
    out.close()
    out2.write(str(missed)+'\n')
    out2.close()


def pool_percentage(FileName, suffix):
    """Calculate the percentage of each pool that a given read represents.
    """
## Calculate the total number of reads present in each of the groups
    groups = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0}
    with open(FileName) as InFile:
        for line in InFile:
            clust_list = line.strip('\n').split()
            try:
                cur_group = int(clust_list[2])
            except IndexError:
                continue
            abundance = int(clust_list[0])
            groups[cur_group] += abundance

# Create a new export file that exports the number of reads, relative abundance, group, and sequence as a tsv.
    out3 = open('clustabund'+suffix+'.txt','w')
    with open(FileName) as InFile2:
        for line in InFile2:
            clust_list = line.strip().split()
            try:
                cur_group = int(clust_list[2])
            except IndexError:
                continue
            abundance = int(clust_list[0])
            group_total = float(groups[cur_group])
            rel_abundance = abundance/group_total
            export_string = '%s\t%f\t%s\t%s\n' % (clust_list[0], rel_abundance, clust_list[2], clust_list[1])
            out3.write(export_string)
    out3.close()



def analyze_missed(FileName, suffix):
    """ Analyzes the missed sequences for patterns to determine the most common reasons for failing the filter.
    """
    import re
    seq_start = 'TCTGGGTCGACT'
    seq_end = 'TACATGGTACG'
    match = 0
    start_match = 0
    end_match = 0
    no_match = 0
    start_truncated = 0
    end_truncated = 0

    out = open('debug_data'+suffix+'.txt','w')
    with open(FileName) as InFile:
        for line in InFile:
            line = line.strip('\n')
            start = line.find(seq_start)
            end = line.find(seq_end)
            start_barcode = '--'
            end_barcode = '--'
            if start != -1 and end != -1:
                match += 1
                start_barcode = line[start-2:start]
                if len(start_barcode) < 2:
                    start_truncated += 1
                end_barcode = line[end+len(seq_end):end+len(seq_end)+2]
                if len(end_barcode) < 2:
                    end_truncated += 1
            elif start != -1 and end == -1:
                start_match += 1
                start_barcode = line[start-2:start]
                if len(start_barcode) < 2:
                    start_truncated += 1
            elif end != -1 and start == -1:
                end_match += 1
                end_barcode = line[end+len(seq_end):end+len(seq_end)+2]
                if end_barcode < 2:
                    end_truncated += 1
            else:
                no_match += 1
            export_string_i = '%s\t%s\n' % (start_barcode, end_barcode)
            out.write(export_string_i)

    out.write('#\t#\tmatch\tstart_match\tend_match\tno_match\tstart_truncated\tend_truncated\n')
    export_string = '#\t#\t%d\t%d\t%d\t%d\t%d\t%d\n' % (match, start_match, end_match, no_match, start_truncated, end_truncated)
    out.write(export_string)
    out.close()



def concat_paired_fastas(FileName1, FileName2, suffix):
    """ Processes reads directly from fasta files instead of the assembled amplicons from BBMerge.
    """
    import re
    from itertools import izip
    out = open('concatReads'+suffix+'.txt','w')
    with open(FileName1) as InFile1, open(FileName2) as InFile2:
        for f, r in izip(InFile1, InFile2):
            forward = f.strip('\n')
            forward = forward[0:40]
            reverse = r.strip('\n')
            reverse = reverse[0:40]
            reverse = reverse_complement(reverse)
            read = forward+reverse
            out.write(read+'\n')
    out.close()



def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Prefix for the output processed files')
    parser.add_argument('seqs', help='Text file containing reads')
    parser.add_argument('-p', '--process', action='store_true', help='Process the reads from the start')
    parser.add_argument('-s', '--sort', action='store_true', help='Call a bash script for sorting the reads (not for use on cluster')
    parser.add_argument('-r', '--rel_abund', action='store_true', help='Calculates the relative abundance given a an already sorted and clustered file')
    parser.add_argument('-d', '--debug', action='store_true', help='Analyze the missed sequences file for patterns')
    parser.add_argument('-t', '--test', action='store_true', help='Run the code in test mode')
    parser.add_argument('-f', '--concat_fastas', action='store_true', help='Concatenates paired fasta files for testing purposes')
    parser.add_argument('-b', '--barcode', action='store_true', help='Keep barcodes in the processed reads or remove everything but the variable region')
    args = parser.parse_args()

    seqfile = args.seqs
    suffix = args.filename

    if args.concat_fastas is True:
        concat_paired_fastas('test_r1.fasta','test_r2.fasta',suffix)

    if args.test is True:
        match_barcodes('AGTCTGGGTCGACTTACATGGTACGCT')

    if args.process is True:
        process_reads_kyle(seqfile, suffix, args.barcode)
#        sort_seqs()

    if args.rel_abund is True:
        pool_percentage('bbseq_clust_sorted'+suffix+'.txt', suffix)
#        sort_seqs2()

    if args.debug is True:
        analyze_missed('missedSeqs'+suffix+'.txt', suffix)
#        sort_seqs3()


if __name__ == "__main__":
    main()
