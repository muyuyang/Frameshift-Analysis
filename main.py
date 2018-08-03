
# 1. Align
# 2. Calculate statistics
# 3. Find gaps and determine which could be frameshift
# 4. Get gene information from annotation
# 5. Get protein information from uniprot
# Need to assemble everything I have

# Input: Sequence and annotation of both ancestral and derived genome
# in Genbank format

import sys
from Bio import SeqIO
from alignment import align
from gap_find import gap_find
from evaluate_effect import evaluate_effect
from statistics import statistics
from classes import *
from read_files import *
import time


def convert_gb_to_fasta(input_path,output_folder):
    index = input_path.rfind('/')+1
    dot = input_path.rfind('.')
    output_path = output_folder + input_path[index:dot] + '.fasta'
    SeqIO.convert(input_path,'genbank',output_path,'fasta')
    return output_path


    
def get_whole_sequence_fasta(path):
    seq = ' '
    with open(path,'r') as f:
        lines = f.read().splitlines()
    i = 0
    while not lines[i][0].isalpha():
        i += 1
    for j in range(i,len(lines)):
        seq += lines[j].strip()
    return seq

# Put things together!
def main(genome1_path,genome2_path,output_path):
    if output_path[-1] != '/':
        output_path += '/'
    #Convert genbank to fasta
    fasta1_path = convert_gb_to_fasta(genome1_path,output_path)
    fasta2_path = convert_gb_to_fasta(genome2_path,output_path)
    #Adapt to command line path style
    adapted_path1 = fasta1_path.replace(' ','\ ')
    adapted_path2 = fasta2_path.replace(' ','\ ')
    adapted_output_path = output_path.replace(' ','\ ')
    #Align two genomes and write results under output_path
    align(adapted_path1,adapted_path2,adapted_output_path)
    #Hardcoding:((
    alignment_path = output_path + 'align-2.maf'
    alignments = read_alignment(alignment_path)
    annotation1 = read_annotation(genome1_path)
    annotation2 = read_annotation(genome2_path)


    stat_path = output_path + 'statistics.txt'
    stat = statistics(alignments,annotation1,annotation2,stat_path)
    # Find gaps on CDS and determine which could be frameshift. Get information.
    frameshifts = gap_find(alignments,annotation1,annotation2,output_path)
    originalSeq1 = get_whole_sequence_fasta(fasta1_path)
    originalSeq2 = get_whole_sequence_fasta(fasta2_path)
    feature_path = output_path + 'features.txt'
    evaluate_effect(frameshifts,originalSeq1,originalSeq2,feature_path)


# genome1_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/K-12_annotation.gb'
# genome2_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/O157_H7_annotation.gb'
# output_path = '/Users/muyuyang/Desktop/Frith Lab/TEST/'


genome1_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/data/Oryza sativa.gb'
genome2_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/data/Oryza coarctata.gb'
output_path = '/Users/muyuyang/Desktop/Frith Lab/TEST2/'


main(genome1_path,genome2_path,output_path)

# if __name__ == '__main__':
#     genome1_path = sys.argv[1]
#     genome2_path = sys.argv[2]
#     output_path = sys.argv[3]
    #options...



