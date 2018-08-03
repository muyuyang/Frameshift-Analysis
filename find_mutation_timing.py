from ete3 import Tree
import os
import sys
from Bio import SeqIO
from alignment import align
from gap_find import gap_find
from classes import *
from read_files import *

# Plan:
# 1. Read in a phylogenetic tree in Newick format(Biopython)
# 2. Add features needed(?)
# 3. Decide which species to compare
# 4. Determine where the mutation happens -> between which two species

def read_tree(path):
    tree = Tree(path)
    return tree

def find_common_mutation(frameshifts1,frameshifts2):
    info1 = set([(F.start1,F.length,F.seq) for F in frameshifts1])
    info2 = set([(F.start1,F.length,F.seq) for F in frameshifts2])
    common = list(info1.intersection(info2))
    newCommon = [(start1,length,'Insertion') if seq == 1 else (start1,length,'Deletion') for (start1,length,seq) in common  ]
    return newCommon


def ancestry(frameshifts1,frameshifts2):
    common = find_common_mutation(frameshifts1,frameshifts2) 
    # print(frameshifts1)
    # print(frameshifts2)
    common.sort()
    print(common)
    return common


def convert_gb_to_fasta(input_path,output_folder):
    index = input_path.rfind('/')+1
    dot = input_path.rfind('.')
    output_path = output_folder + input_path[index:dot] + '.fasta'
    if not os.path.exists(output_path):
        SeqIO.convert(input_path,'genbank',output_path,'fasta')
    return output_path


def pairwise(genome1_path,genome2_path,output_path):
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
    # Find gaps on CDS and determine which could be frameshift. Get information.
    frameshifts = gap_find(alignments,annotation1,annotation2,output_path)
    return frameshifts


def write_results(path,names,common):
    with open(path,'a') as f:
        f.write(str(names)+'\n')
        f.write(str(common)+'\n\n')


def find_mutation_timing(data_path,output_path,tree_path):
    tree = read_tree(tree_path)
    print(tree)
    # A list of leaves
    species = [node for node in tree.get_leaves()]
    names = [node.name for node in species]
    write_path = output_path + 'common.txt'
    with open(write_path,'w') as f:
        f.write('')
    # TODO: Get sequence data, make alignment and find all frameshifts

    # TODO: Choose A, B and C
    # Current thought: B is the closest possible leaf to A. C is an outgroup 
    # but should still be close.
    # Find all shared mutations between AB and AC to get all frameshifts that 
    # happen during common ancestor of A and B -> A

    # frameshifts = dict() # Store frameshifts that have been found

    for A in species:
        sorted_leaves = sorted(species,
            key=lambda leaf:leaf.get_distance(A,topology_only=True))
        B = sorted_leaves[1]
        C = sorted_leaves[3]
        print(A.name,B.name,C.name)
        genome1_path = data_path + A.name + '.gb'
        genome2_path = data_path + B.name + '.gb'
        genome3_path = data_path + C.name + '.gb'
        frameshifts1 = pairwise(genome1_path,genome2_path,output_path)
        frameshifts2 = pairwise(genome1_path,genome3_path,output_path)
        common = ancestry(frameshifts1,frameshifts2)
        write_results(write_path,(A.name,B.name,C.name),common)



data_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/data/'
output_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/results/'
tree_path = 'oryza_tree.nw'
find_mutation_timing(data_path,output_path,tree_path)

