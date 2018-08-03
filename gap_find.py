# Uses genbank format

from classes import Alignment,Frameshift,Annotation
from read_files import *


# Find frameshifts(gaps whose length is not a multiple of 3)
# Write gaps as (a,b,x,y), where a is the start position in seq1, b is the
# start position in seq2, x is the length of gap, and y indicates where the gap 
# is (1 if it's on the seq1 and 2 if it's on seq2)
def find_frameshifts(alignments):
    frameshifts = []
    for A in alignments:
        gapNumber1 = 0
        gapNumber2 = 0
        gap1 = [(i,j,1) for (i,j) in A.gap1]
        gap2 = [(i,j,2) for (i,j) in A.gap2]
        gaps = sorted(gap1+gap2)
        for (start,length,seq) in gaps:
            if length % 3 != 0:
                F = Frameshift(start+A.start1-gapNumber1,
                                    start+A.start2-gapNumber2,length,seq)
                F.alignment = A
                frameshifts.append(F)
            if seq == 1:
                gapNumber1 += length
            elif seq == 2:
                gapNumber2 += length
    return frameshifts

def find_individual_annotation(pos,length,annotation):
    start = 0
    end = len(annotation)-1
    info = ''
    length = len(annotation)
    while length > 0 and start <= end:
        mid = (start+end) // 2
        A = annotation[mid]
        if A.is_in(pos,length):
            return A
        elif pos < A.start[0]:
            end = mid - 1
        else:
            start = mid + 1
    return None

# def find_individual_annotation_slow(pos,length,annotation):
#     # Speed up -> Use binary search
#     for A in annotation:
#         if A.is_in(pos,length):
#             return A

def find_annotation(F,annotation1,annotation2):
    A1 = find_individual_annotation(F.start1,F.length,annotation1)
    A2 = find_individual_annotation(F.start2,F.length,annotation2)
    F.annotation_seq1 = A1
    F.annotation_seq2 = A2
    info1 = None if A1 == None else A1.get_info()
    info2 = None if A2 == None else A2.get_info()
    return (info1,info2)



def write_results(path,frameshifts,info1,info2):
    with open(path,'w') as f:
        for i in range(len(frameshifts)):
            if info1[i] != None or info2[i] != None:
                F = frameshifts[i]
                f.write(str((F.start1,F.start2,F.length,F.seq)))
                f.write('\n')
                f.write('On seq 1: ')
                f.write('\n')
                f.write(str(info1[i]))
                if info1[i] == None:
                    f.write('\n')
                f.write('On seq 2: ')
                f.write('\n')
                f.write(str(info2[i]))
                if info2[i] == None:
                    f.write('\n')
                f.write('\n\n')


def gap_find(alignments,annotation1,annotation2,output_path):
    frameshifts = find_frameshifts(alignments)
    # print(frameshifts)
    (info1,info2) = ([],[])
    for F in frameshifts:
        (i1,i2) = find_annotation(F,annotation1,annotation2)
        info1.append(i1)
        info2.append(i2)
    write_results(output_path+'results.txt',frameshifts,info1,info2)
    return frameshifts

# alignment_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast Sativa vs Coarctata/os-oc-2.maf'
# annotation1_path = '/Users/muyuyang/Desktop/Frith Lab/data/oryza chloroplast/Oryza_sativa_chloroplast.gb'
# annotation2_path = '/Users/muyuyang/Desktop/Frith Lab/data/oryza chloroplast/Oryza_sativa_chloroplast.gb'
# output_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast Sativa vs Coarctata/os-oc-result-type.txt'

alignment_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/results/align-2.maf'
annotation1_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/data/Oryza minuta.gb'
annotation2_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/data/Oryza sativa.gb'
output_path = '/Users/muyuyang/Desktop/Frith Lab/TEST_TIMING/results/testingtesting.txt'
alignments = read_alignment(alignment_path)
annotation1 = read_annotation(annotation1_path)
annotation2 = read_annotation(annotation2_path)
gap_find(alignments,annotation1,annotation2,output_path)

# To-do:
#   1. Compare with outgroup to decide which is ancestral. (Decide where the mutation happens)
#   2. Sort mutations into different categories. 

