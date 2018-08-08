from .classes import Alignment,Annotation
from .read_files import *


def is_silent_mutation(codon1,codon2):
    codon1 = codon1.upper()
    codon2 = codon2.upper()
    # Don't count as substitution
    if len(codon1) < 3 or len(codon2) < 3:
        return True
    table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
            'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
            'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
            'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
            'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
            'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
            'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
            'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
            'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
            'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
            'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
            'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
            'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
            'GGG': 'G', 'TAA': 'STOP', 'TAG' : 'STOP' , 'TGA' : 'STOP'}
    return table[codon1] == table[codon2]


# Without knowing protein sequence

# For each CDS in annotation of reference genome:
# Find region of reference genome in alignment -> Find region in derived genome
# -> Count differences caused by substitution and frameshifts

# Find a table that matches position in original seq to position in aligned seq
def find_position_match(alignments):
    # D[pos] = (A,i) -> A is the alignment, i is the pos on aligned seq
    D1 = dict()
    D2 = dict()
    for A in alignments:
        pos1 = A.start1
        pos2 = A.start2
        for i in range(len(A.seq1)):
            if A.seq1[i] != '-':
                D1[pos1] = (A,i)
                pos1 += 1
            if A.seq2[i] != '-':
                D2[pos2] = (A,i)
                pos2 += 1
    return (D1,D2)

def flip(s):
    newS = ''
    for c in s:
        if c == 'A':
            newS += 'T'
        elif c == 'T':
            newS += 'A'
        elif c == 'C':
            newS += 'G'
        elif c == 'G':
            newS += 'C'
        elif c == '-':
            newS += '-'
    return newS

# Get the aligned sequences of the CDS A
def get_aligned_sequence(A,D1,D2):
    seq1 = ''
    seq2 = ''
    for i in range(len(A.start)):
        # Start and end in the original sequence
        (start,end) = (A.start[i]-1,A.end[i])
        # Find start and end in the alignment sequence
        try:
            (start_A,start_pos) = D1[start]
            (end_A,end_pos) = D1[end]
        except: # Not in aligned region
            return (None,None)
        if start_A == end_A:
            seq1 += start_A.seq1[start_pos:end_pos]
            seq2 += start_A.seq2[start_pos:end_pos]
        else:
            seq1 += start_A.seq1[start_pos:] + end_A.seq1[:end_pos]
            seq2 += start_A.seq2[start_pos:] + end_A.seq2[:end_pos]
    # print(seq1)
    # print(seq2)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    if A.strand == '-':
        seq1 = flip(seq1[::-1])
        seq2 = flip(seq2[::-1])
    return (seq1,seq2)

# Delete gaps in sequence s
def delete_gaps(s):
    newS = ''
    for c in s:
        if c != '-':
            newS += c
    return newS

def write_statistics(path,total_sub,total_fs,in_region,not_in_aligned_region):
    with open(path,'w') as f:
        f.write('Total number of substitutions: '+str(total_sub))
        f.write('\n')
        f.write('Total number of amino acids affected by frameshift: '+\
                str(total_fs))
        f.write('\n')
        f.write('Total number of CDS within the aligned regions: '+\
            str(in_region))
        f.write('\n')
        f.write('Total number of CDS outside the aligned regions: '+\
                str(not_in_aligned_region))


def statistics(alignments,annotation1,annotation2,path):
    (D1,D2) = find_position_match(alignments)
    substitutions = []
    deletions = []
    insertions = []
    frameshifts = []
    not_in_aligned_region = 0
    for A in annotation1:
        (seq1,seq2) = get_aligned_sequence(A,D1,D2)
        if seq1 == None:
            not_in_aligned_region += 1
            continue
        i = 0
        j = 0
        # I'm calculating how much do mutations affect the original seq
        substitution = 0
        deletion = 0 # Frameshift not included
        insertion = 0
        frameshift = 0 # Number of aa affected by frameshift
        while i < len(seq1) and j < len(seq2):
            if '-' not in seq1[i:i+3] and '-' not in seq2[j:j+3]:
                codon1 = seq1[i:i+3] 
                codon2 = seq2[j:j+3]
                if not is_silent_mutation(codon1,codon2):
                    substitution += 1
                i += 3
                j += 3
            elif seq1[i] == '-':
                gap = 0
                while i < len(seq1) and seq1[i] == '-':
                    gap += 1
                    i += 1
                    j += 1
                if gap % 3 == 0:
                    insertion += gap // 3
                else:
                    #FRAMESHIFT!
                    frameshift = len(delete_gaps(seq1[i:])) // 3
                    break
            elif seq2[j] == '-':
                gap = 0
                while j < len(seq2) and seq2[j] == '-':
                    gap += 1
                    i += 1
                    j += 1
                if gap % 3 == 0:
                    deletion += gap // 3
                else:
                    #FRAMESHIFT!!
                    frameshift = len(delete_gaps(seq1[i:])) // 3
                    break
            else:
                i += 1
                j += 1
        substitutions.append(substitution)
        deletions.append(deletion)
        insertions.append(insertion)
        frameshifts.append(frameshift)
    # print(substitutions)
    # print(deletions)
    # print(insertions)
    # print(frameshifts)

    total_sub = sum(substitutions)
    total_fs = sum(frameshifts)
    in_region = len(substitutions)
    # print(total_sub,total_fs,in_region,not_in_aligned_region)
    write_statistics(path,total_sub,total_fs,in_region,not_in_aligned_region)


# alignments = read_alignment('/Users/muyuyang/Desktop/Frith Lab/TEST/align-2.maf')
# annotation1 = read_annotation('/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/K-12_annotation.gb')
# annotation2 = read_annotation('/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/K-12_annotation.gb')



# statistics(alignments,annotation1,annotation2)























