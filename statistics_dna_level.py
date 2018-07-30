from classes import Alignment,Annotation

def read_alignment(path):
    with open(path) as f:
        data = f.read()
    lines = data.splitlines()
    i = 0
    alignments = []
    # The alignments in maf format follows:
        # Line 1: a score=x mishap=y
        # Line 2: s geneID start length + a sequence
        # Line 3: s geneID start length + b sequence
    while i < len(lines):
        if 'score' in lines[i]:
            A = Alignment()
            info1 = lines[i+1].split(' ')
            info2 = lines[i+2].split(' ')
            while '' in info1:
                info1.remove('')
            while '' in info2:
                info2.remove('')
            (A.name1,A.start1,A.end1,A.seq1) = (info1[1],int(info1[2]),\
                                        int(info1[2])+int(info1[3])-1,info1[-1])
            (A.name2,A.start2,A.end2,A.seq2) = (info2[1],int(info2[2]),\
                                        int(info2[2])+int(info2[3])-1,info2[-1])
            # print(A.name1,A.start1,A.end1,A.seq1)
            alignments.append(A)
        i += 1
    return alignments

def read_annotation(path):
    with open(path) as f:
        data = f.read().splitlines()
    annotations = []
    i = 0
    # Stop when it starts to read sequences
    while i < len(data) and data[i].replace(' ','') != 'ORIGIN': 
        # print(data[i][0:5],data[i][5])
        if len(data[i]) > 6 and data[i][0:5] == '     ' and data[i][5] != ' ':
            A = Annotation()
            A.type = data[i][:21].replace(' ','')
            A.strand = '-' if 'complement' in data[i] else '+'
            number = data[i][21:].replace('complement','')
            if data[i+1].replace(' ','')[0] != '/':
                number += data[i+1].replace(' ','')
            number = number.replace('(','')
            number = number.replace(')','')
            number = number.replace('join','')
            # print(number)
            number = number.split(',')
            for seg in number:
                pair = seg.split('..')
                A.start.append(int(pair[0]))
                A.end.append(int(pair[1]))
            info = data[i] + '\n'
            i += 1
            while data[i][5] == ' ':
                info += data[i]+'\n'
                i += 1
            A.attribute = info
            # These two types are not very informative.
            if A.type != 'region' and A.type != 'source':
                annotations.append(A)
            if A.type == 'source':
                genome_size = A.end[0]
        else:
            i += 1
    return (annotations,genome_size)


# Represent regions of CDS as a 0/1 list. If L[i] = 1, i is on CDS.
# Makes it easier for counting
def represent_as_list(annotations,genome_size):
    L = [0] * genome_size
    for A in annotations:
        if A.type == 'CDS':
            for j in range(len(A.start)):
                if A.strand == '+':
                    for i in range(A.start[j]+1,A.end[j]):
                        L[i] = 1
                else:
                    for i in range(A.start[j],A.end[j]-1):
                        L[i] = 1
    return L


def count(alignments,annotation1,annotation2,g1_size,g2_size):
    L1 = represent_as_list(annotation1,g1_size)
    L2 = represent_as_list(annotation2,g2_size)
    gap = 0
    substitution = 0
    frameshift = 0
    for A in alignments:
        pos1 = A.start1
        pos2 = A.start2
        # number of consecutive gaps
        gap1 = 0
        gap2 = 0
        for i in range(len(A.seq1)):  # Think about a more elegant way
            # If on CDS
            if L1[pos1] == 1 or L2[pos2] == 1:
                if A.seq1[i] == '-':
                    gap += 1
                    gap1 += 1
                    frameshift += gap2 if gap2 % 3 != 0 else 0
                    gap2 = 0
                elif A.seq2[i] == '-':
                    gap += 1
                    gap2 += 1
                    frameshift += gap1 if gap1 % 3 != 0 else 0
                    gap1 = 0
                else:
                    if A.seq1[i].lower() != A.seq2[i].lower():
                        substitution += 1
                    frameshift += gap1 if gap1 % 3 != 0 else 0
                    frameshift += gap2 if gap2 % 3 != 0 else 0
                    # print(gap1,gap2)
                    (gap1,gap2) = (0,0)
            else:
                (gap1,gap2) = (0,0)
            if A.seq1[i] == '-':
                pos2 += 1
            elif A.seq2[i] == '-':
                pos1 += 1
            else:
                pos1 += 1
                pos2 += 1

    return (substitution,gap,frameshift)


# Input: paths to alignment file and annotation files
# Output: numbers

def statistics(alignment_path,annotation1_path,annotation2_path):
    alignments = read_alignment(alignment_path)
    (annotation1,g1_size) = read_annotation(annotation1_path)
    (annotation2,g2_size) = read_annotation(annotation2_path)
    (substitution,gap,frameshift) = count(alignments,annotation1,annotation2,g1_size,g2_size)
    print('Substitution:',substitution)
    print('Gap:',gap)
    print('Gaps causing frameshift:',frameshift)



alignment_path = '/Users/muyuyang/Desktop/Frith Lab/results/E coli O157_H7/k12-o157-2.maf'
annotation1_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/K-12_annotation.gb'
annotation2_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/O157_H7_annotation.gb'
statistics(alignment_path,annotation1_path,annotation2_path)

