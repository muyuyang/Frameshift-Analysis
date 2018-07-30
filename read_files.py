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
            A.find_gap()
            alignments.append(A)
        i += 1
    return alignments


# What information: type, start, end, strand(+/-), attributes
def read_annotation(path):
    with open(path) as f:
        data = f.read().splitlines()
    annotations = []
    i = 0
    # Stop when it starts to read sequences
    while i < len(data) and data[i].replace(' ','') != 'ORIGIN': 
        # print(data[i][0:5],data[i][5])
        if len(data[i]) > 6 and data[i][0:5] == '     ' and data[i][5] != ' ':
            typee = data[i][:21].replace(' ','')
            if typee == 'CDS':
                A = Annotation()
                A.type = typee
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
                # WE DON'T WANT TO INCLUDE PSEUDOGENES
                if 'pseudogene' not in info:
                    annotations.append(A)
            else:
                i += 1
        else:
            i += 1
    return annotations
