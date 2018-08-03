class Alignment(object):
    def __init__(self):
        self.seq1 = None
        self.seq2 = None
        self.start1 = None
        self.end1 = None
        self.start2 = None
        self.end2 = None
        self.name1 = None
        self.name2 = None

    # Find sequence around pos 
    def find_seq(self,pos):
        pos -= self.start1
        i = 0 #Count on aligned sequence
        j = 0 #Count on original sequence
        while j < pos:
            if self.seq1[i] != '-':
                j += 1
            i += 1
        start = i - 50 if i > 50 else 0
        end = i + 50 if self.start1+i+50 < self.end1 else self.end1
        print(self.seq1[start:end])
        print(self.seq2[start:end])

    def is_in(self,pos):
        return self.start1 <= pos and pos < self.end1




# Read file in maf format and store alignments 
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


def find_region_alignment(input_path,pos):
    alignments = read_alignment(input_path)
    for A in alignments:
        if A.is_in(pos):
            A.find_seq(pos)


input_path = '/Users/muyuyang/Desktop/Frith Lab/TEST2/align-2.maf'
pos = 34068

find_region_alignment(input_path,pos)





