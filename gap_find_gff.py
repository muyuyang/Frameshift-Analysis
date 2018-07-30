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
        self.gap1 = None
        self.gap2 = None

    # Convert the list of gap index to a list of (a,b), where a is the start 
    # index of gap in this sequence and b is the length of gap.
    def compress_gap(self,index):
        if index == []:
            return []
        gap = [(index[0],1)]
        i = 1
        for i in range(1,len(index)):
            (last,length) = gap[-1]
            if index[i] == last + length:
                gap = gap[:-1] + [(last,length+1)]
            else:
                gap += [(index[i],1)]
        return gap

    def find_gap(self):
        gap_index1 = [i for i in range(len(self.seq1)) if self.seq1[i] == '-']
        gap_index2 = [i for i in range(len(self.seq2)) if self.seq2[i] == '-']
        # Write gaps as (a,b), where a is the start index and b is the length
        self.gap1,self.gap2 = (self.compress_gap(gap_index1),self.compress_gap(gap_index2))

class Segment(object): # Find a better name later
    def __init__(self):
        self.type = None
        self.start = None
        self.end = None
        self.strand = None
        self.attribute = None

    def is_in(self,pos):
        if self.strand == '+':
            return pos > self.start and pos <= self.end
        else:
            return pos >= self.start and pos < self.end

    def get_info(self):
        return self.attribute


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
            A.find_gap()
            alignments.append(A)
        i += 1
    return alignments

# Find frameshifts(gaps whose length is not a multiple of 3)
# Write gaps as (a,b,x,y), where a is the start position in seq1, b is the
# start position in seq2, x is the length of gap, and y indicates where the gap is
# (1 if it's on the seq1 and 2 if it's on seq2)
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
                frameshifts.append((start+A.start1-gapNumber1,
                                    start+A.start2-gapNumber2,length,seq))
            if seq == 1:
                gapNumber1 += length
            elif seq == 2:
                gapNumber2 += length
    return frameshifts

# What information: type, start, end, strand(+/-), attributes
def read_annotation(path):
    with open(path) as f:
        data = f.read()
    annotations = []
    for line in data.splitlines()[3:]:
        if line != '':
            info = line.split('\t')
            A = Segment()
            # print(info)
            (A.type,A.start,A.end,A.strand,A.attribute) = \
            (info[2],int(info[3]),int(info[4]),info[6],info[8])
            annotations.append(A)
    return annotations

def find_position(frameshifts,annotation1,annotation2):
    info1 = []
    info2 = []
    for (pos1,pos2,length,seq) in frameshifts:
        I1 = []
        I2 = []
        for A in annotation1:
            if A.is_in(pos1):
                I1.append(A.type + ' ' + A.get_info())
                # I1.append(A.type)
        for A in annotation2:
            if A.is_in(pos2):
                I2.append(A.type + ' ' + A.get_info())
                # I2.append(A.type)
        info1.append(I1)
        info2.append(I2)
    return (info1,info2)


def write_results(path,frameshifts,info1,info2):
    with open(path,'w') as f:
        for i in range(len(frameshifts)):
            # Only print frameshifts on CDS!
            cds1 = [info for info in info1[i] if 'CDS' == info[:3]]
            cds2 = [info for info in info2[i] if 'CDS' == info[:3]]
            if len(cds1) != 0 and len(cds2) != 0:
                f.write(str(frameshifts[i]))
                f.write('\n')
                f.write('On seq 1: ')
                f.write('\n')
                for info in info1[i]:
                    f.write(info)
                    f.write('\n')
                f.write('On seq 2: ')
                f.write('\n')
                for info in info2[i]:
                    f.write(info)
                    f.write('\n')
                f.write('\n')

def gap_find(alignment_path,annotation1_path,annotation2_path,output_path):
    alignments = read_alignment(alignment_path)
    frameshifts = find_frameshifts(alignments)
    # print(frameshifts)
    annotation1 = read_annotation(annotation1_path)
    annotation2 = read_annotation(annotation2_path)
    (info1,info2) = find_position(frameshifts,annotation1,annotation2)
    write_results(output_path,frameshifts,info1,info2)

# alignment_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza sativa chloroplast vs zizania latifolia chloroplast/os-zl-2.maf'
# annotation1_path = '/Users/muyuyang/Desktop/Frith Lab/data/oryza chloroplast/Oryza_sativa_chloroplast_annotation.gff3'
# annotation2_path = '/Users/muyuyang/Desktop/Frith Lab/data/zizania chloroplast/Zizania_latifolia_chloroplast_annotation.gff3'
# output_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza sativa chloroplast vs zizania latifolia chloroplast/os-zl-result-type.txt'

alignment_path = '/Users/muyuyang/Desktop/Frith Lab/results/E coli O157_H7/k12-o157-2.maf'
annotation1_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/K-12_annotation.gff3'
annotation2_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/O157_H7_annotation.gff3'
output_path = '/Users/muyuyang/Desktop/Frith Lab/results/E coli O157_H7/k12-o157-result.txt'


gap_find(alignment_path,annotation1_path,annotation2_path,output_path)

# To-do:
#   1. Compare with outgroup to decide which is ancestral. (Decide where the mutation happens)
#   2. Sort mutations into different categories. 

