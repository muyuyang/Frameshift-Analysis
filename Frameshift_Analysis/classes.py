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
        self.gap1,self.gap2 = (self.compress_gap(gap_index1),\
                            self.compress_gap(gap_index2))

class Annotation(object): # Find a better name later
    def __init__(self):
        self.type = None
        self.start = []
        self.end = []
        self.strand = None
        self.attribute = None

    def is_in(self,pos,length):
        for i in range(len(self.start)):
            if self.strand == '+':
                if pos > self.start[i] and pos <= self.end[i]:
                    return True
            else:
                if pos >= self.start[i] and pos < self.end[i]:
                    return True
        return False

    def get_info(self):
        return self.attribute

class Frameshift(object):
    def __init__(self,start1,start2,length,seq):
        self.start1 = start1
        self.start2 = start2
        self.length = length #length of gap
        self.seq = seq
        self.alignment = None  #the alignment region it is in
        self.annotation_seq1 = None #the annotations that include this region 
        self.annotation_seq2 = None
        self.gene_id = None
        self.uniprot_id = None



    #print the actual sequence alignment around the frameshift mutation
    def find_alignment_region(self):
        A = self.alignment_seq1
        pos = self.start1
        pos -= A.start1
        i = 0 #Count on aligned sequence
        j = 0 #Count on original sequence
        while j < pos:
            if A.seq1[i] != '-':
                j += 1
            i += 1
        start = i - 50 if i > 50 else 0
        end = i + 50 if A.start1+i+50 < A.end1 else A.end1
        print(A.seq1[start:end])
        print(A.seq2[start:end])


class Protein(object):
    def __init__(self):
        self.name = None #str
        self.comments = []
        self.keywords = []
        self.features = []

class Feature(object):
    def __init__(self):
        self.type = None
        self.start = None
        self.end = None
        self.description = None


