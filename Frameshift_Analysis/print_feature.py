# A hardcoded genbank feature extraction program

class Segment(object): # Find a better name later
    def __init__(self):
        # A list of start/end positions -> Deal with join
        self.start = []
        self.end = []
        self.strand = None
        self.attribute = None

    def is_in(self,pos):
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

def read_result(path):
    with open(path) as f:
        data = f.read().splitlines()
    result = []
    for i in range(len(data)//4):
        if 'CDS' in data[i*4+1] or 'CDS' in data[i*4+2]:
            tup = data[i*4][1:-1].split(', ')
            result.append((int(tup[0]),int(tup[1]),int(tup[2]),int(tup[3])))
    return result


def read_genbank_annotation(path):
    with open(path) as f:
        data = f.read().splitlines()
    annotations = []
    i = 0
    while i < len(data):
        if data[i][5:8] == 'CDS':
            A = Segment()
            A.strand = '-' if 'complement' in data[i] else '+'
            number = data[i][21:].replace('complement','')
            if data[i+1].replace(' ','')[0] != '/':
                number += data[i+1].replace(' ','')
            number = number.replace('(','')
            number = number.replace(')','')
            number = number.replace('join','')
            print(number)
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
            annotations.append(A)
        i += 1
    return annotations

def find_annotation(pos,annotations):
    # Speed up -> Use binary search
    for A in annotations:
        if A.is_in(pos):
            return A.get_info()

def write_results(path,frameshifts,info1,info2):
    with open(path,'w') as f:
        for i in range(len(frameshifts)):
            f.write(str(frameshifts[i]))
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


def main(annotation1_path,annotation2_path,result_path,output_path):
    annotations1 = read_genbank_annotation(annotation1_path)
    annotations2 = read_genbank_annotation(annotation2_path)
    result = read_result(result_path)
    (info1,info2) = ([],[])
    for (start1,start2,length,seq) in result:
        info1.append(find_annotation(start1,annotations1))
        info2.append(find_annotation(start2,annotations2))
    write_results(output_path,result,info1,info2)

annotation1_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/K-12_annotation.gb'
annotation2_path = '/Users/muyuyang/Desktop/Frith Lab/data/Escherichia coli/O157_H7_annotation.gb'
result_path = '/Users/muyuyang/Desktop/Frith Lab/results/E coli O157_H7/k12-o157-result-type.txt'
output_path = '/Users/muyuyang/Desktop/Frith Lab/results/E coli O157_H7/k12-o157-genbank-result.txt'

main(annotation1_path,annotation2_path,result_path,output_path)





