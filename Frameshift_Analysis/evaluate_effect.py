#Also, I'm not sure if I'm correctly identifying which amino acids are affected 
#by frameshift. My current thought is to first translate both ancestral and 
#derived DNA sequences into protein sequences and find the corresponding 
#position of the deletion/addition in the protein sequence. Then all the 
#downstream amino acid difference could be caused by frameshift. I use the 
#positions of these downstream amino acids to look for possible features in 
#Uniprot. 

#Input: Alignment, All frameshifts, annotation
from .classes import Protein,Feature
from .find_protein_feature import find_protein_feature,find_protein_id

# Translate DNA sequence into protein sequence

def reverse(d):
    d = d.upper()
    newD = ''
    for i in range(len(d)-1,-1,-1):
        if d[i] == 'A':
            newD += 'T'
        elif d[i] == 'T':
            newD += 'A'
        elif d[i] == 'G':
            newD += 'C'
        elif d[i] == 'C':
            newD += 'G'
    assert(len(d) == len(newD))
    return newD

def dna_translation(d,strand):
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
        'GGG': 'G', 'TAA': '', 'TAG' : '' , 'TGA' : ''}
    prot = ''
    d = d.upper()
    if strand == '-':
        d = reverse(d)
    j = 0
    while j < len(d) and d[j:j+3] != 'ATG':
        j += 1
    for i in range(j,len(d),3):
        if i+3 <= len(d):
            codon = d[i:i+3]
            prot += table[codon]
    return prot


# Input: annotation, original genome sequence
def assemble_sequence(A,originalSeq):
    seq = ''
    for i in range(len(A.start)):
        seq += originalSeq[A.start[i]:A.end[i]]
    return seq

# Find position of frameshift mutation on protein sequence
def find_corresponding_point(pos,A):
    nt_before = 0
    if A.strand == '+':
        for i in range(len(A.start)):
            if pos < A.start[i]:
                break
            elif pos < A.end[i]:
                nt_before += pos - A.start[i]
            elif pos >= A.end[i]:
                nt_before += A.end[i] - A.start[i]
    elif A.strand == '-':
        for i in range(len(A.start)-1,-1,-1):
            if pos > A.end[i]:
                break
            elif pos > A.start[i]:
                nt_before += A.end[i] - pos
            elif pos <= A.start[i]:
                nt_before += A.end[i] - A.start[i]
    aa_pos = nt_before // 3
    return aa_pos


def write_feature_results(path,pos,prot1,prot2,info,features):
    with open(path,'a') as f:
        f.write(str(pos)+'\n')
        f.write(prot1+'\n')
        f.write(prot2+'\n')
        f.write('Comments:\n')
        for c in info.comments:
            f.write(c+'\n')
        keywords = ''
        for k in info.keywords:
            keywords += k + ' '
        f.write('Keywords: '+keywords+'\n')
        f.write('Features:'+'\n')
        for feature in features:
            f.write(feature.type+' ')
            f.write(feature.start+' ')
            f.write(feature.end+' ')
            f.write(feature.description)
            f.write('\n')
        f.write('\n')


def evaluate_effect(frameshifts,oriSeq1,oriSeq2,output_path):
    # Clear the text file first
    with open(output_path,'w') as f:
        f.write('')
    mapping = find_protein_id(frameshifts)
    for F in frameshifts:
        A1 = F.annotation_seq1
        A2 = F.annotation_seq2
        (prot1,prot2) = ('','')
        if A1 != None:
            seq1 = assemble_sequence(A1,oriSeq1)
            prot1 = dna_translation(seq1,A1.strand)
        if A2 != None:
            seq2 = assemble_sequence(A2,oriSeq2)
            prot2 = dna_translation(seq2,A2.strand)
        if A1 != None and A2 != None:
            # print(F.start1)
            # We assume that all amino acids after the insertion/deletion point
            # are affected.
            frameshift_pt = find_corresponding_point(F.start1,A1)
            # print(prot1[:frameshift_pt+1] + '||' + prot1[frameshift_pt+1:])
            # print(prot2)
            try:
                uniprot_id = mapping[F.gene_id]
            except:
                continue
            prot_info = find_protein_feature(uniprot_id)
            if prot_info == None:
                continue
            features = prot_info.features
            affected = []
            for i in range(len(features)):
                # We assume that if the mutated amino acid is within 2 nts of 
                # the functional domain, it is likely to affect its function
                if int(features[i].end) + 2 >= frameshift_pt:
                    affected.append(features[i])
            write_feature_results(output_path,F.start1,prot1,prot2,prot_info,\
                                features)












