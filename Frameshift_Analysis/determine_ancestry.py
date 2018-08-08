#Determine which species is ancestral for each mutation by comparing with an outgroup.

#Input is the result txt files. [A vs B] and [A vs C]
#Common mutation -> Mutation happens in A.
#Other mutations -> Cannot conclude that mutation happens in B
#                -> Need to double check alignment

def read_Frameshifts(path):
    F = []
    with open(path) as f:
        data = f.read().splitlines()
    for i in range(len(data)//4):
        # We are only interested in frameshifts on genes
        if 'CDS' in data[4*i+1] or 'gene' in data[4*i+2]: 
            tupl = data[i*4][1:-1].split(', ')
            F.append((int(tupl[0]),int(tupl[1]),int(tupl[2]),int(tupl[3])))
    return F


def find_common_mutation(F1,F2):
    common = []
    for (a1,b1,length1,seq1) in F1:
        for (a2,b2,length2,seq2) in F2:
            if a1 == a2 and length1 == length2 and seq1 == seq2:
                typ = 'Deletion' if seq1 == 1 else 'Addition'
                common.append((a1,length1,typ))
    return common


def ancestry(input1_path,input2_path):
    frameshifts1 = read_Frameshifts(input1_path)
    frameshifts2 = read_Frameshifts(input2_path)
    common = find_common_mutation(frameshifts1,frameshifts2) # -> mutations in common are derived in the branch to A
    # Is it possible to construct a tree here? -> But by choosing 'outgroup' we 
    # are making some assumptions on the evolutionary relationship

    # Wait - Is it better to use MSA???

    # print(frameshifts1)
    # print(frameshifts2)

    # print(common)
    return common


input1_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast Sativa vs Coarctata/os-oc-result.txt'
input2_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast Sativa vs Australiensis/os-oa-result.txt'

# ancestry(input1_path,input2_path)

