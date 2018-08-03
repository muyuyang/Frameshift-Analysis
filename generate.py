import os 
from gap_find import gap_find
from determine_ancestry import ancestry
import random
from classes import Alignment,Frameshift

def last_db(input1_path,db_name,output_path):
    command = 'lastdb -P0 -uMAM4 -R01 '+output_path+db_name+' '+input1_path
    print(command)
    os.system(command)



def last_train(input2_path,output_path,training_name,db_name,e_threshold):
    command = 'last-train -P0 --revsym -E'+str(e_threshold)+' -C2 '+output_path+\
                db_name+' '+input2_path+' > '+output_path+training_name
    print(command)
    os.system(command)

def many_to_one_alignment(input2_path,output_path,db_name,training_name,e_threshold,align_name):
    command = 'lastal -E'+str(e_threshold)+' -C2 -p '+output_path+training_name+\
            ' '+output_path+db_name+' '+input2_path+' | last-split > '+\
            output_path+align_name
    print(command)
    os.system(command)

def one_to_one_alignment(output_path,align1_name,align2_name):
    command = 'maf-swap '+output_path+align1_name+' |'+'\n'+\
             'last-split -m1 |' +'\n'\
             'maf-swap > '+output_path+align2_name
    print(command)
    os.system(command)

def plot(output_path,align_name):
    command = 'last-dotplot '+output_path+align_name+' '+\
            output_path+align_name[:-3]+'png'
    os.system(command)



def main(input1_path,input2_path,output_path,db_name,training_name,align1_name,align2_name):
    last_db(input1_path,db_name,output_path)
    e_threshold = 4 * 10**8
    last_train(input2_path,output_path,training_name,db_name,e_threshold)
    many_to_one_alignment(input2_path,output_path,db_name,training_name,e_threshold,align1_name)
    one_to_one_alignment(output_path,align1_name,align2_name)
    plot(output_path,align2_name)


# species = ['sativa','nivara','rufipogon','minuta','latifolia','alta','australiensis','brachyantha','coarctata']
# species = ['sativa','nivara','minuta','latifolia','alta','australiensis','brachyantha','coarctata']
species = ['sativa','brachyantha']

# for a in range(len(species)):
#     for b in range(a+1,len(species)):
#         i = species[a]
#         j = species[b]
#         input1_path = '/Users/muyuyang/Desktop/Frith\ Lab/data/oryza\ chloroplast/Oryza_%s_chloroplast.fasta' % i
#         input2_path = '/Users/muyuyang/Desktop/Frith\ Lab/data/oryza\ chloroplast/Oryza_%s_chloroplast.fasta' % j
#         output_path = '/Users/muyuyang/Desktop/Frith\ Lab/results/oryza\ chloroplast\ %s\ vs\ %s/' % (i[0].upper()+i[1:],j[0].upper()+j[1:])
#         if not os.path.exists(output_path.replace('\\','')):
#             os.system('mkdir %s' % output_path) 
#             db_name = 'o' + i[0]
#             training_name = 'o%s-o%s.mat' % (i[0],j[0])
#             align1_name = 'o%s-o%s-1.maf' % (i[0],j[0])
#             align2_name = 'o%s-o%s-2.maf' % (i[0],j[0])
#             main(input1_path,input2_path,output_path,db_name,training_name,align1_name,align2_name)



for a in range(len(species)):
    for b in range(a+1,len(species)):
        i = species[a]
        j = species[b]
        align_name = 'o%s-o%s-2.maf' % (i[0],j[0])
        alignment_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast %s vs %s/' % (i[0].upper()+i[1:],j[0].upper()+j[1:])+ align_name
        annotation1_path = '/Users/muyuyang/Desktop/Frith Lab/data/oryza chloroplast/Oryza_%s_chloroplast_annotation.gb' % i
        annotation2_path = '/Users/muyuyang/Desktop/Frith Lab/data/oryza chloroplast/Oryza_%s_chloroplast_annotation.gb' % j
        output_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast %s vs %s/o%s-o%s-result.txt' % (i[0].upper()+i[1:],j[0].upper()+j[1:],i[0],j[0])
        gap_find(alignment_path,annotation1_path,annotation2_path,output_path)


common = []

# Get mutation using all intersections
# for a in range(1,len(species)):
#     for b in range(a+1,len(species)):
#         i = species[a]
#         j = species[b]
#         input1_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast Sativa vs %s/os-o%s-result-type.txt' % (i[0].upper()+i[1:],i[0])
#         input2_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast Sativa vs %s/os-o%s-result-type.txt' % (j[0].upper()+j[1:],j[0])
#         # print(i,j)
#         common.append(ancestry(input1_path,input2_path))
# common = [set(i) for i in common if i != []]
# mutations = set.intersection(*common)
# print(list(mutations))


# Get mutation using only the closest species and an outgroup
# for i in range(len(species)-2):
#     a = species[i]
#     b = species[i+1] # closest species
#     c = species[i+2] # outgroup 
#     # c = species[random.randint(i+2,len(species)-1)]
#     input1_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast %s vs %s/o%s-o%s-result-type.txt' % (a[0].upper()+a[1:],b[0].upper()+b[1:],a[0],b[0])
#     input2_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast %s vs %s/o%s-o%s-result-type.txt' % (a[0].upper()+a[1:],c[0].upper()+c[1:],a[0],c[0])
#     mutation = ancestry(input1_path,input2_path)
#     print(a,mutation)


# for i in range(2,len(species)):
#     a = species[i]
#     b = species[i-1] # closest species
#     c = species[i-2] # outgroup 
#     # c = species[random.randint(i+2,len(species)-1)]
#     input1_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast %s vs %s/o%s-o%s-result-type.txt' % (b[0].upper()+b[1:],a[0].upper()+a[1:],b[0],a[0])
#     input2_path = '/Users/muyuyang/Desktop/Frith Lab/results/oryza chloroplast %s vs %s/o%s-o%s-result-type.txt' % (c[0].upper()+c[1:],a[0].upper()+a[1:],c[0],a[0])
#     mutation = ancestry(input1_path,input2_path)
#     print(a,mutation)






# Store in a form?




