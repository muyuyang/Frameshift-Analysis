import os 

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



def align(input1_path,input2_path,output_path):
    db_name = 'db'
    training_name = 'align.mat'
    last_db(input1_path,db_name,output_path)
    e_threshold = 4 * 10 ** 8
    last_train(input2_path,output_path,training_name,db_name,e_threshold)
    align1_name = 'align-1.maf'
    many_to_one_alignment(input2_path,output_path,db_name,training_name,e_threshold,align1_name)
    align2_name = 'align-2.maf'
    one_to_one_alignment(output_path,align1_name,align2_name)
    plot(output_path,align2_name)

