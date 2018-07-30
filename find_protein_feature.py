import urllib.request,urllib.parse
from classes import Protein,Feature


# Get the accession number from Genbank formatted information
# Input: information in Genbank format
# Output: Accession number
def extract_accession(info):
    index = info.find('GeneID:')
    if index != -1:
        index_start = index + len('GeneID:')
        index_end = info.find('"',index_start)
        return info[index_start:index_end]
    else:
        return None

#Revised from https://www.uniprot.org/help/api_idmapping
def id_mapping(gene_id):
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from':'P_ENTREZGENEID',
    'to':' ',
    'format':'tab',
    'query':gene_id
    }
    data = urllib.parse.urlencode(params).encode("utf-8")
    request = urllib.request.Request(url, data)
    response = urllib.request.urlopen(request)
    page = response.read().decode('utf-8')
    mapping = dict()
    for line in page.splitlines():
        from_to = line.split('\t')
        gene_id = from_to[0]
        uniprot_id = from_to[1]
        mapping[gene_id] = uniprot_id
    return mapping



# Input: A list L of strings
# Output: L with all empty strings removed
def delete_empty_string(L):
    newL = [string for string in L if string != '']
    return newL

def is_useful_comments(s):
    useful = ['CAUTION','DISEASE','FUNCTION','PATHWAY']
    for keyword in useful:
        if keyword in s:
            return True
    return False

def is_useful_feature(s):
    useful = ['CHAIN','REPEAT','INIT_MET','PROPEP','PEPTIDE','COILED','COMPBIAS',
            'VARIANT','HELIX','TURN','STRAND']
    for keyword in useful:
        if keyword in s:
            return False
    return True

def find_protein_id(frameshifts):
    query = ''
    for F in frameshifts:
        A = F.annotation_seq1
        if A != None:
            info = A.attribute
            gene_id = extract_accession(info)
            query += gene_id + ' '
            F.gene_id = gene_id
    query = query.strip()
    mapping = id_mapping(query)
    return mapping



# Input: (str1,str2) -> str1 is the genbank annotation in reference genome
#                       str2 is the genbank annotation in derived genome
# Output: Protein object -> protein information of the given CDS in reference genome
        # None if none is found
def find_protein_feature(uniprot_id):
    # print(uniprot_id)
    if uniprot_id == None:
        return None
    url = 'https://www.uniprot.org/uniprot/%s.txt' % uniprot_id
    try:
        response = urllib.request.urlopen(url)
    except urllib.error.HTTPError:
        return None
    data = response.read().decode('utf-8')
    # print(data)
    prot = Protein()
    lines =  data.splitlines()
    i = 0
    while i < len(lines):
        lineType = lines[i][:2]
        if lineType == 'DE':
            index = lines[i].find('RecName: Full=')
            if index != -1:
                prot.name = lines[i][index+len('RecName: Full='):-1]
            i += 1
        elif lineType == 'CC':
            if lines[i][5:8] == '-!-' and is_useful_comments(lines[i]):
                message = lines[i][9:]
                i += 1
                while i < len(lines) and lines[i][5:8] == '   ':
                    message += lines[i][9:] + ' '
                    i += 1
                prot.comments.append(message)
            else:
                i += 1
        elif lineType == 'KW':
            while lines[i][:2] == 'KW': 
                for kw in lines[i][5:].split('; '):
                    prot.keywords.append(kw.strip(';').strip('.'))
                i += 1
        elif lineType == 'FT':
            if lines[i][5] != ' ' and is_useful_feature(lines[i]):
                feature = Feature()
                ft = delete_empty_string(lines[i][5:34].split(' '))
                feature.type = ft[0]
                feature.start = ft[1]
                feature.end = ft[2]
                feature.description = lines[i][34:]
                i += 1
                while i < len(lines) and lines[i][5] == ' ':
                    feature.description += lines[i][34:]+' '
                    i += 1
                prot.features.append(feature)
            else:
                i += 1
        else:
            i += 1
    # print('Name:',prot.name)
    # print('Comments:',prot.comments)
    # print('Keywords:',prot.keywords)
    # print('Features:')
    # for ft in prot.features:
    #     print(ft.type,ft.start,ft.end,ft.description)
    return prot


