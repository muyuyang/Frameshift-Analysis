

# By extracting from 
def main(info1,info2):
    # Starting index of the translated sequence 
    index1 = info1.find('translation=') + 13
    index2 = info2.find('translation=') + 13
    print(info1[index1:-1].replace(' ','').replace('\n',''))
    print(info2[index2:-1].replace(' ','').replace('\n',''))



info1 = '''     gene            212331..213629
                     /gene="tilS"
                     /locus_tag="b0188"
                     /gene_synonym="ECK0187; JW0183; mesJ; yaeN"
                     /db_xref="EcoGene:EG13220"
                     /db_xref="GeneID:944889"

     CDS             212331..213629
                     /gene="tilS"
                     /locus_tag="b0188"
                     /gene_synonym="ECK0187; JW0183; mesJ; yaeN"
                     /EC_number="6.1.1.5"
                     /function="enzyme; tRNA modification"
                     /codon_start=1
                     /transl_table=11
                     /product="tRNA(Ile)-lysidine synthetase"
                     /protein_id="NP_414730.1"
                     /db_xref="ASAP:ABE-0000638"
                     /db_xref="UniProtKB/Swiss-Prot:P52097"
                     /db_xref="EcoGene:EG13220"
                     /db_xref="GeneID:944889"
                     /translation="MTLTLNRQLLTSRQILVAFSGGLDSTVLLHQLVQWRTENPGVAL
                     RAIHVHHGLSANADAWVTHCENVCQQWQVPLVVERVQLAQEGLGIEAQARQARYQAFA
                     RTLLPGEVLVTAQHLDDQCETFLLALKRGSGPAGLSAMAEVSEFAGTRLIRPLLARTR
                     GELVQWARQYDLRWIEDESNQDDSYDRNFLRLRVVPLLQQRWPHFAEATARSAALCAE
                     QESLLDELLADDLAHCQSPQGTLQIVPMLAMSDARRAAIIRRWLAGQNAPMPSRDALV
                     RIWQEVALAREDASPCLRLGAFEIRRYQSQLWWIKSVTGQSENIVPWQTWLQPLELPA
                     GLGSVQLNAGGDIRPPRADEAVSVRFKAPGLLHIVGRNGGRKLKKIWQELGVPPWLRD
                     TTPLLFYGETLIAAAGVFVTQEGVAEGENGVSFVWQKTLS"'''

info2 = '''     gene            215674..216969
                     /gene="tilS"
                     /locus_tag="ECs0190"
                     /db_xref="GeneID:913916"

     CDS             215674..216969
                     /gene="tilS"
                     /locus_tag="ECs0190"
                     /note="Ligates lysine onto the cytidine present at
                     position 34 of the AUA codon-specific tRNA(Ile) that
                     contains the anticodon CAU; ATP-dependent; responsible for
                     modifying the wobble-base of the CAU anticodon of tRNAIle
                     such that it exhibits proper recognition of the AUA codon
                     rather than the AUG codon and is in turn properly
                     recognized by isoleucyl-tRNA synthetase"
                     /codon_start=1
                     /transl_table=11
                     /product="tRNA(Ile)-lysidine synthetase"
                     /protein_id="NP_308217.1"
                     /db_xref="GeneID:913916"
                     /translation="MTLTLNRQLLTSRQILVAFSGGLDSTVLLHQLVQWRTENPGGAL
                     RAIHVHHGLSANADAWVTHCENVCQQWQVPLVVERVQLAQEGLGIEAQARQARYQAFA
                     RTLLPGEVLVTAQHLDDQCETFLLALKRGSGPAGLSAMAEVSEFAGTRLIRPLLARTR
                     GELAQWALAHGLRWIEDESNQDDSYDRNFLRLRVVPLLQQRWPHFAETTARSAALCAE
                     QESLLDELLADDLAHCQSPQGTLQIVPMLAMSDARRAAIIRRWLAGQNAPMPSRDALV
                     RIWQEVALAREDASPCLRLGAFEIRRYQSQLWWIKSVTGQSENIVPWQTWLQPLELPA
                     GQGSVQLNAGGDIRPPRADEAVSVRFKAPGLLHIVGRNGGRKLKKIWQELGVPPWLRD
                     TTPLLFYGETLIAAAGGFVTQEGVAEGENGISFVWQRNA"'''

# main(info1,info2)
