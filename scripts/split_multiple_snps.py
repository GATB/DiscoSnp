import sys

""" 
Given a discoSnp fasta file result, splits the bubbles into bubbles containing a unique SNP. Upper and lower case of each SNP is STRICTLY the same excepted at the variant position. Thus the phasing is lost here.
This is a simplification needed for validation purposes

Examples: 
    >SNP_higher_path_701|P_1:30_C/T|high|nb_pol_1|left_unitig_length_226|right_unitig_length_164|left_contig_length_475|right_contig_length_1438|C1_31|C2_0|Q1_0|Q2_0|G1_0/0:6,98,624|G2_1/1:644,101,6|rank_1
taaaaatgaacatattaAGTACTCACAATTCATCTCTTCCCTGTCACCTGACTGCATCTGACTGCATACCACGTACACcag
>SNP_lower_path_701|P_1:30_C/T|high|nb_pol_1|left_unitig_length_226|right_unitig_length_164|left_contig_length_475|right_contig_length_1438|C1_0|C2_32|Q1_0|Q2_0|G1_0/0:6,98,624|G2_1/1:644,101,6|rank_1
taaaaatgaacatattaAGTACTCACAATTCATCTCTTCCCTGTCACTTGACTGCATCTGACTGCATACCACGTACACcagg
becomes: 
>SNP_higher_path_701_1|P_1:30_C/T
AGTACTCACAATTCATCTCTTCCCTGTCACCTGACTGCATCTGACTGCATACCACGTACA
>SNP_lower_path_701_1|P_1:30_C/T
AGTACTCACAATTCATCTCTTCCCTGTCACTTGACTGCATCTGACTGCATACCACGTACA

AND:
>SNP_higher_path_711|P_1:30_T/G,P_2:32_G/C|high|nb_pol_2|left_unitig_length_18|right_unitig_length_61|left_contig_length_649|right_contig_length_203|C1_0|C2_40|Q1_0|Q2_0|G1_1/1:464,74,5|G2_0/0:6,125,804|rank_1
tcattcactcattcgctcactcattcactcattcattctttactcattcactcattcactcactcactgattcactcattcattcattcactcacacattcacttattcactcacttactcattcactcactcactcacgtccagggcctgcttggtgtccactgctgctgcaggcaactggcttcagttgggggaatgtatttctagtttacctttccctaaaggtgtggggtagggacgatggtcctctcctgagtctggcacagtgctgcagtctctgcccttgtgtacccgaaaacctgaaggtcaaggtcactgctgcaaccgatgccctgctgggaaagagaaagccttgggctgcctcagctttggcctccacagggtctctgcggtgcccttcagggcaccagagaaggtcctttcctgctacccaggctcccagtgcctgggctgagcctttgccccatccacccctcactgcctgtctcacgccggggcctttgcacccgccgatccctctgccgggagcaccctgctgcgcctccttgctccagagaagcctcccctgcccccacccccacccgtgagctgctccctgtgattttctgtggcaccacctgtgtcactatgcaggtcacagtaacagcaATGCGCTGGCCCCGTCCAGACTTCCTGTCTTTGGCTCCATAGGACGGTCCCCAGGAAGAAGAAccctgtctacccggcccagttcccccttctccagcaggcctgggacagagccagagcagataggcaaggctccgtccccaaactcattctgaggaaaccaccacagccatcctcctggggtcagcctggggcttttcccgggcatggcttatatccacagaaatgagaactgcgccaggcgcggtggctcacacctgtaatcc
>SNP_lower_path_711|P_1:30_T/G,P_2:32_G/C|high|nb_pol_2|left_unitig_length_18|right_unitig_length_61|left_contig_length_649|right_contig_length_203|C1_23|C2_0|Q1_0|Q2_0|G1_1/1:464,74,5|G2_0/0:6,125,804|rank_1
tcattcactcattcgctcactcattcactcattcattctttactcattcactcattcactcactcactgattcactcattcattcattcactcacacattcacttattcactcacttactcattcactcactcactcacgtccagggcctgcttggtgtccactgctgctgcaggcaactggcttcagttgggggaatgtatttctagtttacctttccctaaaggtgtggggtagggacgatggtcctctcctgagtctggcacagtgctgcagtctctgcccttgtgtacccgaaaacctgaaggtcaaggtcactgctgcaaccgatgccctgctgggaaagagaaagccttgggctgcctcagctttggcctccacagggtctctgcggtgcccttcagggcaccagagaaggtcctttcctgctacccaggctcccagtgcctgggctgagcctttgccccatccacccctcactgcctgtctcacgccggggcctttgcacccgccgatccctctgccgggagcaccctgctgcgcctccttgctccagagaagcctcccctgcccccacccccacccgtgagctgctccctgtgattttctgtggcaccacctgtgtcactatgcaggtcacagtaacagcaATGCGCTGGCCCCGTCCAGACTTCCTGTCTGTCGCTCCATAGGACGGTCCCCAGGAAGAAGAAccctgtctacccggcccagttcccccttctccagcaggcctgggacagagccagagcagataggcaaggctccgtccccaaactcattctgaggaaaccaccacagccatcctcctggggtcagcctggggcttttcccgggcatggcttatatccacagaaatgagaactgcgccaggcgcggtggctcacacctgtaatcc
becomes 
>SNP_higher_path_711_1|P_1:30_T/G
ATGCGCTGGCCCCGTCCAGACTTCCTGTCTTTGGCTCCATAGGACGGTCCCCAGGAAGAAG
>SNP_lower_path_711_1|P_1:30_T/G
ATGCGCTGGCCCCGTCCAGACTTCCTGTCTGTGGCTCCATAGGACGGTCCCCAGGAAGAAG
>SNP_higher_path_711_2|P_2:32_G/C
GCGCTGGCCCCGTCCAGACTTCCTGTCTTTGGCTCCATAGGACGGTCCCCAGGAAGAAGAA
>SNP_lower_path_711_2|P_2:32_G/C
GCGCTGGCCCCGTCCAGACTTCCTGTCTTTCGCTCCATAGGACGGTCCCCAGGAAGAAGAA
"""
file = open(sys.argv[1])

def get_maj_seq(seq):
    """
    returns the upper case letters from a sequence @seq
    """
    res=""
    for l in seq:
        if l>='A' and l<='Z': 
            res+=l
    return res

delta=30
count=0
while True: 
    com1=file.readline()
    if not com1: break
    count+=1
    com1=com1.rstrip()               #>SNP_higher_path_30|P_1:30_A/T,P_2:55_G/T|low|nb_pol_2|left_unitig_length_6|right_unitig_length_0
    seq1=get_maj_seq(file.readline().rstrip())
    com2=file.readline().rstrip()
    seq2=get_maj_seq(file.readline().rstrip())
    str_nb_pol=com1.split('|')[3]
    assert str_nb_pol.startswith("nb_pol")
    nb_pol=int(com1.split('|')[3].split('_')[-1])
    if com1.startswith(">INDEL"):
        print (com1)
        print (seq1)
        print (com2)
        print (seq2)
        continue
    
    all_pos_pol=com1.split('|')[1].split(',')
    splited_com1=com1.split('|')
    splited_com2=com2.split('|')
    for pol in range(nb_pol): 
        pos_pol=int(all_pos_pol[pol].split(':')[1].split('_')[0])
        print (com1)
        alt_allele=all_pos_pol[pol].split(':')[1].split('_')[1].split('/')[1]
        subseq1=seq1[pos_pol-delta:pos_pol+delta+1]                                     # get the upper case sequence, including the current variant and containing the sequence surrounding this variant
        subseq2=seq1[pos_pol-delta:pos_pol]+alt_allele+seq1[pos_pol+1:pos_pol+delta+1]  # get the upper case sequence, including the current variant and containing the sequence from sequence1 containing this variant. 
        subcom1=splited_com1[0]+'_'+str(pol+1)+'|'+all_pos_pol[pol]
        subcom2=splited_com2[0]+'_'+str(pol+1)+'|'+all_pos_pol[pol]
        print (subcom1)
        print (subseq1)
        print (subcom2)
        print (subseq2)
        
        
        
    
    
    
    