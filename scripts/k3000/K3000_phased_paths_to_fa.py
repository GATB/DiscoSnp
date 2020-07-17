#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Creation of a FA file from set of phased facts, as obtained by the ILP model:
#ID #path ("variant id"_"distance between the end of the previous bubble and the start of the current one) #abundance
1       1054l_0;-826l_-26;-914l_-31;-434h_-36;-369h_-15;675h_-5;-543l_-17;-386l_-27;326l_-1;88h_-10;359l_-21;542l_-28;-540h_-36;-153h_-2;382l_-3;1077h_-37;40h_56;413l_-40;-744h_-32;-510h_2;-636l_-22;-1076l_-22;646h_39;897l_-32;968h_-34;1130h_18;578l_-28;812l_-8;-547h_-10;412h_-14;435h_-3;113h_-16;-975l_-35;-807l_8;-1094l_-8;-447h_-40;-802l_-23;-779l_-19;1066h_-25;-429l_11;701h_-35;180l_-30;-563h_6;67h_5;169h_-36;-864l_-21;-219l_-14;282h_5;567h_-31;453l_-24;740l_-35;-255l_-29;723h_-12;1050h_-24;-204h_17;303h_-12;-237l_-34;236h_-30;-408h_-32;-102l_-6;-1128l_-29;-474h_-12;1092h_-14;-491h_-35;889h_-35;45l_-8;1182l_-40;477l_-28;-172l_-13;1060h_-13;221h_-26;157h_-31;1091l_-17;-1021l_-2;846h_-33;1014h_-36;964h_-38;   162

Outputs: 
>path_ID|abundance|info about starting and ending positions of each allele on the global sequence (SP = Sequence positions), and the Bubble Positions (BP). For each allele in the fact we store the distance between the bubble start (upper case letter and the end of the previous bubble (also upper case letter). We add the length of the bubble (upper case letter).
ACGTN* sequence
'''


__author__ = "Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr"


import sys
import K3000_common as kc




    
def generate_sequence_paths(sequences, compacted_fact_file_name, hamming_max=1):
    '''
    Given a set of indexed sequences, the k value and some compacted facts, prints the corresponding fasta seequences
    Enables also to validate the compacted facts. 
    '''
    non_writen_list=""
    mfile = open(compacted_fact_file_name)
    nb_non_writen=0
    for line in mfile: 
        # 1       1054l_0;-826l_-26;-914l_-31;-434h_-36;-369h_-15;675h_-5;-543l_-17;-386l_-27;326l_-1;88h_-10;359l_-21;542l_-28;-540h_-36;-153h_-2;382l_-3;1077h_-37;40h_56;413l_-40;-744h_-32;-510h_2;-636l_-22;-1076l_-22;646h_39;897l_-32;968h_-34;1130h_18;578l_-28;812l_-8;-547h_-10;412h_-14;435h_-3;113h_-16;-975l_-35;-807l_8;-1094l_-8;-447h_-40;-802l_-23;-779l_-19;1066h_-25;-429l_11;701h_-35;180l_-30;-563h_6;67h_5;169h_-36;-864l_-21;-219l_-14;282h_5;567h_-31;453l_-24;740l_-35;-255l_-29;723h_-12;1050h_-24;-204h_17;303h_-12;-237l_-34;236h_-30;-408h_-32;-102l_-6;-1128l_-29;-474h_-12;1092h_-14;-491h_-35;889h_-35;45l_-8;1182l_-40;477l_-28;-172l_-13;1060h_-13;221h_-26;157h_-31;1091l_-17;-1021l_-2;846h_-33;1014h_-36;964h_-38;   162

        
        line = line.strip().split() 
        path_id=int(line[0])
        abundance=int(line[2])

        toprint, _, _, full_seq = kc.line2seq(line[1], sequences, False, hamming_max)
        if toprint:
            # print(f">path_{path_id}|{abundance}|{header}|{bubble_facts_position_start_stops}\n{full_seq}")# TMI
            print(f">path_{path_id}|{abundance}\n{full_seq}")
        else: 
            nb_non_writen+=1
            non_writen_list+=f"path_{path_id}|{abundance}, "
            
    if nb_non_writen>0:
        sys.stderr.write(f"Warning, {nb_non_writen} facts were removed as their sequence concatenation were not coherent or because they contained non coherent predictions\n")
        sys.stderr.write(f"Removed paths: {non_writen_list[:-2]}\n")
    mfile.close()





def main():
    '''
    Creation of a FA file from a compacted fact int file. 
    '''
    sequences=kc.index_sequences(sys.argv[1]) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    sequences=kc.index_sequences(sys.argv[2], sequences) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    generate_sequence_paths(sequences, sys.argv[3])
    



if __name__ == "__main__":
     main()
