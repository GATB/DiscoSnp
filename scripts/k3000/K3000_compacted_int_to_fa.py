#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Creation of a FA file from a compacted fact int file. 
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
import K3000_common as kc


def index_sequences(fa_file_name):
    sequences = {}
    mfile = open(fa_file_name)
    while True: 
        line1 = mfile.readline()
        if not line1: break
        line1 = line1.strip()
        line2 = mfile.readline().strip()
        line3 = mfile.readline().strip()
        line4 = mfile   .readline().strip()
        
        if not line1.startswith(">SNP"): continue
        
        #line1: 
        #>SNP_higher_path_9|P_1:30_A/C|high|nb_pol_1|left_unitig_length_152|right_unitig_length_3|C1_25|Q1_63|G1_0/1:399,14,359|rank_0
            # key is 9h
            # one need to store the left unitig length. Here 152
            # one need to store the right unitig length. Here 3
            # note that the position of the left_unitig_length field is always the same with or without multiple snps.
        line1 = line1.split('|')
        snp_id = line1[0].split('_')[-1] # from SNP_higher_path_9 to 9
        
        left_unitig_len = int(line1[4].split('_')[-1])
        right_unitig_len = int(line1[5].split('_')[-1])
        
       
        sequences[snp_id] = [left_unitig_len, right_unitig_len, line2, line4] #sequences[snp_id] = [left_unitig_len, right_unitig_len, upperseq, lowerseq] 
        
    mfile.close()
    return sequences
    
def generate_sequence_paths(sequences, k, compacted_fact_file_name):

    mfile = open(compacted_fact_file_name)
    
    for line in mfile: 
        # 38772_0;-21479_1;27388_3;-494_28;-45551_36;-11894_10;-50927_7;-66981_10;29405_22;34837_1;20095_5;
        print(">"+line.strip())
        line=line.split(';')
        previous_snp_to_end=0
        
 

        full_seq = ""
        for i,int_snp_id_d in enumerate(line[:-1]): 
            
            ## _________--------X---------______________________________      previous snp (or full sequence we dont care)
            ##                   <------- previous_snp_to_end --------->
            ##                   <- k-1 -><------ right_unitig_len ---->      previous snp  -> previous_snp_to_end = k-1+right_unitig_len of the previous snp
            ##              __--------X---------_________________________________  new SNP
            ##                   <---->  shift between snps             <------->  to be written
            
            
            ##              __--------X---------_________________________________ (new SNP)
            ##                   <---------------------------------------------->   shift between snps + k-1 + right_unitig_len 
            ##                                   - (minus)
            ##                   <------------------------------------->            previous_snp_to_end + k-1
            ## _________--------X---------______________________________      previous snp 
            ##  (k-1)+ight_unitig_len new caracters to write -all already written
            ##  to be written = (shift between snps) + k-1 + right_unitig_len - previous_snp_to_end - (k-1)
            ## If a SNP is reversed, we reverse complement the sequence and change "right“ unitig for "left" unitig
        
            allele_id = kc.unitig_id2snp_id(kc.allele_value(int_snp_id_d))
            snp_id = allele_id[:-1]

            higher=True
            if allele_id[-1] == 'l': 
                higher=False
            
            forward=True
            if allele_id[0] == '-':
                forward=False
                snp_id = snp_id[1:]
            
            if higher: seq = sequences[snp_id][2]
            else: seq = sequences[snp_id][3]
            if forward: 
                start_to_snp = sequences[snp_id][0]+k-1
                snp_to_stop = sequences[snp_id][1]+k-1
            else: 
                seq=kc.get_reverse_complement(seq)
                start_to_snp = sequences[snp_id][1]+k-1
                snp_to_stop = sequences[snp_id][0]+k-1
            
            #treat first snp apart
            if i==0: 
                full_seq+=seq
                
            else:
                to_be_written = int(kc.distance_string_value(int_snp_id_d))+ snp_to_stop - previous_snp_to_end
                if to_be_written>0:                                 # an overelap exists
                    full_seq+=seq[-to_be_written:]
                else:                                               # no overlap
                    for v in range(-to_be_written):
                        full_seq+="N"
                # print(seq[-to_be_written:], end='')
            # print (i,int_snp_id_d,seq)
            previous_snp_to_end = snp_to_stop
            
        print(full_seq)
            
            
            
        
    mdile.close()

def main():
    '''
    Creation of a FA file from a compacted fact int file. 
    '''
    sequences=index_sequences(sys.argv[1]) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    k = kc.determine_k(sys.argv[1])
    generate_sequence_paths(sequences, k, sys.argv[2])
    



if __name__ == "__main__":
     main()
