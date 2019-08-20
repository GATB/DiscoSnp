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
    
    
def hamming (s1, s2):
    res=0
    if len(s1) != len(s2): return 100000
    for i in range(len(s1)):
        if s1[i].upper()!=s2[i].upper() and s1[i].upper()!='N' and s2[i].upper()!='N':
            res+=1
    return res
            

def check_overlap(s1,s2, int_snp_id_d, allele_id):
    # print()
    # if hamming(s1,s2) >= 5:
   #      print("s1",s1)
   #      print("s2",s2,int_snp_id_d,allele_id)
    # assert hamming(s1,s2) < 5, ""+str(hamming(s1,s2))
    return  hamming(s1,s2) < 5
    # assert s1.upper() == s2.upper()
    # print()
    
def generate_sequence_paths(sequences, k, compacted_fact_file_name):

    mfile = open(compacted_fact_file_name)
    nb_non_writen=0
    for line in mfile: 
        # 38772_0;-21479_1;27388_3;-494_28;-45551_36;-11894_10;-50927_7;-66981_10;29405_22;34837_1;20095_5;
        header = ">"+line.strip()
        # print(line)         #DEBUG
        line=line.split(';')
        previous_bubble_ru=0
        full_seq = ""
        toprint = True 
        for i,int_snp_id_d in enumerate(line[:-1]): 
            
            # print("#################@ i #######################@",i)
            
            ## _________--------X---------______________________________      previous snp (or full sequence we dont care)
            ##                            <---- previous_bubble_ru ---->
            ##                            <------ right_unitig_len ---->      previous snp  -> previous_bubble_ru = k-1+right_unitig_len of the previous snp
            ##                                __--------X---------_________________________________  new SNP
            ##                            <---->  shift between snps    <------------------------->  to be written
            ##                                  <---upper_case---><----------------ru------------->
            
            
            ## shift + uppercase + ru = to_be_written + previous_bubble_ru
            ## ->
            ## to_be_written = shift + uppercase + ru - previous_bubble_ru
            ## If a SNP is reversed, we reverse complement the sequence and change "right“ unitig for "left" unitig
        
            allele_id = kc.unitig_id2snp_id(kc.allele_value(int_snp_id_d))
            # print(int_snp_id_d,allele_id)
            snp_id = allele_id[:-1]

            higher=True
            if allele_id[-1] == 'l': 
                higher=False
            
            forward=True
            if allele_id[0] == '-':
                forward=False
                snp_id = snp_id[1:]
            try:
                if higher: seq = sequences[snp_id][2]
                else: seq = sequences[snp_id][3]
                if forward: 
                    lu = sequences[snp_id][0]
                    ru = sequences[snp_id][1]
                else: 
                    seq=kc.get_reverse_complement(seq)
                    lu = sequences[snp_id][1]
                    ru = sequences[snp_id][0]

                len_upper_case = len(seq)-lu-ru # len sequence - len left unitig - len right unitig
                
                #treat first snp apart
                if i==0: 
                    full_seq+=seq
                    previous_bubble_ru = ru
                    # print("full_seq =",full_seq)
                    
                else:
                    to_be_written = len_upper_case + ru + int(kc.distance_string_value(int_snp_id_d)) - previous_bubble_ru
                    #DEBUG
                    # print("to_be_written =",to_be_written)
                    #
                    # print("len seq =",len(seq))
                    # print("previous_bubble_ru =",previous_bubble_ru)
                    # print("len_upper_case =", len_upper_case)
                    # print("start_to_end  =", len_upper_case+ru)
                    # print("shift =", int(kc.distance_string_value(int_snp_id_d)))
                    #
                    # print("full_seq =",full_seq)
                    # print("seq =",seq)
                    # if to_be_written>0: print("add ", seq[-to_be_written:])
                    # else: print("add nothing")
                    
                    if to_be_written>0:                                 # an overlap exists
                        if to_be_written<=len(seq):                     # the to_be_written part is smaller or equal to the length of the new sequence, we concatenate the to_be_written part.
                            # check that the overlap is correct. Fake read happen with reads containing indels. We could try to retreive the good overlap, but it'd be time consuming and useless as other reads should find the good shift.
                            p=len(seq)-to_be_written                    # maximal size of the overlap. 
                            if p > len(full_seq):
                            # if the previous sequence is included into the new one, the start on the new seq is shifted
                            #      ------------------   full_seq
                            # ------------------------- seq
                            # <---> = len(seq)-len(full_seq)-to_be_written
                                start_on_seq=len(seq)-len(full_seq)-to_be_written
                                stop_on_seq=start_on_seq+min(p,len(full_seq))
                            else:
                                start_on_seq=0
                                stop_on_seq=start_on_seq+p
                                
                            if not check_overlap(full_seq[-p:],seq[start_on_seq:stop_on_seq], int_snp_id_d, allele_id):    #Fake read (happens with reads containing indels). We could try to retreive the good sequence, but it'd be time consuming and useless as other reads should find the good shift.
                                header+=" ZZZZSHIFTED NOT ASSEMBLED"
                                toprint = False
                                break
                                # full_seq=""
                                
                            full_seq+=seq[-to_be_written:]
                        else:                                           # the to_be_written part is bigger than the length of the new sequence, we fill with Ns and add the full new seq
                            for i in range(to_be_written-len(seq)):
                                full_seq+='N'
                            full_seq+=seq
                        previous_bubble_ru = ru
                            
                    else:                                               # the new seq finishes before the already created full_seq. In this case we need to update the previous ru wrt 
                        ### ----------XXXXXXXX-------------------------- 
                        ###                   <--      pbru           -->
                        ###        -------------------XXXXX-------
                        ###                   <shift >     
                        ###                                <---npbru --->  (next previous_bubble_ru)
                        ### pbru = shift +len(upper) + npbru --> 
                        ### npbru = pbru - shift - len(upper)
                        previous_bubble_ru = previous_bubble_ru-int(kc.distance_string_value(int_snp_id_d))-len_upper_case
                
            except KeyError: # in case a variant is in the phasing file but absent from the disco file. This is due to uncoherent prediction
                header+=" UNCOHERENT"
                full_seq=""
                toprint=False
                break
        if toprint:
            print(header+"\n"+full_seq)
        else: nb_non_writen+=1
            
    if nb_non_writen>0:
        sys.stderr.write("Warning, "+str(nb_non_writen)+" facts were removed as their sequence concatenation were not coherent or because they contained non coherent predictions\n")
        
    mfile.close()

def main():
    '''
    Creation of a FA file from a compacted fact int file. 
    '''
    sequences=index_sequences(sys.argv[1]) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    k = kc.determine_k(sys.argv[1])
    generate_sequence_paths(sequences, k, sys.argv[2])
    



if __name__ == "__main__":
     main()
