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



def index_sequences(fa_file_name, sequences = {}):
    mfile = open(fa_file_name)
    while True: 
        line1 = mfile.readline()
        if not line1: break
        line1 = line1.strip()
        line2 = mfile.readline().strip()
        mfile.readline().strip()                # USELESS
        line4 = mfile.readline().strip()
        
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
    
    
def line2seq(line, sequences, hamming_max=0):
    '''
    Parses a (non paired) fact, represented by ints or not
    Returns a bench of information relative to the line to be printed or not to be printed
    '''
    header = line.strip()+ "|SP:"  # add latter the starting and ending positions of each allele on the global sequence (SP = Sequence positions). Enables to recover the good overlap length in the final GFA file
    bubble_facts_position_start_stops = "BP:" # to the header is also added the Bubble positions. For each allele in the fact we store the distance between the bubble start (upper case letter and the end of the previous bubble (also upper case letter). We add the length of the bubble (upper case letter).
    # EG:
    # ------XXXXXXXXXXXXXXXXXX------  0_18
    #                ------------XXXXXXXXXXXXXX----------- 3:14
    #                         ----------XXXXXXXXXXXXXXXX----------   -7:16 (distance is negative is bubbles overlap 
    # print("\n NEW LINE ",line)         #DEBUG
    line=line.strip(';').split(';')
    previous_bubble_ru=0
    full_seq = ""
    toprint = True
    for i,int_snp_id_d in enumerate(line): 
        
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
        allele_id = int_snp_id_d.split("_")[0]
        # print(f"\n\n *****{int_snp_id_d}, {allele_id} ********")
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
                header+="0_"+str(len(full_seq))+";" # SP==Sequence Positions
                bubble_facts_position_start_stops+="0_"+str(len_upper_case)+";"
                # print("full_seq =",full_seq)
                
            else:
                to_be_written = len_upper_case + ru + int(kc.distance_string_value(int_snp_id_d)) - previous_bubble_ru
                bubble_facts_position_start_stops+=kc.distance_string_value(int_snp_id_d)+"_"+str(len_upper_case)+";"
                # #DEBUG
                # print("to_be_written =",to_be_written)
                # print("len seq =",len(seq))
                # print("previous_bubble_ru =",previous_bubble_ru)
                # print("len_upper_case =", len_upper_case)
                # print("start_to_end  =", len_upper_case+ru)
                # print("shift =", int(kc.distance_string_value(int_snp_id_d)))
                
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
                            
                        if not kc.hamming_near_perfect(full_seq[-p:],seq[start_on_seq:stop_on_seq], hamming_max):    #Fake read (happens with reads containing indels). We could try to retreive the good sequence, but it'd be time consuming and useless as other reads should find the good shift.
                            print(f"Too far\n{full_seq[-p:]}\n{seq[start_on_seq:stop_on_seq]}\n")
                            toprint = False
                            break
                        header+=str(len(full_seq)-len(seq)+to_be_written)    # starting position of the new sequence on the full seq that overlaps the full seq by len(seq)-to_be_written
                        full_seq+=seq[-to_be_written:] 
                        header+="_"+str(len(full_seq))+";"                           # ending position of the new sequence on the full seq. 
                    else:                                           # the to_be_written part is bigger than the length of the new sequence, we fill with Ns and add the full new seq
                        for i in range(to_be_written-len(seq)):
                            full_seq+='N'
                        header+=str(len(full_seq))                           # starting position of the new sequence on the full seq (possibly overlapping the full seq)
                        full_seq+=seq
                        header+="_"+str(len(full_seq))+";"                           # ending position of the new sequence on the full seq. 
                    previous_bubble_ru = ru
                        
                else:                                               # the new seq finishes before the already created full_seq. In this case we need to update the previous ru wrt 
                    ### ----------XXXXXXXX-------------------------- 
                    ###                   <--      pbru           -->
                    ###        -------------------XXXXX-------
                    ###                   <shift >            < tbw > (negative)
                    ###                                <---npbru --->  (next previous_bubble_ru)
                    ### pbru = shift +len(upper) + npbru --> 
                    ### npbru = pbru - shift - len(upper)
                    ### TODO BUG ? (I_ here is alway followed by a negative value in the experiments I made (june 2020))
                    previous_bubble_ru = previous_bubble_ru-int(kc.distance_string_value(int_snp_id_d))-len_upper_case
                    header += "I_"+str(len(full_seq)+to_be_written-len(seq))+"_"+str(len(full_seq)+to_be_written)+";"                                          # this allele is useless we do not store its start and stop positions
                    if not kc.hamming_near_perfect(full_seq[len(full_seq)+to_be_written-len(seq):len(full_seq)+to_be_written], seq, hamming_max): 
                        print(f"Too far:\n{full_seq[len(full_seq)+to_be_written-len(seq):len(full_seq)+to_be_written]}\n{seq}\n")
                        toprint=False
                    # print(full_seq[len(full_seq)+to_be_written-len(seq):len(full_seq)+to_be_written]+"\n"+seq+"\n")
                    
            
        except KeyError: # in case a variant is in the phasing file but absent from the disco file. This is due to uncoherent prediction
            print(f"{snp_id} not in sequences")
            toprint=False
            break
    return toprint, header, bubble_facts_position_start_stops, full_seq

    
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

        toprint, header, bubble_facts_position_start_stops, full_seq = line2seq(line[1], sequences, hamming_max)
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
    sequences=index_sequences(sys.argv[1]) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    sequences=index_sequences(sys.argv[2], sequences) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    generate_sequence_paths(sequences, sys.argv[3])
    



if __name__ == "__main__":
     main()
