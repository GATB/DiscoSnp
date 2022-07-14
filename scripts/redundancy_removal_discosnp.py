#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###################################
# from kissnp output: 
#   check pair of variants that start with the same kmer and keep only one bubble
#   check pair of variants that end with the same kmer and keep only one bubble
import sys

def get_first_kmer(seq,k):
    for i in range (len(seq)):
        if seq[i]>='a' and seq[i]<='z': continue
        return seq[i:i+k]


def get_last_kmer(seq,k):
    variant_seen=False
    for i in range (len(seq)):
        if seq[i]>='a' and seq[i]<='z':
            if not variant_seen: continue # nothin to do, continue
            else: # first position after upper case variant:
                return seq[i-k:i]
        else: variant_seen=True
    return seq[i-k+1:]

def add_id(kmer_to_var_id,kmer,var_id):
    if kmer not in kmer_to_var_id:
        kmer_to_var_id[kmer]=[]
    kmer_to_var_id[kmer].append(var_id)
        
def non_empty_intersection(a,b):
    sa = set(a)
    sb = set(b)
    if len(sa.intersection(sb)) >0: return True
    return False
    

def parse(fafile,k,faout):
    start_kmer_to_var_id={}
    stop_kmer_to_var_id ={}
    i=1
    printed=0
    while (True):
        if i%100==0 : 
            sys.stdout.write ("\r"+str(i)+ " bubbles treated, "+str(printed)+" bubbles non redundant")
            sys.stdout.flush()
        i+=1
        com1=fafile.readline()
        if not com1: break
        com1=com1.rstrip() #>SNP_higher_path_1|P_1:30_T/G|high|nb_pol_1|left_unitig_length_22|right_unitig_length_0
        seq1=fafile.readline().rstrip()
        com2=fafile.readline().rstrip()
        seq2=fafile.readline().rstrip()
        
        # get variant_id: 
        var_id=com1.split("|")[0].split("_")[-1]
        
        # deal with starting kmer
        kmer_start1=get_first_kmer(seq1,k)
        kmer_start2=get_first_kmer(seq2,k)
        if kmer_start1 in start_kmer_to_var_id and kmer_start2 in start_kmer_to_var_id:
            list_var_id_kmer_start1 = start_kmer_to_var_id[kmer_start1]
            list_var_id_kmer_start2 = start_kmer_to_var_id[kmer_start2]
            if non_empty_intersection(list_var_id_kmer_start1,list_var_id_kmer_start2): continue # this variant has already been seen with another context
        
        # deal with ending kmer
        kmer_stop1=get_last_kmer(seq1,k)
        kmer_stop2=get_last_kmer(seq2,k)
        # print (kmer_stop1,kmer_stop2)
        if kmer_stop1 in stop_kmer_to_var_id and kmer_stop2 in stop_kmer_to_var_id:
            list_var_id_kmer_stop1 = stop_kmer_to_var_id[kmer_stop1]
            list_var_id_kmer_stop2 = stop_kmer_to_var_id[kmer_stop2]
            if non_empty_intersection(list_var_id_kmer_stop1,list_var_id_kmer_stop2): continue  # this variant has already been seen with another context
        
        add_id(start_kmer_to_var_id,kmer_start1,var_id)
        add_id(start_kmer_to_var_id,kmer_start2,var_id)
        add_id(stop_kmer_to_var_id, kmer_stop1, var_id)
        add_id(stop_kmer_to_var_id, kmer_stop2, var_id)
        faout.write(com1+"\n")
        faout.write(seq1+"\n")
        faout.write(com2+"\n")
        faout.write(seq2+"\n")
        printed+=1
    sys.stdout.write ("\r"+str(i)+ " bubbles treated, "+str(printed)+" bubbles non redundant\n")

if len(sys.argv) !=4: 
    print ("Script "+ sys.argv[0].split("/")[-1])
    print ("  From a discoSnp .fa output: from all variants that start or stop with the same kmer, keep only one of their occurrences")
    print ("  usage: "+ sys.argv[0].split("/")[-1]+ " input.fa k_value output.fa")
    exit(1)
fafile=open(sys.argv[1],"r")
k=int(sys.argv[2])
faout=open(sys.argv[3],"w")
parse(fafile,k,faout)

