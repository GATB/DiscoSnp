#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Compaction of super reads. Non ambiguous paths from a set of super reads are compacted
Resulting super reads are called MSR for Maximal Super Reads
Common file
@author (except for the 'unique' function) pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
import sorted_list
import os



# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
#https://stackoverflow.com/questions/3160699/python-progress-bar 
def update_progress(progress):
    barLength = 50 # Modify this to change the length of the progress bar
    status = ""
    # if isinstance(progress, int):
    #     progress = float(progress)
    # if not isinstance(progress, float):
    #     progress = 0
    #     status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100,2), status)
    sys.stderr.write(text)
    sys.stderr.flush()
    
def file_size(f):
    old_file_position = f.tell()
    f.seek(0, os.SEEK_END)
    size = f.tell()
    f.seek(old_file_position, os.SEEK_SET)
    return size

def hamming (s1, s2):
    res=0
    if len(s1) != len(s2): return 100000
    for i in range(len(s1)):
        if s1[i].upper()!=s2[i].upper() and s1[i].upper()!='N' and s2[i].upper()!='N':
            res+=1
    return res
            

def check_overlap(s1,s2):
    # print()
    # if hamming(s1,s2) >= 5:
   #      print("s1",s1)
   #      print("s2",s2,int_snp_id_d,allele_id)
    # assert hamming(s1,s2) < 5, ""+str(hamming(s1,s2))
    return  hamming(s1,s2) < 5
    # assert s1.upper() == s2.upper()
    # print()
    
    
    
def get_complement(char):
    complement = {"A" : "T", "T" : "A", "G" : "C", "C" : "G", "a" : "t", "t" : "a", "g" : "c" , "c" : "g", "N":"N"}
    return complement[char]

    
def get_reverse_complement(seq):
    s = ""
    for i in range(len(seq)):
        s = get_complement(seq[i]) + s
    return s

allele_value = lambda x: int(x.split('_')[0])
allele_values = lambda list_: [allele_value(x) for x in list_]
distance_string_value = lambda x: x.split('_')[1]

def generate_header(raw_int_facts):
    # from 204_0;201_-23;336_-85; to 102h;100l;168h;
    res=""
    for raw_int_fact in raw_int_facts.strip(";").split(';'):
        res+=unitig_id2snp_id(allele_value(raw_int_fact))+";"
    return res

def d_list_equal(a_d,b_d):
    a=[allele_value(x) for x in a_d]
    b=[allele_value(x) for x in b_d]
    return a==b
    
def d_list_sup(a_d,b_d):
    a=[allele_value(x) for x in a_d]
    b=[allele_value(x) for x in b_d]
    return a>b
    
    
def d_list_order(a_d,b_d):
    a=[allele_value(x) for x in a_d]
    b=[allele_value(x) for x in b_d]
    if a<b: return -1
    if a==b: return 0
    return 1


def get_reverse_sr(x):
    ''' reverse of a super read x. Example is super read x = ["4_0","2_3","-6_21"], reverse(x) = ["6_0","-2_21","-4_3"] '''
    res =[]
    for i in range (len(x)):
        value = -allele_value(x[-i-1])                          # With i=0, '-6_21' is transformed into 6
        if i==0: 
            value=str(value)+"_0"                               # With i=0, add '_0' to the value
        else:
            value=str(value)+"_"+distance_string_value(x[-i])   # With i=1, add '_21' to the value
        res.append(value)
    return res


def is_palindromic(x):
    ''' return true is a sr x is palindromic, eg [1_0,2_12,-2_13,-1_12]'''
    if len(x)%2 == 1: return False
    for i in range(0,int(len(x)/2)):
        if allele_value(x[i]) != -allele_value(x[-i-1]): return False
    return True
    #
#
# print(get_reverse_sr(["4_0","2_3","-6_21"]))
# print(is_palindromic(["4_0","4_21"]))

def f(variant):
    ''' 
    sVp 
    * s='-' or nothing
    * V=int value
    * p='h' or 'l' (path)

    f(sVp)=s2V+g(p) with g(p)=0 if p='h' and 1 if p='l'
    '''
    s=''
    if variant[0]=='-': 
        s='-'
        V=int(variant[1:-1])
    else: 
        V=int(variant[:-1])
    p=variant[-1]
    odd=0
    if p=='l':
        odd=1
    
    res=(V*2)+odd
    return s+str(res)


        
    

def generate_SR_from_disco_pashing(file_name):
    mfile = open(file_name)
    sl = sorted_list.sorted_list()
    for line in mfile: 
        #9h_0;35100h_34;-42157l_33; -16792l_0;-41270h_70; => 1
        # or
        #-2586h_0;19690h_40; => 2
        if line[0]=='#': 
            # print(line,end='')
            continue
        line=line.strip().split("=>")[0]
        line=line.strip().split(" ")
        for fact in line: 
            facttab=[]
            for variant in fact.split(';')[:-1]:
                facttab.append(f(variant.split('_')[0])+"_"+variant.split('_')[1])
        
            
            facttab = get_canonical(facttab) # store the canonical version of the fact. Btw, afterwards we add all reverse complements. 
            sl.add(facttab)

    sl.unique() # Remove redundancies

    return sl

def generate_SR(file_name):
    ''' Given an input file storing super reads, store them in the SR array'''
    # -10021_0;68561_21;-86758_3;27414_12;
    sr_file = open(file_name, 'r')
    sl = sorted_list.sorted_list()
    for line in sr_file:
        if line[0]==">": continue # compatible with fasta-file format
        line = line.rstrip()[:-1].split(';')
        sr=[]
        for allele_id in line:
            # sr_val=int(unitig_id)
            sr=sr+[allele_id]
        sl.add(sr)
    return sl



def add_reverse_SR(SR):
    ''' For all super reads in SR, we add there reverse in SR
    This double the SR size, unless there are palindromes ([1_0,-1_21] for instance). Those are not added.
    We don't check that this does not create any duplicates'''
    SR_ = sorted_list.sorted_list()
    for sr in SR.traverse():
        if not is_palindromic(sr):
            SR_.add(get_reverse_sr(sr))
    for sr in SR_.traverse():
        SR.add(sr)
    return SR



def colinear(x,X,starting_positions):
    ''' Check that all sr in X are colinear with x
    For each sr in X, one knows its starting position on x, with table starting_positions'''
    for i in range(len(X)):
        other = X[i]
        starting_position = starting_positions[i]
        pos=0
        while True:
            if pos>=len(other) or pos+starting_position>=len(x) : break
            if allele_value(other[pos]) != allele_value(x[pos+starting_position]):          # "non colinear"
                return False
            pos+=1
    return True


def is_canonical(sr):
    ''' return True if the canonical representation of sr is itself'''
    if d_list_sup(sr, get_reverse_sr(sr)):
        return True
    else:
        return False

def get_canonical(sr):
    ''' return the canonical representation of sr'''
    sr_=get_reverse_sr(sr)
    if d_list_sup(sr, sr_):    
        return sr
    else:
        return sr_

def print_maximal_super_reads(SR):
    '''print all maximal super reads as a flat format'''
    for sr in SR.traverse():
        if is_canonical(sr) or is_palindromic(sr):
            if len(sr)==1:
                print (str(allele_value(sr[0]))+";")
            else:
                for unitig_id in sr:
                    print (str(unitig_id)+";", end="")
                print ()



def determine_k(fa_file_name):
    """ given the output disco file, ie discoRes_k_31_c_2_D_0_P_3_b_2_coherent.fa return the k value (ie 31). 
    """
    return int(fa_file_name.split("_k_")[1].split("_")[0]) 

def unitig_id2snp_id(unitig_id):
    sign=1
    if unitig_id<0: sign=-1
    unitig_id=abs(unitig_id)
    res=str(sign*(unitig_id//2))
    if unitig_id%2==0:
        res+="h"
    else:
        res+="l"
    return res


def get_msr_id(msr):
    ''' returns the id of a msr
    WARNING: here msr contains as last value its unique id.
    '''
    return int(msr[-1].split("_")[1])

def get_reverse_msr_id(msr,MSR):
    ''' returns the id of the reverse complement of msr
    1/ find the sequence of the msr_ in SR
    2/ grab its id
    WARNING: here msr contains as last value its unique id.
    '''
    #1/
    # print("reverse msr id of",msr)
    without_id_reverse_msr = get_reverse_sr(msr[:-1])                           # get the msr reverse complement
    # print("rc is", without_id_reverse_msr)
    Y=MSR.get_lists_starting_with_given_prefix(without_id_reverse_msr)          # find the reverse complement in the list.
    # print("Y is", Y)
    for y in Y:                                                                 # should be of size >= 1. One of them is exactly 'without_id_reverse_msr' plus its id.
        if len(y) == len(without_id_reverse_msr)+1:                             # 'y' is 'without_id_reverse_msr' with its node id
            return get_msr_id(y)                                                # 2/
    return None                                                                 # Should not happend

# SR=generate_SR("test.txt")
# SR.unique()
# SR=add_reverse_SR(SR)
# SR.unique()
# for sr in SR.traverse():
#     print (sr)
