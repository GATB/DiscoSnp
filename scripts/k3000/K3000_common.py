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



def get_reverse_sr(x):
    ''' reverse of a super read x. Example is super read x = [4,2,6], reverse(x) = [-6,-2,-4] '''
    return [-b for b in x][::-1]



def is_palindromic(x):
    ''' return true is a sr x is palindromic, eg [1,2,-2,-1]'''
    if len(x)%2 == 1: return False
    for i in range(0,int(len(x)/2)):
        if x[i] != -x[-i-1]: return False
    return True
    

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
                variant=variant.split('_')[0] #remove distance to previous alleles
                facttab.append(int(f(variant)))
        
            
            facttab = get_canonical(facttab) # store the canonical version of the fact. Btw, afterwards we add all reverse complements. 
            sl.add(facttab)
    sl.unique() # Remove redundancies
    return sl

def generate_SR(file_name):
    ''' Given an input file storing super reads, store them in the SR array'''
    # -10021;68561;-86758;27414;
    sr_file = open(file_name, 'r')
    sl = sorted_list.sorted_list()
    for line in sr_file:
        if line[0]==">": continue # compatible with fasta-file format
        line = line.rstrip()[:-1].split(';')
        sr=[]
        for unitig_id in line:
            sr_val=int(unitig_id)
            sr=sr+[sr_val]
        sl.add(sr)
    return sl



def add_reverse_SR(SR):
    ''' For all super reads in SR, we add there reverse in SR
    This double the SR size, unless there are palindromes ([1,-1] for instance). Those are not added.
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
            if other[pos] != x[pos+starting_position]:          # "non colinear"
                return False
            pos+=1


    return True


def is_canonical(sr):
    ''' return True if the canonical representation of sr is itself'''
    if sr>get_reverse_sr(sr):    
        return True
    else:
        return False

def get_canonical(sr):
    ''' return the canonical representation of sr'''
    sr_=get_reverse_sr(sr)
    if sr>sr_:    
        return sr
    else:
        return sr_

def print_maximal_super_reads(SR):
    '''print all maximal super reads as a flat format'''
    for sr in SR.traverse():
        if is_canonical(sr) or is_palindromic(sr):
            if len(sr)==1:
                print (str(sr[0])+";")
            else:
                for unitig_id in sr:
                    print (str(unitig_id)+";", end="")
                print ()






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
    without_id_reverse_msr = get_reverse_sr(msr[:-1])                           # get the msr reverse complement
    Y=MSR.get_lists_starting_with_given_prefix(without_id_reverse_msr)          # find the reverse complement in the list.
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
