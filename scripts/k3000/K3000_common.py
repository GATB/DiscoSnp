#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Compaction of facts. Non ambiguous paths from a set of facts are compacted
@author (except for the 'unique' function) pierre peterlongo pierre.peterlongo@inria.fr
'''

__author__ = "Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr"


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

def hamming_near_perfect (s1, s2, max=0):
    """ returns True is the hamming distance between {s1} and {s2} at most equal to {max} 
    jocker N are authorized.
    s1 and s2 are compared as upper case sequences (eg a==A)
    """ 
    res=0
    
    if len(s1) != len(s2): return False
    for i in range(len(s1)):      
        if s1[i].upper() == 'N' or s2[i].upper()=='N': continue
        if s1[i].upper()!=s2[i].upper():
            res+=1
            if res>max: return False
    return True

    
    
    
def get_complement(char):
    complement = {"A" : "T", "T" : "A", "G" : "C", "C" : "G", "a" : "t", "t" : "a", "g" : "c" , "c" : "g", "N":"N"}
    return complement[char]

    
def get_reverse_complement(seq):
    s = ""
    for i in range(len(seq)):
        s = get_complement(seq[i]) + s
    return s

string_allele_value = lambda x: x.split('_')[0]
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
    return a>=b
    
    
def d_list_order(a_d,b_d):
    a=[allele_value(x) for x in a_d]
    b=[allele_value(x) for x in b_d]
    if a<b: return -1
    if a==b: return 0
    return 1


def get_reverse_fact(x):
    ''' reverse of a fact x. Example is fact x = ["4_0","2_3","-6_21"], reverse(x) = ["6_0","-2_21","-4_3"] '''
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
    ''' return true is a fact x is palindromic, eg [1_0,2_12,-2_13,-1_12]'''
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


def valid_fact(rawfact):
    """
    Checks if a fact is valid.
    A fact containing twice or more the same variant id is not valid. 
    eg: 9h_0;35100h_34;-42157l_33; ok
    eg: 10081h_0;10081l_13; ko
    eg: 10081h_0;10081h_13; ko
    eg: -10081h_0;10081l_13; ko
    """
    id_variants = set()
    for variant in rawfact.rstrip(';').split(';'):
        id_variant = variant.lstrip("-").split("_")[0][:-1] # from -10081h_0 to 10081
        if id_variant in id_variants: return False
        id_variants.add(id_variant)
    return True    

def generate_facts_from_disco_pashing(file_name):
    mfile = open(file_name)
    sl = sorted_list.sorted_list()
    nb_non_valid = 0
    for line in mfile: 
        #9h_0;35100h_34;-42157l_33; -16792l_0;-41270h_70; => 1
        # or
        #-2586h_0;19690h_40; => 2
        if line[0]=='#': 
            # print(line,end='')
            continue
        line=line.strip().split("=>")[0]
        line=line.strip().split()
        for fact in line: 
            if not valid_fact(fact): 
                nb_non_valid+=1
                continue
            facttab=[]
            for variant in fact.split(';')[:-1]:
                facttab.append(f(variant.split('_')[0])+"_"+variant.split('_')[1])
        
            
            facttab = get_canonical(facttab) # store the canonical version of the fact. Btw, afterwards we add all reverse complements. 
            sl.add(facttab)

    sl.unique() # Remove redundancies
    return sl, nb_non_valid

def generate_facts(file_name):
    ''' Given an input file storing facts, store them in the fact array'''
    # -10021_0;68561_21;-86758_3;27414_12;
    sr_file = open(file_name, 'r')
    sl = sorted_list.sorted_list()
    for line in sr_file:
        if line[0]==">": continue # compatible with fasta-file format
        line = line.split()[0].rstrip()[:-1].split(';')
        sr=[]
        for allele_id in line:
            # sr_val=int(unitig_id)
            sr=sr+[allele_id]
        sl.add(sr)
        # print(sr)#DEBUG
    return sl

def detect_input_output_edges(facts):
    '''
    Given all stored facts and ther reverse complements, stores variants that have at least one input edge (in has_input_edge) 
    and those that have at least one output edge (in has_output_edge)

    if a fact is ["4_0","2_3","-6_21"], 
        has_input_edge = {2,-6}
        has_output_edge = {4,2}
    
    
    this means that ["6_0","-2_21","-4_3"] is also stored, hence:
        has_input_edge = {2, -6, -2, -4}
        has_output_edge = {4, 2, 6, -2}



    Returns:
        has_input_edge, has_output_edge
    '''
    has_input_edge= set()
    has_output_edge = set()


    for fact in facts.traverse():
        previous_variant_id = ""
        # print(fact)
        for variant_id in allele_values(fact):
            if previous_variant_id == "": # first value, just store it:
                previous_variant_id = variant_id
                continue
            has_output_edge.add(previous_variant_id)
            has_input_edge.add(variant_id)
            previous_variant_id = variant_id
        # print(f'output {has_output_edge}')
        # print(f'inputput {has_input_edge}')
        # sys.exit()
    return has_input_edge,has_output_edge

def add_reverse_facts(facts):
    ''' For all facts, we add their reverse in the same structure 
    This double the fact size, unless there are palindromes ([1_0,-1_21] for instance). Those are not added.
    We don't check that this does not create any duplicates'''
    facts_ = sorted_list.sorted_list()
    for fact in facts.traverse():
        if not is_palindromic(fact):
            facts_.add(get_reverse_fact(fact))
    for fact in facts_.traverse():
        facts.add(fact)
    return facts



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
    if d_list_sup(sr, get_reverse_fact(sr)):
        return True
    else:
        return False

def get_canonical(fact):
    ''' return the canonical representation of sr'''
    fact_=get_reverse_fact(fact)
    if d_list_sup(fact, fact_):    
        return fact
    else:
        return fact_

def print_maximal_facts(facts):

    '''print all maximal facts as a flat format'''
    for sr in facts.traverse():
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


def get_fact_id(fact):
    ''' returns the id of a fact
    WARNING: here fact contains as last value its unique id.
    '''
    return int(fact[-1].split("_")[1])

def get_reverse_fact_id(fact,facts):
    ''' returns the id of the reverse complement of fact
    1/ find the sequence of the fact_ in facts
    2/ grab its id
    WARNING: here fact contains as last value its unique id.
    '''
    #1/
    # print("reverse fact id of",fact)
    without_id_reverse_fact = get_reverse_fact(fact[:-1])                           # get the fact reverse complement
    # print("rc is", without_id_reverse_fact)
    Y=facts.get_lists_starting_with_given_prefix(without_id_reverse_fact)          # find the reverse complement in the list.
    # print("Y is", Y)
    for y in Y:                                                                 # should be of size >= 1. One of them is exactly 'without_id_reverse_fact' plus its id.
        if len(y) == len(without_id_reverse_fact)+1:                             # 'y' is 'without_id_reverse_fact' with its node id
            return get_fact_id(y)                                                # 2/
    return None                                                                 # Should not happend

# facts=generate_facts("test.txt")
# facts.unique()
# facts=add_reverse_facts(facts)
# facts.unique()
# for sr in facts.traverse():
#     print (sr)

def get_snp_positions(comment_line):
    """ 
    from a comment line, eg:
    >SNP_higher_path_99|P_1:20_C/T,P_2:25_T/C|high|nb_pol_2|left_unitig_length_21|right_unitig_length_14|C1_262|Q1_73|G1_0/1:5710,65,3395|rank_0
    returns the snp positions, eg
    [20, 25]
    """
    positions = []
    for variants_localization in comment_line.split("|")[1].split(','):
        positions.append(int(variants_localization.split(":")[1].split("_")[0]))
    return positions


def index_sequences(fa_file_name, sequences={}):
    mfile = open(fa_file_name)
    while True:
        line1 = mfile.readline()
        if not line1: break
        line1 = line1.strip()
        line2 = mfile.readline().strip()
        mfile.readline()              # USELESS
        line4 = mfile.readline().strip()

        if not line1.startswith(">SNP"): continue

        #line1: 
        #>SNP_higher_path_9|P_1:30_A/C|high|nb_pol_1|left_unitig_length_152|right_unitig_length_3|C1_25|Q1_63|G1_0/1:399,14,359|rank_0
            # key is 9h
            # one need to store the left unitig length. Here 152
            # one need to store the right unitig length. Here 3
            # note that the position of the left_unitig_length field is always the same with or without multiple snps.
        snp_positions = get_snp_positions(line1)
        line1 = line1.split('|')
        snp_id = line1[0].split('_')[-1]  # from SNP_higher_path_9 to 9

        left_unitig_len = int(line1[4].split('_')[-1])
        right_unitig_len = int(line1[5].split('_')[-1])

        sequences[snp_id] = [left_unitig_len, right_unitig_len, line2, line4, snp_positions]  #sequences[snp_id] = [left_unitig_len, right_unitig_len, upperseq, lowerseq, snp_positions] 
        
    mfile.close()
    return sequences

def seq_to_lower_case_except_SNPs(seq, snp_positions):
    """ given a sequence in ACGTN, put everything in lower case, except the variant positions on the upper case bubble"""
    res=""
    position_on_bubble=-1
    for letter in seq:
        if letter.isupper():
            position_on_bubble+=1
        if position_on_bubble!=-1 and position_on_bubble in snp_positions: res+=letter.upper()
        else: res+=letter.lower()
    return res

def update_SNP_positions(seq1, seq2):
    """
    Given a sequence `seq1` with upper case letters showing SNP positions, eg
    acggcgagTg
    update SNPs positions given a a sequence `seq2` of the same length, eg: 
    aGggTgaggg    
    final results is 
    aGggTgagTg
    In case of equality of lower case letter, seq1 is abitrary showen
    raises an error if two distinct upper case letters occur at the same position
    """
    new_seq = ""
    # print(f'Try concatenate \n{seq1} with \n{seq2}')

    for i, letter in enumerate(seq1):
        if seq2[i].isupper():
            # detects errors
            if letter.isupper() and letter != seq2[i]: raise ValueError (f'Cannot concatenate \n{seq1} with \n{seq2}, no solution for position {i} of suffix')
            # no error conserve the previous value as it was upper case
            new_seq+=seq2[i]
        # here either letter is upper, we keep it, or letter is lower, we also keept it. 
        else: new_seq+=letter
    return new_seq



def test_update_SNP_positions():
    seq1 = "atggcgagTg"
    seq2 = "aGggTgtggg"
    res = update_SNP_positions(seq1,seq2)
    assert res == "aGggTgagTg"
    
# test_update_SNP_positions()
    
def line2seq(line, sequences, int_facts_format, hamming_max=3):
    '''
    Parses a (non paired) fact, represented by ints or not
    Returns a bench of information relative to the line to be printed or not to be printed
    '''
    header = line.strip()+ " SP:"  # add latter the starting and ending positions of each allele on the global sequence (SP = Sequence positions). Enables to recover the good overlap length in the final GFA file
    bubble_facts_position_start_stops = "BP:" # to the header is also added the Bubble positions. For each allele in the fact we store the distance between the bubble start (upper case letter and the end of the previous bubble (also upper case letter). We add the length of the bubble (upper case letter).
    # EG:
    # ------XXXXXXXXXXXXXXXXXX------  0_18
    #                ------------XXXXXXXXXXXXXX----------- 3:14
    #                         ----------XXXXXXXXXXXXXXXX----------   -7:16 (distance is negative is bubbles overlap 

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
        if int_facts_format:
            allele_id = unitig_id2snp_id(allele_value(int_snp_id_d))
        else:
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
                len_bubble = len(seq)-lu-ru # len sequence - len left unitig - len right unitig
                snp_positions = sequences[snp_id][4]
            else: 
                seq=get_reverse_complement(seq)
                lu = sequences[snp_id][1]
                ru = sequences[snp_id][0]
                len_bubble = len(seq)-lu-ru # len sequence - len left unitig - len right unitig
                snp_positions = [len_bubble-x-1 for x in sequences[snp_id][4]] # reverse snp positions
                snp_positions.sort() # and sort them in increasing order


            #conserves only SNP positions in upper case: 
            seq = seq_to_lower_case_except_SNPs(seq, snp_positions)
            
            #treat first snp apart
            if i==0: 
                full_seq+=seq
                previous_bubble_ru = ru
                header+="0_"+str(len(full_seq))+";" # SP==Sequence Positions
                bubble_facts_position_start_stops+="0_"+str(len_bubble)+";"
                # print("full_seq =",full_seq)
                
            else:
                to_be_written = len_bubble + ru + int(distance_string_value(int_snp_id_d)) - previous_bubble_ru
                bubble_facts_position_start_stops+=distance_string_value(int_snp_id_d)+"_"+str(len_bubble)+";"
                # #DEBUG
                # print("to_be_written =",to_be_written)
                # print("len seq =",len(seq))
                # print("previous_bubble_ru =",previous_bubble_ru)
                # print("len_bubble =", len_bubble)
                # print("start_to_end  =", len_bubble+ru)
                # print("shift =", int(distance_string_value(int_snp_id_d)))
                
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
                            
                        if not hamming_near_perfect(full_seq[-p:], seq[start_on_seq:stop_on_seq], hamming_max):    #Fake read (happens with reads containing indels). We could try to retreive the good sequence, but it'd be time consuming and useless as other reads should find the good shift.
                            toprint = False
                            break
                        header+=str(len(full_seq)-len(seq)+to_be_written)    # starting position of the new sequence on the full seq that overlaps the full seq by len(seq)-to_be_written
                        try:
                            overlap = update_SNP_positions(full_seq[-p:], seq[start_on_seq:stop_on_seq])
                        except ValueError as err: 
                            toprint = False
                            # sys.stderr.write(f"{format(err)} \n")
                            break
                            
                        full_seq = full_seq[:-p]+overlap+seq[-to_be_written:] 

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
                    previous_bubble_ru = previous_bubble_ru-int(distance_string_value(int_snp_id_d))-len_bubble
                    
                    overlap_start = len(full_seq)+to_be_written-len(seq) #(included)
                    overlap_stop = len(full_seq)+to_be_written #(excluded)
                    
                    header += "I_"+str(overlap_start)+"_"+str(overlap_stop)+";"                                          # this allele is useless we do not store its start and stop positions
                    
                    if not hamming_near_perfect(full_seq[overlap_start:overlap_stop], seq, hamming_max): 
                        # print(f"Too far:\n{full_seq[overlap_start:overlap_stop]}\n{seq}\n")
                        toprint=False
                    # print(full_seq[overlap_start:overlap_stop]+"\n"+seq+"\n")
                    try:
                        overlap = update_SNP_positions(full_seq[overlap_start:overlap_stop], seq)
                    except ValueError as err: 
                        toprint = False
                        # sys.stderr.write(f"{format(err)} \n")
                        break
                    full_seq = full_seq[:overlap_start] + overlap + full_seq[overlap_stop:]
                    
            
        except KeyError: # in case a variant is in the phasing file but absent from the disco file. This is due to uncoherent prediction
            sys.stderr.write(f"{snp_id} not in sequences\n")
            toprint=False
            break
    return toprint, header, bubble_facts_position_start_stops, full_seq