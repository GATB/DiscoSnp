#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
'''
Compaction of super reads. Non ambiguous paths from a set of super reads are compacted
Resulting super reads are called MSR for Maximal Super Reads
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
# import getopt
import K3000_common as kc
# import sorted_list
import argparse



def is_subsequence(x,y,position_suffix):
    pos_on_x = position_suffix
    pos_on_y = 0
    len_x = len(x)
    len_y = len(y)
    while True:
        if pos_on_y == len_y: 
            return True  # All y was read, it is included in x 
        if pos_on_x == len_x: 
            return False # All x was read, thus y is not included in x
        if kc.allele_value(x[pos_on_x]) == kc.allele_value(y[pos_on_y]): 
            pos_on_x += 1
            pos_on_y += 1
        else:
            pos_on_x += 1
            

def remove_y_subsequence_of_x(x_ref,SR):
    ''' remove all y that are subsequence of x
    Do not care about distances. 
    Exemple 3_0,4_1,5_10 is a subsequence of 0_0,1_12,2_13,3_13,4_12,5_123,6_1
    '''
    if len(x_ref) == 1: return # as we removed strict equalities, no read can be included in a read of size one.
    n = len(x_ref)
    # print ("x",x)

    for x in [x_ref, kc.get_reverse_sr(x_ref)]:
        for position_suffix in range(0,n):
            u = x[position_suffix]
            Y = SR.get_lists_starting_with_given_prefix([u])

            if x in Y: Y.remove(x)

            for y in Y:
                if len(y)+position_suffix <= n and is_subsequence(x,y,position_suffix):
                    SR.remove(y)
                    if not kc.is_palindromic(y): 
                        SR.remove(kc.get_reverse_sr(y))

def remove_strict_inclusions(SR):
    ''' remove all super reads strictly included in any other '''
    n = len(SR)
    # to_remove=[False for i in range(n)]
    checked = 0

    for sr in SR.traverse():
        if checked%100 == 0: sys.stderr.write("      Removing inclusions, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\r")
        checked += 1
        remove_y_subsequence_of_x(sr,SR)

    sys.stderr.write("      Removing inclusions, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\n")
    return SR





def right_unique_extention(SR,sr):#, unitig_lengths,k,min_conflict_overlap):
    ''' return the unique  possible right sr extension with the largest overlap
        return None if no right extensions or if more than one possible non colinear extensions
        The returned extension can be sr itself
    '''

    n = len(sr)
    #  **** Get the largest right overlap with sr ****
    for len_u in range(n-1,0,-1):
        u = sr[-len_u:]
        Y = SR.get_lists_starting_with_given_prefix(u)
        if len(Y) == 0: continue              # No y starting with u
        if len(Y) > 1: return None,len_u      # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus not right extensible.
        y = Y[0]                              # We found the largest y right overlapping sr.


        # **** check that all other right extensions are collinear with y.
        # get all y' such that LCSP(sr,y') in [1,len(u)-1]
        # get all y' starting with a suffix of u
        Y = []
        starting_positions = []
        for starting_suffix_position in range(1,len_u):
            suffix_u = u[starting_suffix_position:]

            others = SR.get_lists_starting_with_given_prefix(suffix_u)
            if len(others) >0:
                Y += others
                starting_positions += [starting_suffix_position for zz in range(len(others))]

        if len(starting_positions)>0 and not kc.colinear(y,Y,starting_positions): return None,len_u     # more than a unique right extention for x.
        return y,len_u                                                                                  # A unique maximal right extention for x (unless gready: a unique largest extention, smaller possible extentions under the gready threahold)

    # not any right extention found.
    return None,None


def  fusion    (SR,x):
    '''Main function. For a given super read x, we find y that overlap x with the highest overlap, such that :
    1/ there exists no other y' right overlapping x that is not collinear with y
    2/ there exists no other x' left overlapping y that is not collinear with x
    Once done, we compact x and y, and this is finished for x.
    '''
    y,len_u = right_unique_extention(SR,x)              # Define, if exists, the unique y != x having the largest right overlap with x.
    if y == None: return 0                              # if no unique right extension, finished, x is not right extensible.
    if y == x: return 0                                 # Do not compact x with itself, else, enter an infinite loop
    y_ = kc.get_reverse_sr(y)                           # what are left extentions of y, we right extend its reverse complement.
    if y_ == x: return 0                                # Do not compact x with its own reverse complement.
    xprime_, dontcare = right_unique_extention(SR,y_)   # Define, if exists, the unique xprime_ (!= y_) having the largest right overlap with y_.
    if xprime_ == None: return 0                        # if no unique left extension of the unique right extention of x, finished, x is not right extensible.

    # assert xprime_ == kc.get_reverse_sr(x), "X "+str(x)+" xprime_ "+str(xprime_)+" Y "+str(y)+" Y_"+str(y_)+"\n"

    # ***** FUSION *****
    # 1/ remove x and its reverse complement if not palindromic
    # 2/ if y is not x (x==y is possible if x is repeated 2,2,2,2 for instance or if prefix = suffix (1,2,1 for instance)), remove y and its reverse complement if not palindromic
    # 3/ create the new xy SR and add it (sorted fashion)

    # isthere = SR.contains(debug_id_node) #DEBUG

    # 1
    SR.remove(x)
    # ### DEBUG
    # if isthere and not SR.contains(debug_id_node):
    #     sys.stderr.write("\n\n\nSUPPRESS X \n\n\n"+str(x)+" "+str(y)+"\n\n\n")
    #     sys.exit(0)
    # ### END DEBUG
        
    if not kc.is_palindromic(x):   
        SR.remove(kc.get_reverse_sr(x))                                           
        # ### DEBUG
        # if isthere and not SR.contains(debug_id_node):
        #     sys.stderr.write("\n\n\nSUPPRESS X_ \n\n\n"+str(kc.get_reverse_sr(x))+"\n\n\n")
        #     sys.exit(0)
        # ### END DEBUG
        
    # 2
    if x != y:
        SR.remove(y)
        # ### DEBUG
        # if isthere and not SR.contains(debug_id_node):
        #     sys.stderr.write("\n\n\nSUPPRESS Y \n\n\n"+str(x)+" "+str(y)+"\n\n\n")
        #     sys.exit(0)
        # ### END DEBUG
        if not kc.is_palindromic(y):
            SR.remove(kc.get_reverse_sr(y))
                # ### DEBUG
                # if isthere and not SR.contains(debug_id_node):
                #     sys.stderr.write("\n\n\nSUPPRESS Y_ \n\n\n"+str(x)+" "+str(y)+"\n\n\n")
                #     sys.exit(0)
                # ### END DEBUG
    # 3
    new = x+y[len_u:]
    SR.sorted_add(new)
    if not kc.is_palindromic(new): SR.sorted_add(kc.get_reverse_sr(new))
    
    # if isthere and not SR.contains(debug_id_node):
    #     sys.stderr.write("\n\n\n"+str(x)+" "+str(y)+" "+str(new)+"\n\n\n")
    #     sys.exit(0)
    # we made a compaction, return 1.
    return 1



def compaction(SR):
    ''' Compaction of all sr in SR
    If a sr was not compacted, it will never be compacted.
    If it was compacted, maybe it will be re-compacted later on. However, no need to re-run the fusion on this read as
     - either I could have been compacted on the left and this was done before or this will be done latter or
     - it will be right extended later: the 'new' (the compacted sr) sr is higher in the lexicographic order than the original sr (as it is longer), thus its position is higher in the SR data structure, thus it will be seen again later.
    Note that this may be not true in parallel computations.
    '''
    
    checked = 0
    compacted = 0
    n = len(SR)
    for sr in SR.traverse():
        if checked%100 == 0:
            sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\r")
        checked += 1
        witness = fusion(SR,sr)

        if witness == 1: # a fusion was done
            compacted += 1
    sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\n")
    return SR



def main():
    '''
    Compaction of set of super reads coded as set of ids of unitigs
    '''

    parser = argparse.ArgumentParser(description='Compaction of set of super reads coded as set of ids of unitigs.')
    parser.add_argument("input_file", type=str,
                        help="input file containing dbg paths as a list of unitig ids, eg. on line looks like \"-1;24;198;\"" )



    args = parser.parse_args()
    input_file = str(args.input_file)
    
    sys.stderr.write("  Load phased alleles \r")
    SR = kc.generate_SR_from_disco_pashing(input_file)
    sys.stderr.write("  Load phased alleles.  Done          - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Add reverse complements \r")
    kc.add_reverse_SR(SR)
    sys.stderr.write("  Add reverse complements. Done       - nb SR="+ str(len(SR))+"\n")
    
    sys.stderr.write("  Remove strict inclusions\r")
    SR = remove_strict_inclusions(SR)
    sys.stderr.write("  Remove strict inclusions. Done      - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Compaction of simple paths \r")
    SR = compaction(SR)
    sys.stderr.write("  Compaction of simple paths. Done    - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Remove2 strict inclusions\r")
    SR = remove_strict_inclusions(SR)
    sys.stderr.write("  Remove2 strict inclusions. Done      - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Compaction2 of simple paths \r")
    SR = compaction(SR)
    sys.stderr.write("  Compaction2 of simple paths. Done    - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Remove3 strict inclusions\r")
    SR = remove_strict_inclusions(SR)
    sys.stderr.write("  Remove3 strict inclusions. Done      - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Compaction3 of simple paths \r")
    SR = compaction(SR)
    sys.stderr.write("  Compaction3 of simple paths. Done    - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Print canonical compacted phased alleles\n")
    kc.print_maximal_super_reads(SR)

if __name__ == "__main__":
     main()
