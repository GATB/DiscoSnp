#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Creation of a GFA file from a set of compacted maximal super reads
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
import K3000_common as kc


def show_right_edges (MSR,x,id_x):#,unitigs,k):
    ''' Main function. For a given super read x, we find y that overlap x, and we print the links in a GFA style:
    L	11	+	12	-	overlap size
    Note that one treat x only if its canonical.
    Four cases :
    1/ x overlaps y, with y canonical. One prints x + y + blabla
    2/ x overlaps y_, with y_ non canonical. In this case, y_ overlaps x. One of the two solutions has to be chosen. We chose min(idx,idy) (with idx,idy being the ids of the MSR x,y in SR) One searches the id of y, and one prints x + y - blabla.
    3/ x_ overlaps y. same as 2.
    4/ x_ overlaps y_. We do nothing, this case is treated when the entry of the function is y that thus overlaps x.
    WARNING: here x and each msr in MSR contain as last value its unique id.
    '''
    x=x[:-1]                                # remove the x_id from the x msr
    if not kc.is_canonical(x): return
    n=len(x)

    # CASES 1 AND 2
    strandx='+'
    # print ("x is", x)
    for len_u in range(1,n): # for each possible x suffix
        u=x[-len_u:]
        # print ("u is", u)
        Y=MSR.get_lists_starting_with_given_prefix(u)
        # if x in Y: Y.remove(x)            # we remove x itself from the list of y : note that it should not occur.
        # print (x)
        # assert(x not in Y)
        if len(Y)==0: continue              # No y starting with u
        for y in Y:
            # detect the y strand
            # CASE 1/
            if kc.is_canonical(y[:-1]):     # remove the last value that corresponds to the node id
                strandy ='+'
                # id_y=indexed_nodes.index(y)                               # get the id of the target node
                id_y=kc.get_msr_id(y)                                           # last value is the node id, here the id of the target node
            # CASE 2/
            else:
                strandy='-'
#                id_y = indexed_nodes.index(kc.get_reverse_msr(y))
                id_y=kc.get_reverse_msr_id(y,MSR)                                # find the reverse of list y in MSR to grab its id.
                if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
            # print the edges
            print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(len_u)+"M")


    # CASES 3 AND 4
    strandx='-'
    x_=kc.get_reverse_sr(x)
    for len_u in range(1,n): # for each possible x suffix
        u=x_[-len_u:]
        Y=MSR.get_lists_starting_with_given_prefix(u)
        # assert(x_ not in Y)
        if len(Y)==0: continue  # No y starting with u
        for y in Y:
            if kc.is_canonical(y[:-1]): # CASE 3
                strandy ='+'
               # id_y=indexed_nodes.index(y)                                # get the id of the target node
                id_y=kc.get_msr_id(y)                                           # last value is the node id, here the id of the target node
                # we determine min(id_x,id_y)
                if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
                print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(len_u)+"M") # note that strand x is always '-' and strandy is always '+' in this case.

#            else: continue # CASE 4, nothing to do.




def print_GFA_edges(MSR):#,unitigs,k):
    '''print each potiential edge in GFA format. Note that each edge is printed in one unique direction, the other is implicit
    WARNING: here each msr in MSR contains as last value its unique id.
    '''
    for msr in MSR.traverse():
        x_id = kc.get_msr_id(msr)                                         # last value is the node id
        if x_id%100==0: sys.stderr.write("\t%.2f"%(100*x_id/len(MSR))+"%\r")
        show_right_edges(MSR,msr,x_id)#,unitigs,k)
    sys.stderr.write("\t100.00%\n")


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

def print_GFA_nodes_as_ids(MSR):#, unitigs, k):
    '''print canonical unitigs ids
    WARNING: here each msr in MSR contains as last value its unique id.
    '''
    for msr in MSR.traverse():
        node_id = kc.get_msr_id(msr)                        # last value is the node id
        msr = msr[:-1]                                      # remove the last value that corresponds to the node id
        if not kc.is_canonical(msr):                       continue
        print ("S\t"+str(node_id)+"\t", end="")

        for unitig_id in msr:                       
            print (unitig_id2snp_id(unitig_id)+";", end="")
        print ()

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))


def main():
    '''
    Creation of a GFA file from a set of compacted maximal super reads
    '''
    # SR=[[1,3,4],[14],[4,6],[-6,-1,-2],[4,6,7,9],[9,10,11,12],[4,6],[-13, -12, -11]]
    
    MSR=kc.generate_SR(sys.argv[1])


    kc.add_reverse_SR(MSR)
    MSR.sort()
    MSR.index_nodes()                          # This adds a final value to each sr, providing its node id.
    sys.stderr.write("Print GFA Nodes\n")
    print_GFA_nodes_as_ids(MSR)
    sys.stderr.write("Print GFA Edges\n")
    print_GFA_edges(MSR)#,unitigs,k)


if __name__ == "__main__":
     main()
