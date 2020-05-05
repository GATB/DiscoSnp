#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Creation of a GFA file from a set of compacted maximal facts
'''

__author__ = "Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr"


import sys
import K3000_common as kc


def show_right_edges (facts,x,id_x):
    ''' Main function. For a given fact x, we find y that overlap x, and we print the links in a GFA style:
    L	11	+	12	-	overlap size
    Note that one treat x only if its canonical.
    Four cases :
    1/ x overlaps y, with y canonical. One prints x + y + blabla
    2/ x overlaps y_, with y_ non canonical. In this case, y_ overlaps x. One of the two solutions has to be chosen. We chose min(idx,idy) (with idx,idy being the ids of the facts x,y in SR) One searches the id of y, and one prints x + y - blabla.
    3/ x_ overlaps y. same as 2.
    4/ x_ overlaps y_. We do nothing, this case is treated when the entry of the function is y that thus overlaps x.
    WARNING: here x and each fact in facts contain as last value its unique id.
    '''
    x=x[:-1]                                # remove the x_id from the x fact
    if not kc.is_canonical(x): return
    n=len(x)

    # CASES 1 AND 2
    strandx='+'
    # print ("x is", x)
    for len_u in range(1,n): # for each possible x suffix
        u=x[-len_u:]
        # print ("u is", u)
        Y=facts.get_lists_starting_with_given_prefix(u)
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
                id_y=kc.get_fact_id(y)                                           # last value is the node id, here the id of the target node
            # CASE 2/
            else:
                strandy='-'
#                id_y = indexed_nodes.index(kc.get_reverse_fact(y))
                id_y=kc.get_reverse_fact_id(y,facts)                                # find the reverse of list y in facts to grab its id.
                # assert kc.is_canonical(full_y[:-1]), "not canonical "+str(full_y[-1])+" full = "+str(full_y)
                # print("y is ", y)
                # print("id_y is", id_y)
                if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
            # print the edges
            print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(len_u)+"M")


    # CASES 3 AND 4
    strandx='-'
    x_=kc.get_reverse_fact(x)
    for len_u in range(1,n): # for each possible x suffix
        u=x_[-len_u:]
        Y=facts.get_lists_starting_with_given_prefix(u)
        # assert(x_ not in Y)
        if len(Y)==0: continue  # No y starting with u
        for y in Y:
            if kc.is_canonical(y[:-1]): # CASE 3
                strandy ='+'
               # id_y=indexed_nodes.index(y)                                # get the id of the target node
                id_y=kc.get_fact_id(y)                                           # last value is the node id, here the id of the target node
                # we determine min(id_x,id_y)
                if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
                print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(len_u)+"M") # note that strand x is always '-' and strandy is always '+' in this case.

#            else: continue # CASE 4, nothing to do.




def print_GFA_edges(facts):
    '''print each potiential edge in GFA format. Note that each edge is printed in one unique direction, the other is implicit
    WARNING: here each fact in facts contains as last value its unique id.
    '''
    for fact in facts.traverse():
        x_id = kc.get_fact_id(fact)                                         # last value is the node id
        if x_id%100==0: sys.stderr.write("\t%.2f"%(100*x_id/len(facts))+"%\r")
        show_right_edges(facts,fact,x_id)
    sys.stderr.write("\t100.00%\n")


def check_fact(fact, fact_int):
    # fact ['49648_0', '67994_-20', '20000_23']
    # fact_int 49648_0;67994_-20;20000_23; SP:0_166;126_261;178_444; BP:0_83;-20_72;23_61;
    for i,allele_id in enumerate(fact_int.split()[0].split(";")[:-1]):
        if fact[i] != allele_id:
            sys.stderr.write("Not corresponding fact and fact_int:\n")
            sys.stderr.write(str(fact)+"\n")
            sys.stderr.write(fact_int+"\n")
            sys.exit(0)
            
def index_nodeid_to_distance(facts, compacted_fact_int_file_name): 
    """ returns an index nodeid -> distances
    """
    compacted_fact_int_file=open(compacted_fact_int_file_name)
    nodeid_to_distance = {}
    for fact_line in compacted_fact_int_file.readlines():
        # 49648_0;67994_-20;20000_23; SP:0_166;126_261;178_444; BP:0_83;-20_72;23_61;
        s_fact_line = fact_line.strip().split()
        node_as_list = [node for node in s_fact_line[0].split(";")[:-1]]
        # print(s_fact_line[0], node_as_list)
        if not  kc.is_canonical(node_as_list):                       continue
        node_id = facts.get_node_id(node_as_list)
        assert node_id not in nodeid_to_distance, "node "+str(node_id)+" already in truc, with value "+ nodeid_to_distance[node_id]
        nodeid_to_distance[node_id]=s_fact_line[1]+"\t"+s_fact_line[2]
    compacted_fact_int_file.close() 
    return nodeid_to_distance
    

def contains_extreme_variant(fact, has_input_edge, has_output_edge):
    '''
    Returns True if a fact contains at least one extreme variant
    Else returns false
    '''
    for variant_id in kc.allele_values(fact):
        if variant_id not in has_input_edge or variant_id not in has_output_edge:
            return True
    return False
    
def print_GFA_nodes_as_ids(facts, compacted_fact_int_file_name, has_input_edge, has_output_edge):
    '''print canonical unitigs ids
    WARNING: here each fact in facts contains as last value its unique id.
    '''
    nodeid_to_distance = index_nodeid_to_distance(facts, compacted_fact_int_file_name)
    # compacted_fact_int_file=open(compacted_fact_int_file_name)
    for fact in facts.traverse():
        # fact_int=compacted_fact_int_file.readline().strip() # 49648_0;67994_-20;20000_23; SP:0_166;126_261;178_444; BP:0_83;-20_72;23_61;
        node_id = kc.get_fact_id(fact)                        # last value is the node id
        fact = fact[:-1]                                      # remove the last value that corresponds to the node id
        if not kc.is_canonical(fact):                       continue
        print ("S\t"+str(node_id)+"\t", end="")
        for unitig_id in fact:                       
            print (kc.unitig_id2snp_id(kc.allele_value(unitig_id))+";", end="")
        # check_fact(fact, fact_int)
        # assert str(node_id) in nodeid_to_distance, nodeid_to_distance
        print ("\t"+nodeid_to_distance[str(node_id)], end="")
        if contains_extreme_variant(fact, has_input_edge, has_output_edge):
            print("\tEV:1")
        else:
            print("\tEV:0")

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))


def check(facts): 
    """ debug purpose """ 
    for fact in facts.traverse(): 
        assert kc.get_reverse_fact_id(fact,facts) != None, fact
   
    
def main():
    '''
    Creation of a GFA file from a set of facts
    For eahc fact: indicates if it contains at least an extreme variant
    Field EV
     * EV:1 if it contains at least an extreme variant
     * EV:0 else
    '''
    # SR=[[1,3,4],[14],[4,6],[-6,-1,-2],[4,6,7,9],[9,10,11,12],[4,6],[-13, -12, -11]]
    
    facts=kc.generate_facts(sys.argv[1])


    kc.add_reverse_facts(facts)
    has_input_edge, has_output_edge = kc.detect_input_output_edges(facts)
    facts.sort()
    facts.index_nodes()                          # This adds a final value to each sr, providing its node id.
    # check(facts)
    sys.stderr.write("Print GFA Nodes\n")
    print_GFA_nodes_as_ids(facts, sys.argv[1], has_input_edge, has_output_edge)
    sys.stderr.write("Print GFA Edges\n")
    print_GFA_edges(facts)


if __name__ == "__main__":
     main()
