import sys
import K3000_common as kc

"""
Given an input phased variants, select and output only "working zones"
A working zone is defined by
1/ Coverage XXX todo
2/ Topology: conserve only parts defined with X patterns. 
A node is linked to at most two other nodes.
If a node n1 is linked to two nodes n2 and n3, then n2 and n3 have at most two sons, n4 and n5
"""

def get_reverse_fact(x):
    ''' reverse of a fact x. Example  x = "4l_0;2h_3;-6h_21", reverse(x) = "6h_0;-2h_21;-4l_3"'''
    res =""
    x=x.strip(";").split(";")
    for i in range (len(x)):
        original_id = x[-i-1]
        if original_id[0]=="-": reversed_id = original_id[1:]
        else:                   reversed_id = "-"+original_id
        if i==0: 
            value=kc.string_allele_value(reversed_id)+"_0"                               # With i=0, add '_0' to the value
        else:
            value=kc.string_allele_value(reversed_id)+"_"+kc.distance_string_value(x[-i])   # With i=1, add '_21' to the value
        res+=value+";"
    return res


def store_nodes(fact_file_name):
    """
    Store nodes with their min, max and mean coverages
    """
    nodes = {}      # key= (int) node id, value variant_id
    with open(fact_file_name) as fact_file:
        for line in fact_file.readlines():
            if not line: break
            line=line.strip()
            if line[0]=="#": continue
            " -1145l_0;4780h_-60;-248h_-63;-7035l_-60;4077l_-59;-6265l_-60;-2978h_-69;4958h_-59;-4221h_-60;-1376l_-60;-4664l_-59;-2207h_-60;-6257l_-60;-5529l_-59;-5106h_-59;1310h_-59;-4514h_-59;3933h_-60;-5727h_-60;2325l_-60;-3976l_-59;-3363h_-60;-716l_-59;-3526h_-51;-5819h_-57;2927l_-59;-6146l_-60;-4374h_-53; => 1"
            coverage = int(line.split()[-1])
            for fact_id in range(len(line.split())-2): # In case of pairend facts, deal with the two facts (in this case len((line.split()) is 4
                original_fact = line.split()[fact_id]
                for fact in original_fact, get_reverse_fact(original_fact):
                    # print ("fact", fact)
                    for variant in fact.strip(";").split(";"):
                        # print(variant, fact)
                        signed_variant_id=variant.split("_")[0]
                        if signed_variant_id not in nodes: nodes[signed_variant_id] = 0
                        assert len(signed_variant_id)>0,line+"\n"+str(line.split()[fact_id])
                        nodes[signed_variant_id] += coverage
        return nodes
        
def f(tuple, el):
    return tuple[0]==el
    
def isin(list, el):
    for tuple in list: 
        if f(tuple, el): return True
    return False
    
def exists_edge(edges, node_id_1, node_id_2):

    if node_id_1 in edges: 
        if isin(edges[node_id_1], node_id_2): return True
    if node_id_2 in edges: 
        if isin(edges[node_id_2], node_id_1): return True
    return False
    
    

def store_ordered_edges(fact_file_name, nodes):
    """
    Store edges
    """
    edges = {}      # key = (signed int) node id, value = set of (signed int) target node ids with their distance (signed int). This distance may be negative as it reflects the position of the second variant wrt the end of the bubble sequence of the first. 
                    # Hence a key is a couple (signed_int, signed int).
    # rev_edges = {}  # key = (signed int) node id, value = set of (signed int) target node ids (idem)
    with open(fact_file_name) as fact_file:
        for line in fact_file.readlines():
            if not line: break
            line=line.strip()
            if line[0]=="#": continue
            # print(line)
            # print(len(line.split()))
            for fact_id in range(len(line.split())-2): # In case of pairend facts, deal with the two facts (in this case len((line.split()) is 4
                original_fact = line.split()[fact_id]
                for used_fact in original_fact, get_reverse_fact(original_fact):
                    # sys.stderr.write(f"used fact {used_fact}\n")
                    comapped_variants = used_fact.strip(";").split(";")
                    # print(comapped_variants)
                    for i in range(len(comapped_variants)-1):
                        from_node_id = comapped_variants[i].split("_")[0]
                        to_node_id   = comapped_variants[i+1].split("_")[0]
                        distance = int(comapped_variants[i+1].split("_")[1])
                        
                        # add edge:
                        if from_node_id not in edges: edges[from_node_id] = set()
                        edges[from_node_id].add((to_node_id, distance))
                        # if to_node_id not in rev_edges: rev_edges[to_node_id] = set()
                        # rev_edges[to_node_id].add((from_node_id, distance))
                    # sys.stderr.write(f"{edges}\n")
    ## Sort each list wrt distance: 

    for from_variant_id in edges:
        edges[from_variant_id] = sorted(edges[from_variant_id], key=lambda x : x[1])
    
    
    return edges
    
    
def remove_transitive_redundant(edges):
    nb_removed = 0
    nb_total = 0
    nb_removed = 0
    nb_total = 0
    # new_edges = {}
    for from_variant_id in edges:
        # if an ordered list contains more than 2 children (A->B and A->C for instance), remove A->C if B->C exists somewhere else
        updated_list = set()
        updated_list.add(edges[from_variant_id][0]) # the first one stays, whatever. 
        for id_in_edges in range(len(edges[from_variant_id])-1):
            nb_total+=1
            cur_var = edges[from_variant_id][id_in_edges]
            next_var = edges[from_variant_id][id_in_edges+1]
            if cur_var[1]==next_var[1]: # same distance to start, we keep it
                updated_list.add(edges[from_variant_id][id_in_edges+1])
                continue                        
            if not exists_edge(edges, cur_var[0], next_var[0]): # this is a not a transitive edge, we keep it
                updated_list.add(edges[from_variant_id][id_in_edges+1])
            else:
                # sys.stderr.write(f"removed  {from_variant_id} -> {edges[from_variant_id][id_in_edges+1]}\n")
                nb_removed+=1
        updated_list = sorted(updated_list, key=lambda x : x[1])
        # sys.stderr.write(f"avant {from_variant_id} -> {edges[from_variant_id]}\n")
        edges[from_variant_id]=updated_list
        # sys.stderr.write(f"apres {from_variant_id} -> {new_edges[i][from_variant_id]}\n")
    sys.stderr.write(f"Removed {nb_removed} transitive redundant edges among {nb_total}\n")
    return nb_removed
    # edges = new_edges
    
def keep_only_consecutive_edges(edges):
    nb_total = 0
    nb_removed = 0
    for from_variant_id in edges:
        # if an ordered list contains more than 2 children (A->B and A->C for instance), conserve the closest one
        updated_list = set()
        updated_list.add(edges[from_variant_id][0]) # the first one stays, whatever. 
        next_id = 0
        while True: # add all soon with the same distance
            next_id += 1
            if next_id == len(edges[from_variant_id]): break
            if edges[from_variant_id][next_id][1] == edges[from_variant_id][0][1]:
                updated_list.add(edges[from_variant_id][next_id])
            else: break
        nb_total += len(edges[from_variant_id])
        nb_removed += len(edges[from_variant_id])-next_id
        updated_list = sorted(updated_list, key=lambda x : x[1])
        # sys.stderr.write(f"avant {from_variant_id} -> {edges[from_variant_id]}\n")
        edges[from_variant_id]=updated_list
        # sys.stderr.write(f"apres {from_variant_id} -> {new_edges[i][from_variant_id]}\n")
    sys.stderr.write(f"Removed {nb_removed} among {nb_total}\n")
        
        
  
def remove_distances(edges):
    ## remove distances: 
    for from_variant_id in edges:
        id_only_list = set()
        for to_variant_id, distance in edges[from_variant_id]:
            id_only_list.add(to_variant_id)
        edges[from_variant_id] = id_only_list
        
    
def remove_tips(edges, nodes):
    iodegree = {}   # key = signed variant id, value = (At least one incoming edge, at least one ougoing edge)
    ## Stores for all variants, if they have at least one incoming edge and at least one ougoing edge
    for source_id, target_ids in edges.items():
        if source_id not in nodes: continue                                         # were removed
        if source_id not in iodegree:       iodegree[source_id] = [False,True]
        else:                               iodegree[source_id][1] = True
        for target_id in target_ids:
            if target_id not in nodes: continue                                     # were removed
            if target_id not in iodegree:   iodegree[target_id] = [True,False]
            else:                           iodegree[target_id][0] = True
    
    nb_tips_removed = 0
    to_remove = set()
    for node in nodes: 
        if node not in iodegree: 
            to_remove.add(node)
            nb_tips_removed+=1
            continue
        if iodegree[node] != [True, True]:
            to_remove.add(node)
            nb_tips_removed+=1
            continue
    for node_to_remove in to_remove: 
        del nodes[node_to_remove]
    sys.stderr.write(f"Removed {nb_tips_removed} tips\n")
    return nb_tips_removed
    # sys.exit()
    
def working_zones(edges, nodes):
    """ 
    if a node has 3 or more children or 3 or more ancestors, remove all of them
    """
    removed=0
    for from_variant_id in edges: 
        if len(edges[from_variant_id]) > 2:
            for to_variant_id in edges[from_variant_id]:
                if to_variant_id in nodes: 
                    removed+=1
                    del nodes[to_variant_id]
    sys.stderr.write(f"Working zones, children removed {removed} non bi-allelic nodes\n")
    
    
    
    ### ancestors. We first have to store them
    
    input_edges = {} # key = son, value = list of ancestors
    for source_id, target_ids in edges.items():
        if source_id not in nodes : continue
        for target_id in target_ids:
            if target_id not in nodes: continue
            if target_id not in input_edges: input_edges[target_id] = set()
            input_edges[target_id].add(source_id)
            
    removed=0
    for from_variant_id in input_edges: 
        if len(input_edges[from_variant_id]) > 2:
            for to_variant_id in input_edges[from_variant_id]:
                if to_variant_id in nodes: 
                    removed+=1
                    del nodes[to_variant_id]
    sys.stderr.write(f"Working zones, ancestor removed {removed} non bi-allelic nodes\n")
            
    

def existing_fact(fact, nodes, edges):
    """
    Given a fact as: 
    -1145l_0;4780h_-60;-248h_-63;-7035l_-60;4077l_-59;-6265l_-60;-2978h_-69;4958h_-59;-4221h_-60;-1376l_-60;
    returns True if 
    1/ all nodes exist in the gfa
    2/ all edges exist in the gfa
    else return False
    """
    
    s_fact = fact.strip(";").split(";")
    
    ### Nodes
    for shift_node in s_fact:
        
        if shift_node.split('_')[0] not in nodes: 
            # sys.stderr.write(f"{shift_node.split('_')[0]} not in nodes\n")
            return False
        
    ### Edges
    for i in range(len(s_fact)-1):
        source = s_fact[i].split('_')[0]
        target = s_fact[i+1].split('_')[0]
        if source in edges and target in edges[source]: continue
        if target in edges and source in edges[target]: continue
        return False
    return True

def print_as_facts(nodes, edges, fact_file_name):
    """
    Store nodes with their min, max and mean coverages
    """
    with open(fact_file_name) as fact_file:
        for line in fact_file.readlines():
            if not line: break
            line=line.strip()
            if line[0]=="#": continue
            " -1145l_0;4780h_-60;-248h_-63;-7035l_-60;4077l_-59;-6265l_-60;-2978h_-69;4958h_-59;-4221h_-60;-1376l_-60;-4664l_-59;-2207h_-60;-6257l_-60;-5529l_-59;-5106h_-59;1310h_-59;-4514h_-59;3933h_-60;-5727h_-60;2325l_-60;-3976l_-59;-3363h_-60;-716l_-59;-3526h_-51;-5819h_-57;2927l_-59;-6146l_-60;-4374h_-53; => 1"
            coverage = int(line.split()[-1])
            printed_something = False
            for fact_id in range(len(line.split())-2): # In case of pairend facts, deal with the two facts (in this case len((line.split()) is 4
                original_fact = line.split()[fact_id]
                if existing_fact(original_fact, nodes, edges):
                    print(original_fact, ' ', end='')
                    printed_something = True
            if printed_something:
                print ("=> ", line.split()[-1])
            
    

    
def print_as_gfa(nodes, edges):
    # cov_threshold=10
    # index nodes:
    node_id_to_variant_id = {} # key = node_id (1,2, ..., n) value = variant_id (eg. -14h)
    variant_id_to_node_id = {} # key = variant_id (eg. -14h) value = node_id (eg. 2)
    node_id=1
    for node in nodes:
        node_id_to_variant_id[node_id] = node
        variant_id_to_node_id[node] = node_id
        node_id+=1
        
    
    nb_nodes = 0
    for node_id in node_id_to_variant_id:
        nb_nodes+=1
        # if nodes[node_id_to_variant_id[node_id]] < cov_threshold: continue
        print("S\t"+str(node_id)+"\t"+str(node_id_to_variant_id[node_id])+"\tRC:i:"+str(nodes[node_id_to_variant_id[node_id]]))
    
    nb_edges = 0
    for from_variant_id in edges:
        if from_variant_id not in nodes: continue # had been removed 
        # if nodes[from_variant_id] < cov_threshold: continue # too low cove
        for to_variant_id in edges[from_variant_id]:
            if to_variant_id not in nodes: continue # had been removed
            # if nodes[to_variant_id] < cov_threshold: continue # too low cove
            nb_edges+=1
            print("L\t"+str(variant_id_to_node_id[from_variant_id])+"\t+\t"+str(variant_id_to_node_id[to_variant_id])+"\t+\t0M")
    
    sys.stderr.write(f"Printed {nb_nodes} nodes and {nb_edges} edges\n")

def main():
    outgfa = False
    if len(sys.argv) > 2 and sys.argv[2]=="--gfa": 
        outgfa=True
        sys.stderr.write(f"Output data as a gfa graph\n")
    nodes = store_nodes(sys.argv[1])
        
    # print (nodes)
    edges = store_ordered_edges(sys.argv[1], nodes)
    # while remove_transitive_redundant(edges) > 0:
 #        True
    keep_only_consecutive_edges(edges)
    
    
    
    remove_distances(edges)
    remove_tips(edges, nodes)
    working_zones(edges, nodes)
    remove_tips(edges, nodes)
    if outgfa:
        print_as_gfa (nodes, edges)
    else:
        print_as_facts(nodes, edges, sys.argv[1])
    

if __name__ == '__main__':
    main()