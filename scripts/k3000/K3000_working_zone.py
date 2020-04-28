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
    
def exists_edge(edges, rev_edges, node_id_1, node_id_2):
    for current_edges in edges, rev_edges: 
        if node_id_1 in current_edges: 
            if isin(current_edges[node_id_1], node_id_2): return True
        if node_id_2 in current_edges: 
            if isin(current_edges[node_id_2], node_id_1): return True
    return False
    
    

def store_ordered_edges(fact_file_name):
    """
    Store edges
    Conserve only the closest edges. In case of equalities of minimal distance edges, conserve the closest ones. 
    """
    edges = {}      # key = (signed int) node id, value = set of (signed int) target node ids with their distance (signed int). This distance may be negative as it reflects the position of the second variant wrt the end of the bubble sequence of the first. 
                    # Hence a key is a couple (signed_int, signed int).
    rev_edges = {}  # key = (signed int) node id, value = set of (signed int) target node ids (idem)
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
                    comapped_variants = used_fact.strip(";").split(";")
                    # print(comapped_variants)
                    for i in range(len(comapped_variants)-1):
                        from_node_id = comapped_variants[i].split("_")[0]
                        to_node_id   = comapped_variants[i+1].split("_")[0]
                        distance = int(comapped_variants[i+1].split("_")[1])
                        
                        # add edge:
                        if from_node_id not in edges: edges[from_node_id] = set()
                        edges[from_node_id].add((to_node_id, distance))
                        if to_node_id not in rev_edges: rev_edges[to_node_id] = set()
                        rev_edges[to_node_id].add((from_node_id, distance))
                    # print(edges)
    ## Sort each list wrt distance: 
    for current_edges in edges, rev_edges: 
        for from_variant_id in current_edges:
            current_edges[from_variant_id] = sorted(current_edges[from_variant_id], key=lambda x : x[1])
    sys.stderr.write(f"sorted edges {edges}\n")
    ## remove transitive edges:
    nb_removed = 0
    nb_total = 0
    new_edges = [{},{}]
    i=-1
    for current_edges in edges, rev_edges: 
        i+=1
        # new_edges[i] = {}
        for from_variant_id in current_edges:
            # if an ordered list contains more than 2 children (A->B and A->C for instance), remove A->C if B->C exists somewhere else
            updated_list = set()
            updated_list.add(current_edges[from_variant_id][0]) # the first one stays, whatever. 
            for id_in_current_edges in range(len(current_edges[from_variant_id])-1):
                nb_total+=1
                cur_var = current_edges[from_variant_id][id_in_current_edges]
                next_var = current_edges[from_variant_id][id_in_current_edges+1]
                if cur_var[1]==next_var[1]: # same distance to start, we keep it
                    updated_list.add(current_edges[from_variant_id][id_in_current_edges+1])
                    continue                        
                if not exists_edge(edges, rev_edges, cur_var[0], next_var[0]): # this is a not a transitive edge, we keep it
                    updated_list.add(current_edges[from_variant_id][id_in_current_edges+1])
                else:
                    sys.stderr.write(f"removed  {from_variant_id} -> {current_edges[from_variant_id][id_in_current_edges+1]}\n")
                    nb_removed+=1
            updated_list = sorted(updated_list, key=lambda x : x[1])
            # sys.stderr.write(f"avant {from_variant_id} -> {current_edges[from_variant_id]}\n")
            new_edges[i][from_variant_id]=updated_list
            # sys.stderr.write(f"apres {from_variant_id} -> {new_edges[i][from_variant_id]}\n")
    sys.stderr.write(f"Removed {nb_removed} among {nb_total}\n")
    edges = new_edges[0]
    rev_edges = new_edges[1]
    
    ## remove distances: 

    for current_edges in edges, rev_edges: 
        for from_variant_id in current_edges:
            id_only_list = set()
            # print("HO avant", current_edges[from_variant_id])
            for to_variant_id, distance in current_edges[from_variant_id]:
                id_only_list.add(to_variant_id)
            current_edges[from_variant_id] = id_only_list
            # print("HO apres", current_edges[from_variant_id])
    
    
    ## filter edges.
    # for current_edges in edges, rev_edges:
   #      for from_variant_id in current_edges:
   #          # print(from_variant_id, current_edges[from_variant_id])
   #          ## determine min dist
   #          min_dist= sys.maxsize
   #          for to_variant_id, distance in current_edges[from_variant_id]:
   #              print (from_variant_id, to_variant_id, distance)
   #              if distance < min_dist: min_dist = distance
   #          print()
   #          ## conserve only edges equal to min dist:
   #          to_variants = set()
   #          for to_variant_id, distance in current_edges[from_variant_id]:
   #              if distance == min_dist:
   #                  to_variants.add(to_variant_id)
   #          current_edges[from_variant_id] = to_variants
   #          # print(from_variant_id, current_edges[from_variant_id])
        
    
    
    return edges,rev_edges
    
    
def store_edges(fact_file_name):
    """
    Store edges
    """
    edges = {}      # key = (signed int) node id, value = set of (signed int) target node ids.
    rev_edges = {}  # key = (signed int) node id, value = set of (signed int) target node ids.
    with open(fact_file_name) as fact_file:
        for line in fact_file.readlines():
            if not line: break
            line=line.strip()
            if line[0]=="#": continue
            # print(line)
            # print(len(line.split()))
            for fact_id in range(len(line.split())-2): # In case of pairend facts, deal with the two facts (in this case len((line.split()) is 4
                comapped_variants = line.split()[fact_id].strip(";").split(";")
                for i in range(len(comapped_variants)-1):
                    from_node_id = comapped_variants[i].split("_")[0]
                    to_node_id   = comapped_variants[i+1].split("_")[0]
                    
                    # add edge:
                    if from_node_id not in edges: edges[from_node_id] = set()
                    edges[from_node_id].add(to_node_id)
                    if to_node_id not in rev_edges: rev_edges[to_node_id] = set()
                    rev_edges[to_node_id].add(from_node_id)
                    
    return edges,rev_edges
    
def remove_nodes_with_high_io_degree(edges, rev_edges, nodes, max_io=2):
    """
    remove any node that has two or more childre and remove all the corresponding children
    """
    for current_edges in edges, rev_edges:
        for from_variant_id in current_edges: 
            if len(current_edges[from_variant_id]) > max_io:
                if from_variant_id in nodes: del nodes[from_variant_id]
                for to_variant_id in current_edges[from_variant_id]:
                    if to_variant_id in nodes: del nodes[to_variant_id]
                
    
    
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
    
        
    for node_id in node_id_to_variant_id:
        # if nodes[node_id_to_variant_id[node_id]] < cov_threshold: continue
        print("S\t"+str(node_id)+"\t"+str(node_id_to_variant_id[node_id])+"\tRC:i:"+str(nodes[node_id_to_variant_id[node_id]]))
    

    for from_variant_id in edges:
        if from_variant_id not in nodes: continue # had been removed 
        # if nodes[from_variant_id] < cov_threshold: continue # too low cove
        for to_variant_id in edges[from_variant_id]:
            if to_variant_id not in nodes: continue # had been removed
            # if nodes[to_variant_id] < cov_threshold: continue # too low cove
            print("L\t"+str(variant_id_to_node_id[from_variant_id])+"\t+\t"+str(variant_id_to_node_id[to_variant_id])+"\t+\t0M")

def main():
    nodes = store_nodes(sys.argv[1])
    # print (nodes)
    edges,rev_edges = store_ordered_edges(sys.argv[1])
    remove_nodes_with_high_io_degree(edges, rev_edges, nodes, 2000)
    print_as_gfa (nodes, edges)
    

if __name__ == '__main__':
    main()