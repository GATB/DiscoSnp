import networkx as nx

import sys
# import getopt
import K3000_common as kc
# import sorted_list
import argparse

import negative_binomial as nb

def store_graph(gfa_file_name):
    DG = nx.DiGraph()
    gfa_file = open(gfa_file_name)

    while(True):
        gfa_line=gfa_file.readline()
        if not gfa_line: break
        if gfa_line[0]=='H': continue       #Header
        if gfa_line[0]=='S':                #node
            split_gfa_line=gfa_line.strip().split()
            # print("node", split_gfa_line)
            node_id = split_gfa_line[1]
            node_coverage = split_gfa_line[-1].split(':')[-1] #  RC:i:24 -> 24
            DG.add_node(node_id, coverage=node_coverage)
            # print("added", node_id, node_coverage)
        
        if gfa_line[0]=='L':                #edge
            split_gfa_line=gfa_line.strip().split()
            # print("edge", split_gfa_line)
            source = split_gfa_line[1]
            target = split_gfa_line[3]
            overlap = split_gfa_line[5]
            if overlap == "-2M": continue   #forbiden link
            DG.add_edge(source, target,type=overlap)
    gfa_file.close()
    return DG
    
    
def remove_outsider_nodes_from_cc(DG,cc):
    """ for each CC compute the r,p values of the coverage distributions of the nodes
    Remove nodes that are too far from the mean distribution
    This creates new CC
    """
    return                                  # does nothing for now, waiting for amin.
    if len(cc)<8: return                    # TODO: parameter
    coverages = []
    for node_id in cc: 
        print(DG.nodes[node_id]['coverage'])
        coverages.append(int(DG.nodes[node_id]['coverage']))
    print("for",coverages)
    n,p=nb.neg_bin_fit(coverages)
    E=n/p
    V=n*(1-p)/(p*p)
    print("n",n,"p",p,"E",E,"V",V)
    #TODO remove outsider nodes from the graph
    
def assign_cc(DG,max_cc_size):
    CC=list(nx.connected_components(nx.Graph(DG)))
    for cc in CC:
        remove_outsider_nodes_from_cc(DG,cc)
    
    # recompute CC after having removed outsiders from original CCs
    # assign each node to its cc_id
    CC=list(nx.connected_components(nx.Graph(DG)))
    cc_id=0
    for cc in CC:
        cc_id +=1
        if len(cc) > max_cc_size:
            # print ("remove nodes from too large", len(cc) )
            for node_id in cc:
                DG.remove_node(node_id)
            continue
        for node_id in cc:
            DG.nodes[node_id]['cc_id'] = cc_id


def remove_cc_with_cycles(DG):
    # remove pairend links and unitig links (unoriented)
    edges_to_remove = []
    for edge in DG.edges.data():
        if edge[2]['type'] == '-1M' or edge[2]['type'] == '0M' :
            edges_to_remove.append(edge)
            
    for edge in edges_to_remove:
        DG.remove_edge(edge[0],edge[1])
    cycles = list(nx.simple_cycles(DG))
    # print("cycles", cycles)
    
    G=nx.Graph(DG)
    for nodes in cycles:
        first_node_in_cycle = nodes[0]              # get the first node, the other ones in the cycle are in the same CC
        # remove the whole CC:
        CC_with_cycle = nx.node_connected_component(G, first_node_in_cycle)
        for node_id in CC_with_cycle:
            if node_id in DG.nodes():
                DG.remove_node(node_id)
    

def main():
    '''
    Compaction of set of super reads coded as set of ids of unitigs
    '''

    parser = argparse.ArgumentParser(description='Compaction of set of super reads coded as set of ids of unitigs.')
    parser.add_argument("input_file", type=str,
                        help="gfa file " )


    max_cc_size = 1000
    args = parser.parse_args()
    gfa_file_name = str(args.input_file)
    DG = store_graph(gfa_file_name)
    assign_cc(DG,max_cc_size)
    # print(DG.nodes.data())
    remove_cc_with_cycles(DG)
    
    if '1' in DG.nodes(): 
        print("is in DG", DG.nodes['1']['cc_id'])
    
if __name__ == "__main__":
     main()
