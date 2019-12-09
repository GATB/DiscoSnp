import sys
import K3000_gfa_post_treatment as gpt # :)



"""
Dat file: 
AcorusCalamus_multigraph.dat

data;
param Title := AcorusCalamus;
param start := "3_0_F";
set V :=
3_0_F
1_0_R
0_0_R
0_0_F
1_0_F
1_1_F
3_0_R
1_1_R
;
param w :=
3_0_F	83691
1_0_R	26025
0_0_R	18476
0_0_F	18476
1_0_F	26025
1_1_F	26025
3_0_R	83691
1_1_R	26025
;
set reverse :=
3_0_F	3_0_R
1_0_R	1_0_F
0_0_R	0_0_F
0_0_F	0_0_R
1_0_F	1_0_R
1_1_F	1_1_R
3_0_R	3_0_F
1_1_R	1_1_F
;
set Edges :=
3_0_F	1_1_F	overlaps
3_0_F	1_1_F	links
3_0_F	1_0_F	overlaps
3_0_F	1_0_F	links
1_0_R	3_0_F	overlaps
1_0_R	3_0_F	links
1_0_R	3_0_R	overlaps
1_0_R	3_0_R	links
0_0_R	1_0_R	overlaps
0_0_R	1_0_R	links
0_0_R	1_1_R	overlaps
0_0_R	1_1_R	links
0_0_F	1_0_R	overlaps
0_0_F	1_0_R	links
0_0_F	1_1_R	overlaps
0_0_F	1_1_R	links
1_0_F	0_0_R	overlaps
1_0_F	0_0_R	links
1_0_F	0_0_F	overlaps
1_0_F	0_0_F	links
1_1_F	0_0_R	overlaps
1_1_F	0_0_R	links
1_1_F	0_0_F	overlaps
1_1_F	0_0_F	links
3_0_R	1_1_F	overlaps
3_0_R	1_1_F	links
3_0_R	1_0_F	overlaps
3_0_R	1_0_F	links
1_1_R	3_0_F	overlaps
1_1_R	3_0_F	links
1_1_R	3_0_R	overlaps
1_1_R	3_0_R	links
;
param l :=
3_0_F	1_1_F	overlaps	-99
3_0_F	1_1_F	links	0
3_0_F	1_0_F	overlaps	-99
3_0_F	1_0_F	links	0
1_0_R	3_0_F	overlaps	-99
1_0_R	3_0_F	links	0
1_0_R	3_0_R	overlaps	-99
1_0_R	3_0_R	links	0
0_0_R	1_0_R	overlaps	-99
0_0_R	1_0_R	links	0
0_0_R	1_1_R	overlaps	-99
0_0_R	1_1_R	links	0
0_0_F	1_0_R	overlaps	-99
0_0_F	1_0_R	links	0
0_0_F	1_1_R	overlaps	-99
0_0_F	1_1_R	links	0
1_0_F	0_0_R	overlaps	-99
1_0_F	0_0_R	links	0
1_0_F	0_0_F	overlaps	-99
1_0_F	0_0_F	links	0
1_1_F	0_0_R	overlaps	-99
1_1_F	0_0_R	links	0
1_1_F	0_0_F	overlaps	-99
1_1_F	0_0_F	links	0
3_0_R	1_1_F	overlaps	-99
3_0_R	1_1_F	links	0
3_0_R	1_0_F	overlaps	-99
3_0_R	1_0_F	links	0
1_1_R	3_0_F	overlaps	-99
1_1_R	3_0_F	links	0
1_1_R	3_0_R	overlaps	-99
1_1_R	3_0_R	links	0
;
"""

def print_header():
    print("data;")
    #print("param Title := ???;")
    #print("param start := \"???\";")


def print_nodes(gfa_file_name, DG=None):
    print("#id of forward nodes")
    print("set V :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        if not line: break
        line=line.strip()
        if line[0]=="S":
            "S       0       28175h;10031h;12786h;-41223l;-26670h; SP:0_426;383_541;427_586;542_661;587_731; BP:0_93;-17_61;54_61;-16_61;14_84;      FC:i:64 RC:i:21"
            node_id = line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            print("p"+node_id)  # 'p' stands for "plus strand"
    print(";")
    gfa_file.close()
    
def print_nodes_weight(gfa_file_name, DG=None):
    print("#id of forward nodes with their coverage. Here (3 sept 2019) coverage refers to the read coverage of the less covered allele of all alleles of the fact")
    print("param w :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            node_id=line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            "S       0       28175h;10031h;12786h;-41223l;-26670h; SP:0_426;383_541;427_586;542_661;587_731; BP:0_93;-17_61;54_61;-16_61;14_84;      FC:i:64 RC:i:21"
            print("p"+node_id+"\t"+line.split()[6].split(":")[-1])
    print(";")
    gfa_file.close()
    

def print_reverse(gfa_file_name,DG=None):
    print("#for each forward node, indicates the id of the reverse version")
    print("set reverse :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            "S       0       28175h;10031h;12786h;-41223l;-26670h; SP:0_426;383_541;427_586;542_661;587_731; BP:0_93;-17_61;54_61;-16_61;14_84;      FC:i:64 RC:i:21"
            node_id=line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            print("p"+node_id+"\t"+"m"+line.split()[1])  # 'p' stands for "plus strand"
    print(";")
    gfa_file.close()
    
    
def print_nodes_connected_components(gfa_file_name, DG):
    print("#id of forward nodes with their connected component id. ")
    print("param c :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            "S       0       28175h;10031h;12786h;-41223l;-26670h; SP:0_426;383_541;427_586;542_661;587_731; BP:0_93;-17_61;54_61;-16_61;14_84;      FC:i:64 RC:i:21"
            node_id=line.split()[1]
            if node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            cc_id = DG.nodes[node_id]['cc_id']
            print("p"+node_id+"\t"+str(cc_id))
            print("m"+node_id+"\t"+str(cc_id))
    print(";")
    gfa_file.close()
    
    
def print_edges(gfa_file_name, DG=None):
    print("#set of edges. Four types of edges, 1/ \"overlaps\" edges, that show an overlap between facts and 2/ \"links\" edges, that represent facts linked by paired reads (distanace unknown) and 3/ \"successive\" edges that represent two successive facts (without phasing) and 4/ \"incompatible\" edges, no path should contain two nodes linked by such an edge ")
    print("set Edges :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="L":
            "L      1       -       29384   +       8M"
            "to"
            "m1	p29384	overlaps"
            source_id = line.split()[1]
            target_id = line.split()[3]
            if DG:
                if source_id not in DG.nodes() or target_id not in DG.nodes(): continue # one of these nodes was removed during gfa post treatment
            overlap_len = int(line.split()[5].rstrip("M"))
            sign_source="p"
            if line.split()[2]=='-': sign_source="m"
            sign_target="p"
            if line.split()[4]=='-': sign_target="m"
            type="overlaps"
            if overlap_len==0: type="links"
            if overlap_len==-1: type="successive"
            if overlap_len==-2: type="incompatibles"
            print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type)  
            # do not print other edges (unitig linked and snp linked)
    print(";")
    gfa_file.close()
    
    
def print_edges_content(gfa_file_name, DG=None):
    print("#overlap length of each edge. For an \"overlaps\" edge, it indicates the number of common variants. For any other edge type (links, successive, or incompatibles), this is set to zero")
    print("param l :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="L":
            "L      1       -       29384   +       8M"
            "to"
            "m1	p29384	overlaps 8"

            source_id = line.split()[1]
            target_id = line.split()[3]
            if DG:
                if source_id not in DG.nodes() or target_id not in DG.nodes(): continue # one of these nodes was removed during gfa post treatment
            overlap_len = int(line.split()[5].rstrip("M"))
            sign_source="p"
            if line.split()[2]=='-': sign_source="m"
            sign_target="p"
            if line.split()[4]=='-': sign_target="m"
            type="overlaps"
            if overlap_len==0: type="links"
            if overlap_len==-1: type="successive"
            if overlap_len==-2: type="incompatibles"
            print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(max(0,overlap_len)))  
            # do not print other edges (unitig linked and snp linked)
    print(";")
    gfa_file.close()
    
    

def main(gfa_file_name):
    '''
    Creation of a DAT file from the graph_plus.gfa GFA file 
    Usage: 
        python ~/workspace/gatb-discosnp/scripts/k3000/K3000_gfa_to_dat.py graph_plus.gfa > graph_diploid.dat
    '''
    # Store the information as a graph. 
    # This enables 
    #   to compute connected components
    #   to remove cycles
    #   to remove too large cc
    DG = None
    if True: 
        max_cc_size=1000
        DG = gpt.store_graph(gfa_file_name)
        gpt.assign_cc(DG,max_cc_size)
        gpt.remove_cc_with_cycles(DG)
    
    print_header()
    print_nodes(gfa_file_name,DG)
    print_nodes_weight(gfa_file_name,DG)
    print_reverse(gfa_file_name,DG)
    if DG:
        print_nodes_connected_components(gfa_file_name, DG)
    print_edges(gfa_file_name,DG)
    print_edges_content(gfa_file_name,DG)



if __name__ == "__main__":
     main(sys.argv[1])
