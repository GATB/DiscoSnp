import sys
import K3000_gfa_post_treatment as gpt # :)



"""
#id of forward and reverse nodes
set V :=
p0
p1
p2
m0
m1
m2
;
param number_loci
p0	3
p1	2
p2	2
m0 	3
m1	2
m2	2
;
#id of nodes with their coverage. Here (3 sept 2019) coverage refers to the read coverage of the less covered allele of all alleles of the fact
param w :=
p0	5626
m0  5626
p1	10152
m1	10152
p2	1142
m2	1142
;
#for each forward node, indicates the id of the reverse version
set reverse :=
p0	m0
p1	m1
p2	m2
;
#set of edges. Four types of edges, 1/ "overlaps" edges, that show an overlap between facts and 2/ "links" edges, that represent facts linked by paired reads (distanace unknown) and 3/ "successive" edges that represent two successive facts (without phasing) and 4/ "incompatible" edges, no path should contain two nodes linked by such an edge 
set Edges :=
p0	m55	overlaps #TODO (plus tard): il faut aussi doubler les arêtes et donc générer p55	m0	overlaps (pas d'ordre à respecter au sein de chaque champ)
p0	m53	overlaps
p1	m55	overlaps
p303	p305	links
p303	p310	links
p303	p304	links
p0	p54	incompatibles
p0	p1	incompatibles
p0	p2	incompatibles
p243	p341	successive
p244	p341	successive
p320	m388	successive
;
#overlap length of each edge. For an "overlaps" edge, it indicates the number of common variants. For any other edge type (links, successive, or incompatibles), this is set to zero
param l :=
p0	m55	overlaps	3
p0	m53	overlaps	3
p1	m55	overlaps	3
p303	p305	links	0
p303	p310	links	0
p303	p304	links	0
p0	p54	incompatibles	0
p0	p1	incompatibles	0
p0	p2	incompatibles	0
p243	p341	successive	0
p244	p341	successive	0
p320	m388	successive	0
;
#Coverage of the pairend links. Note that overlap links do not have any coverage (just computed from fact overlaps). Also incompatible and successive links do not have coverage, by definition. 
param pairend_end_links_coverage :=
p0	m55	overlaps	0
p0	m53	overlaps	0
p1	m55	overlaps	0
p303	p305	links	12
p303	p310	links	3
p303	p304	links	7
p0	p54	incompatibles	0
p0	p1	incompatibles	0
p0	p2	incompatibles	0
p243	p341	successive	0
p244	p341	successive	0
p320	m388	successive	0

"""

def print_header():
    print("#haplotype number")
    print("param mult := XXX")


def print_nodes(gfa_file_name, DG=None):
    print("#id of plus (p) and minus (m) nodes")
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
            print("m"+node_id)  # 'p' stands for "minus strand"
    print(";")
    gfa_file.close()
    
def print_nodes_number_loci(gfa_file_name, DG=None):
    print("#id of nodes with their number of variants")
    print("param length :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            node_id=line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            "S       0       28175h;10031h;12786h;-41223l;-26670h; SP:0_426;383_541;427_586;542_661;587_731; BP:0_93;-17_61;54_61;-16_61;14_84;      FC:i:64 RC:i:21"
            print("p"+node_id+"\t"+str(len(line.split()[2].strip(";").split(";"))))
            print("m"+node_id+"\t"+str(len(line.split()[2].strip(";").split(";"))))
    print(";")
    gfa_file.close()

def print_nodes_weight(gfa_file_name, DG=None):
    print("#id of nodes with their coverage. Here coverage refers to the read coverage of the less covered allele of all alleles of the fact")
    print("param ab :=")
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
    
    
def print_nodes_weight_phased_alleles(gfa_file_name, DG=None):
    print("#id of nodes with their coverage. Here coverage refers to the total number of reads that phased at least two alleles of the fact ")
    print("param ab_comapped :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            node_id=line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            "S       0       28175h;10031h;12786h;-41223l;-26670h; SP:0_426;383_541;427_586;542_661;587_731; BP:0_93;-17_61;54_61;-16_61;14_84;      FC:i:64 RC:i:21"
            print("p"+node_id+"\t"+line.split()[5].split(":")[-1])
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
    print("param comp :=")
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
            if type == "overlaps": 
                print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type)  
                # print the reverse of these nodes
                if sign_source == "m": sign_source="p"
                else: sign_source = "m"
                if sign_target == "m": sign_target="p"
                else: sign_target = "m"
                print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type)
                
                
                
            
            else: # All those nodes are non oriented - need all possible combinations
                for sign_source in "m","p":
                    for sign_target in "m", "p":
                        print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type)
                        if type=="successive":
                            print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type)
                
            
                
            
            # do not print other edges (unitig linked and snp linked)
    print(";")
    gfa_file.close()
    
def print_edge_coverages(gfa_file_name, DG=None):
    print ("#Coverage of the pairend links. Note that overlap links do not have any coverage (just computed from fact overlaps). Also incompatible and successive links do not have coverage, by definition. ")
    print ("param pairend_end_links_coverage :=")
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
            coverage=0
            if overlap_len==0: 
                type="links"
                coverage = int(line.split()[6].split(":")[-1])
            if overlap_len==-1: type="successive"
            if overlap_len==-2: type="incompatibles"
            # print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(coverage))
#             # do not print other edges (unitig linked and snp linked)
#             # print the reverse of these nodes
#             if sign_source == "m": sign_source="p"
#             else: sign_source = "m"
#             if sign_target == "m": sign_target="p"
#             else: sign_target = "m"
#             print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(coverage))
#
            if type == "overlaps": 
                print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(coverage))  
                # print the reverse of these nodes
                if sign_source == "m": sign_source="p"
                else: sign_source = "m"
                if sign_target == "m": sign_target="p"
                else: sign_target = "m"
                print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(coverage))
            
            else: # All those nodes are non oriented - need all possible combinations
                for sign_source in "m","p":
                    for sign_target in "m", "p":
                        print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(coverage))  
                        if type=="successive":
                            print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(coverage))  
            
            
            
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
            # print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(max(0,overlap_len)))
           #  # do not print other edges (unitig linked and snp linked)
           #
           #  # print the reverse of these nodes
           #  if sign_source == "m": sign_source="p"
           #  else: sign_source = "m"
           #  if sign_target == "m": sign_target="p"
           #  else: sign_target = "m"
           #  print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(max(0,overlap_len)))
           #           
            if type == "overlaps": 
                print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(max(0,overlap_len)))  
                # print the reverse of these nodes
                if sign_source == "m": sign_source="p"
                else: sign_source = "m"
                if sign_target == "m": sign_target="p"
                else: sign_target = "m"
                print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(max(0,overlap_len)))
            
            else: # All those nodes are non oriented - need all possible combinations
                for sign_source in "m","p":
                    for sign_target in "m", "p":
                        print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(max(0,overlap_len)))  
                        if type=="successive":
                            print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(max(0,overlap_len)))  
           
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
    print_nodes_number_loci(gfa_file_name, DG=None)
    print_nodes_weight(gfa_file_name,DG)
    print_nodes_weight_phased_alleles(gfa_file_name, DG)
    print_reverse(gfa_file_name,DG)
    if DG:
        print_nodes_connected_components(gfa_file_name, DG)
    print_edges(gfa_file_name,DG)
    print_edge_coverages(gfa_file_name, DG)
    print_edges_content(gfa_file_name,DG)



if __name__ == "__main__":
     main(sys.argv[1])
