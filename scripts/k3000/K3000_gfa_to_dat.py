import sys
import K3000_gfa_post_treatment as gpt # :)
# from deprecated import deprecated


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
p0	m55	overlaps 
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
    print("# *** note about abundances ***")
    print("#  This file contains four values of abundances per fact")
    print("#   - ab_min: the coverage of the less covered variant of the fact")
    print("#   - ab_max: the coverage of the most covered variant of the fact")
    print("#   - ab_avg: the average coverage of the variants of the fact")
    print("#   - ab_comapped: total number of reads i/ mapped on the fact and ii/ mapping at least two variants of the fact")
    print("#  Full example")
    print("#  *  *  *  *  *  *      (variant positions on the fact)")
    print("#  --------              (read mapping at least two variants)")
    print("#    --------            (read mapping at least two variants)")
    print("#    -----               (read mapping at least two variants)")
    print("#       --------         (read mapping at least two variants)")
    print("#             ----       (read mapping at least two variants)")
    print("#          --------      (read mapping at least two variants)")
    print("# --                     (read mapping one variant)")
    print("# --                     (read mapping one variant)")
    print("#    --                  (read mapping one variant)")
    print("#       --               (read mapping one variant)")
    print("#       --               (read mapping one variant)")
    print("#          --            (read mapping one variant)")
    print("#             --         (read mapping one variant)")
    print("#                --      (read mapping one variant)")
    print("#  3  4  6  4  4  2  (for each variant how many reads mapped on it)")
    print("#  Here: ab_comapped=6, ab_min=2, ab_max=6, ab_avg=3.84")
    print("# *** End note about abundances ***")
    
    print("param mult := XXX")
    
oposite_sign = lambda x: "m" if x == "p" else "p"

def print_nodes(gfa_file_name, DG=None):
    print("#id of plus (p) and minus (m) nodes")
    print("set V :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        if not line: break
        line=line.strip()
        if line[0]=="S":
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
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
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
            print("p"+node_id+"\t"+str(len(line.split()[2].strip(";").split(";"))))
            print("m"+node_id+"\t"+str(len(line.split()[2].strip(";").split(";"))))
    print(";")
    gfa_file.close()

def print_nodes_ab_min(gfa_file_name, DG=None):
    print("#id of nodes with their ab_min coverage")
    print("param ab_min :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            node_id=line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
            print("p"+node_id+"\t"+line.split()[6].split(":")[-1])
            assert(line.split()[6].split(":")[0]=="min")
    print(";")
    gfa_file.close()
    
def print_nodes_ab_max(gfa_file_name, DG=None):
    print("#id of nodes with their ab_max coverage")
    print("param ab_max :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            node_id=line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
            print("p"+node_id+"\t"+line.split()[7].split(":")[-1])
            assert(line.split()[7].split(":")[0]=="max")
    print(";")
    gfa_file.close()
    
def print_nodes_ab_avg(gfa_file_name, DG=None):
    print("#id of nodes with their ab_avg coverage")
    print("param ab_avg :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            node_id=line.split()[1]
            if DG and node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
            print("p"+node_id+"\t"+line.split()[8].split(":")[-1])
            assert(line.split()[8].split(":")[0]=="mean")
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
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
            print("p"+node_id+"\t"+line.split()[5].split(":")[-1])
            assert(line.split()[5].split(":")[0]=="FC")
    print(";")
    gfa_file.close()
    

def print_reverse(gfa_file_name,DG=None):
    print("#for each forward node, indicates the id of the reverse version")
    print("set reverse :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
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
            "S       0       1234h;-1354h;977h;577h; SP:0_54;53_96;64_117;84_128;    BP:0_41;6_41;-25_41;-26_41;     FC:i:579        min:410 max:467 mean:441.5      AC:436;467;453;410;"
            node_id=line.split()[1]
            if node_id not in DG.nodes(): continue       # this node was removed during gfa post treatment
            cc_id = DG.nodes[node_id]['cc_id']
            print("p"+node_id+"\t"+str(cc_id))
            print("m"+node_id+"\t"+str(cc_id))
    print(";")
    gfa_file.close()
    
    
# @deprecated(reason="Not used anymore, redundant with function `already_used` generic for all kind of links")
# def printable_successive(printed_successive, source_id, target_id):
#
#     """
#     1/ checks of source_id <-> target_id not already in printed_successive. If already inside, do nothing and returns false
#     2/ [else] add target_id in source_id key and return true
#     """
#     # source is always the smallest, avoid to test both directions.
#     if target_id < source_id:
#         tmp = source_id
#         source_id = target_id
#         target_id = tmp
#     if source_id in printed_successive and target_id in printed_successive[source_id]: return False
#
#     if source_id not in printed_successive: printed_successive[source_id] = []
#     printed_successive[source_id].append(target_id)
#     return True
    

def already_used(edges,sign_source,source_id,sign_target,target_id,type):
    """ detects if the edge had already been used 
    This is necessary to avoid duplicated edges. 
    For instance if the gfa contains: 
    L	571	+	829	-	-1M
    ...
    L	829	+	571	-	-1M
    Then, the line is duplicated as it contains the same piece of information. 
    Hence we have to store the fact that the edge was already printed with: L	571	+	829	-	-1M for instance. 
    We do this by storing for the lowest node id.
    """
    if type not in edges: edges[type] = {}
    if source_id > target_id:
        if sign_source == "m": sign_source="p"
        else: sign_source = "m"
        if sign_target == "m": sign_target="p"
        else: sign_target = "m"
        tmp = source_id
        source_id = target_id
        target_id = tmp
        
    source = sign_source+source_id
    target = sign_target+target_id
    already = True
    if source not in edges[type]:
        edges[type][source] = set()
    if target not in edges[type][source]:
        edges[type][source].add(target)
        already = False
    
    if type == "successive":
        sign_source = oposite_sign(sign_source)
        sign_target = oposite_sign(sign_target)
        source = sign_source+source_id
        target = sign_target+target_id
        if source not in edges[type]:
            edges[type][source] = set()
        if target not in edges[type][source]:
            edges[type][source].add(target)
            already = False
            
    return already




def print_edges(gfa_file_name, DG=None):
    print("#set of edges. Four types of edges, 1/ \"overlaps\" edges, that show an overlap between facts and 2/ \"links\" edges, that represent facts linked by paired reads (distanace unknown) and 3/ \"successive\" edges that represent two successive facts (without phasing) and 4/ \"incompatible\" edges, no path should contain two nodes linked by such an edge ")
    printed_edges = {}
    print("set Edges :=")
    gfa_file = open(gfa_file_name)
    # printed_successive = {} # key = id, target = [ids]. Used to retain which successive links ad been writen, avoiding redundancies
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
            if already_used(printed_edges,sign_source,source_id,sign_target,target_id,type): continue
            if type == "overlaps" or type=="successive" or type=="links": 
                print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type)  
                # print the reverse of these nodes
                sign_source = oposite_sign(sign_source)
                sign_target = oposite_sign(sign_target)
                print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type)
            else: # All those nodes are non oriented - need all possible combinations
                for sign_source in "m","p":
                    for sign_target in "m", "p":
                        print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type)
    print(";")
    gfa_file.close()
    
def print_edge_coverages(gfa_file_name, DG=None):
    print ("#Coverage of the pairend links. Note that overlap links do not have any coverage (just computed from fact overlaps). Also incompatible and successive links do not have coverage, by definition. ")
    print ("param pairend_end_links_coverage :=")
    printed_edges = {}
    gfa_file = open(gfa_file_name)
    # printed_successive = {} # key = id, target = [ids]. Used to retain which successive links ad been writen, avoiding redundancies
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

            if already_used(printed_edges,sign_source,source_id,sign_target,target_id,type): continue
            if type == "overlaps"  or type=="successive" or type=="links": 
                print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(coverage))  
                # print the reverse of these nodes
                sign_source = oposite_sign(sign_source)
                sign_target = oposite_sign(sign_target)
                print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(coverage))
            
            else: # All those nodes are non oriented - need all possible combinations
                for sign_source in "m","p":
                    for sign_target in "m", "p":
                        print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(coverage))  
    print(";")
    gfa_file.close()
    
def print_edges_content(gfa_file_name, DG=None):
    print("#overlap length of each edge. For an \"overlaps\" edge, it indicates the number of common variants. For any other edge type (links, successive, or incompatibles), this is set to zero")
    print("param l :=")
    printed_edges = {}
    gfa_file = open(gfa_file_name)
    # printed_successive = {} # key = id, target = [ids]. Used to retain which successive links ad been writen, avoiding redundancies
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
            if already_used(printed_edges,sign_source,source_id,sign_target,target_id,type): continue
            if type == "overlaps" or type=="successive" or type=="links": 
                print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(max(0,overlap_len)))  
                # print the reverse of these nodes
                sign_source = oposite_sign(sign_source)
                sign_target = oposite_sign(sign_target)
                print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(max(0,overlap_len)))
            
            else: # All those nodes are non oriented - need all possible combinations
                # if type=="successive" and  not printable_successive(printed_successive, source_id, target_id): continue
                # if type=="links":
                #     sign_source = "p"
                #     sign_target = "p"
                #     print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(max(0,overlap_len)))
                #     continue
                for sign_source in "m","p":
                    for sign_target in "m", "p":
                        print(sign_target+target_id+"\t"+sign_source+source_id+"\t"+type+"\t"+str(max(0,overlap_len)))  
                        # if type=="successive":              # in case of successive edges, one needs to indicate all possibilities p/m and forward and reverse
                        #     print(sign_source+source_id+"\t"+sign_target+target_id+"\t"+type+"\t"+str(max(0,overlap_len)))
           
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

    sys.stderr.write("Storing the graph\n")
    DG = gpt.store_graph(gfa_file_name) # store the graph without paired edges

    sys.stderr.write("Detecting connected components\n")
    gpt.assign_cc(DG)

    sys.stderr.write("Removing connected compoenents with cycles\n")
    gpt.remove_cc_with_cycles(DG)
    
    sys.stderr.write("Printing out the dat file\n")
    print_header()
    print_nodes(gfa_file_name,DG)
    print_nodes_number_loci(gfa_file_name, DG)
    print_nodes_ab_min(gfa_file_name,DG)
    print_nodes_ab_max(gfa_file_name,DG)
    print_nodes_ab_avg(gfa_file_name,DG)
    print_nodes_weight_phased_alleles(gfa_file_name, DG)
    print_reverse(gfa_file_name,DG)
    if DG:
        print_nodes_connected_components(gfa_file_name, DG)
    print_edges(gfa_file_name,DG)
    print_edge_coverages(gfa_file_name, DG)
    print_edges_content(gfa_file_name,DG)



if __name__ == "__main__":
     main(sys.argv[1])
