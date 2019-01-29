import sys
import getopt
import linecache

def store_connected_component_nodes(connected_components_file_name, connected_component_id): # connected_component_id is the line corresponding to the cc of interest (one-based)
    connected_componant_nodes=set()
    cc=linecache.getline(connected_components_file_name, connected_component_id).strip().split(" ")
    for value in cc:
        connected_componant_nodes.add(value)
    return connected_componant_nodes
    

def store_abundances(fa_file_name,set_id,connected_component_nodes, node_coverages=None):
    fa_file = open(fa_file_name)
    pos_coverage_determined=False
    pos_coverage=-1
    if not node_coverages:
        node_coverages={}                    #snp id (991h) -> coverage in the right read set 
    for oline in fa_file:               #>SNP_higher_path_991|P_1:30_C/G|high|nb_pol_1|C1_38|C2_0|Q1_0|Q2_0|G1_0/0:6,119,764|G2_1/1:664,104,6|rank_1
        if oline[0] != '>': continue
        line=oline.rstrip().split('|')
        variant_id=line[0].split('_')[-1]    #here 991
        if variant_id not in connected_component_nodes: continue
        variant_id+=line[0].split('_')[1][0] #'h' or 'l'
        if not pos_coverage_determined:
            for pos_coverage in range(len(line)): 
                if line[pos_coverage][0]=='C':
                    value=line[pos_coverage][1:].split('_')[0]
                    if value==set_id:
                        pos_coverage_determined=True
                        break
            #if not pos_coverage_determined:
            assert pos_coverage_determined, "Set id "+ str(set_id)+ " not findable in header like "+ oline.rstrip()
        node_coverages[variant_id]=int(line[pos_coverage].split('_')[1]) # get the right coverage corresponding to the searche read set
    fa_file.close()
    return node_coverages
    
    
def store_existing_nodes(phased_alleles_file_name,connected_component_nodes):
    phased_alleles_file = open(phased_alleles_file_name)
    existing_ids = set()
    for line in phased_alleles_file:
        line=line.strip()
        if line[0]=='#': continue       # ./fof_unpaired.txt : comment
        line = line.split(' ')          # -1000183h_0;3224567l_287; -10183h_0;267l_287; => 1
        for i in range(len(line)-2):     # either only one pass or two passed if the line contains pairend data
            all_infos = line[i][:-1].split(';')
            for all_info in all_infos:
                if all_info[0]=='-':
                    all_info=all_info[1:]
                all_info = all_info.split('_')[0]
                if all_info[:-1] in connected_component_nodes:
                    existing_ids.add(all_info)    
    phased_alleles_file.close()
    return existing_ids

def print_existing_nodes(existing_ids):
    # print ("param n := "+str(len(existing_ids))+";")
    print ("set V:=")
    for node_id in existing_ids:
        print (node_id)
    

def print_existing_nodes_and_coverage(existing_ids, node_coverages):
    print ("param s:=")
    for node_id in existing_ids:
        print (node_id, node_coverages[node_id])
    
def line_in_connected_component(line, connected_component_nodes):
    # -1000183h_0;3224567l_287; -10183h_0;267l_287; => 1
    for i in range(len(line)-2):     # either only one pass or two passed if the line contains pairend data
        all_infos = line[i][:-1].split(';')
        for all_info in all_infos:
            current = all_info.split("_")[0][:-1]
            if current[0]=='-': current=current[1:]
            if current in connected_component_nodes: return True
    return False
    
    
def transform_to_mp(id_var):

    if id_var[0]=='-':  return 'm'+id_var[1:]
    else:               return 'p'+id_var

def print_existing_edges(phased_alleles_file_name,connected_component_nodes):
    print ("set E:=")
    phased_alleles_file = open(phased_alleles_file_name)
    edges={}
    for line in phased_alleles_file:
        line=line.strip()
        if line[0]=='#': continue       # ./fof_unpaired.txt : comment
        line = line.split(' ')          # -1000183h_0;3224567l_287; -10183h_0;267l_287; => 1
        if not line_in_connected_component(line, connected_component_nodes): continue
        path_abundance = line[-1]
        for i in range(len(line)-2):     # either only one pass or two passed if the line contains pairend data
            all_infos = line[i][:-1].split(';')
            previous=all_infos[0].split("_")[0]
            for all_info in all_infos[1:]:
                current = all_info.split("_")[0]
                if (current)<previous: 
                    key=current
                    value=previous
                else: 
                    key=previous
                    value=current
                    if key[0]=='-':     key=key[1:]
                    else:               key='-'+key
                    if value[0]=='-':   value=value[1:]
                    else:               value='-'+value
                    
                if key not in edges: 
                    edges[key]=set()
                edges[key].add(value)
                previous=current
    phased_alleles_file.close()
    
    for key in edges:
        for value in edges[key]:
            print (transform_to_mp(key)+" "+transform_to_mp(value))
    

def print_paths(phased_alleles_file_name, connected_component_nodes, print_fact_id=False):
    phased_alleles_file = open(phased_alleles_file_name)
    # NOTE: this function does not print the full path when pairend data are used. In this case it output two distinct paths and looses the link.
    fact_id=1
    printed_fact_id=""
    if print_fact_id:   print ("param l:=")
    else:               print ("set facts:=")
    for line in phased_alleles_file:
        line=line.strip()
        if line[0]=='#': continue       # ./fof_unpaired.txt : comment
        line = line.split(' ')          # -1000183h_0;3224567l_287; -10183h_0;267l_287; => 1
        if not line_in_connected_component(line, connected_component_nodes): continue
        path_abundance = line[-1]
        for i in range(len(line)-2):     # either only one pass or two passes if the line contains pairend data
            all_infos = line[i][:-1].split(';')

            if print_fact_id: 
                printed_fact_id=str(fact_id)+" "
                fact_id+=1
            for i in range(len(all_infos)-1):
                print(printed_fact_id+transform_to_mp(all_infos[i].split("_")[0])+" "+transform_to_mp(all_infos[i+1].split("_")[0])+" "+path_abundance)
                
        
    if print_fact_id: print(";\nparam nb_facts := "+str(fact_id-1))
    phased_alleles_file.close()
    

        
        
def usage():
    usage= """
    #########################################
    from_passed_nodes_to_dot_for_flot.py
    #########################################
    
    -h --help : print this message
    ** Mandatory ** 
    --coherent_file:            <file>.fa:  coherent fa file from discoSnp    
    --uncoherent_file:          <file>.fa:  uncoherent fa file from discoSnp 
    --set_id:                   int:        id of the read set for which one has the phasing information. If set to 'i': corresponds to Ci in the coverage .fa file(s) [Mandatory]
    --connected_components_file <file>.txt: file containing connected components from SNPs. Computed with "sh from_phased_alleles_to_clusters.sh phased_alleles_read_set_id_i.txt"
    --connected_component_id    int:        id (one-based) of the connected component to consider, 1=first line from the 'connected_components_file', and so on. 
    --phased_alleles_file       <file>.txt: file from discoSnp (-A option) named phased_alleles_read_set_id_i.txt with 'i' being set_id

    
    ----Classical usage example----
     python3 from_passed_nodes_to_dot_for_flot.py --uncoherent_file disco_specrep_k_31_c_3_D_0_P_10_b_2_uncoherent.fa --coherent_file disco_specrep_k_31_c_3_D_0_P_10_b_2_coherent.fa --set_id 2 --phased_alleles_file phased_alleles_read_set_id_2.txt  --connected_components_file connected_components_phased_alleles_read_set_id_2.txt   --connected_component_id $i
    """
    print(usage)
    
###OPTIONS 
def main():
    use_uncoherent=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:u:C:i:s:p:",["help","coherent_file=","uncoherent_file=","connected_components_file=","connected_component_id=","set_id=","phased_alleles_file=","keep_useless_SNPs"])
        if not opts:
            usage()
            sys.exit(2)
    except getopt.GetoptError as e:
        print(e)
        usage()
        sys.exit(2)
    for opt, arg in opts : 
        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-c","--coherent_file"):
            coherent_fa_file_name = arg
        elif opt in ("-u","--uncoherent_file"):
            use_uncoherent=True
            uncoherent_fa_file_name = arg
        elif opt in ("-s","--set_id"):
            set_id = arg
        elif opt in ("-C","--connected_components_file"):
            connected_components_file_name = arg
        elif opt in ("i","--connected_component_id"):
            connected_component_id = int(arg)
        elif opt in ("-p","--phased_alleles_file"):
            phased_alleles_file_name = arg
        else:
            print("Unkwnown option {} ".format(opt))
            usage()
            sys.exit(2)
            
            
    
    
    connected_component_nodes = store_connected_component_nodes(connected_components_file_name, connected_component_id)
    existing_ids=store_existing_nodes(phased_alleles_file_name,connected_component_nodes)
    print_existing_nodes(existing_ids)
    print (';')
    print_existing_edges(phased_alleles_file_name,connected_component_nodes)
    print(';')
    node_coverages=store_abundances(coherent_fa_file_name,set_id,connected_component_nodes)                      # Store abundances from coherent snps
    node_coverages=store_abundances(uncoherent_fa_file_name,set_id,connected_component_nodes,node_coverages)      # Add abundances from uncoherent snps

    
    print_paths(phased_alleles_file_name,connected_component_nodes)
    print(';')

    print_existing_nodes_and_coverage(existing_ids, node_coverages)
    print(';')
        
    print_paths(phased_alleles_file_name,connected_component_nodes, True)
    print(';')
    
if __name__ == "__main__":
    main()
        