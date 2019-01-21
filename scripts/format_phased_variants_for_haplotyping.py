import sys
import getopt
#global usage: 
### first create connected components from disco (-A option)
#sh from_phased_alleles_to_clusters.sh phased_alleles_read_set_id_1.txt # creates file connected_components_phased_alleles_read_set_id_1.txt
### them from the .fa file, the id of the set your interested in (e.g. 1 for phased_alleles_read_set_id_1.txt, this will correspond to C1 coverage in the fa file), the file containing the connected components, and the phased_alleles_read_set_id_X.txt file, generate the fact file





def remove_non_variable_snps_from_set(my_set):
    return my_set # Useless now as input is non ghost in the pipeline
    my_set_bis={}
    for allele_id in my_set:
        if allele_id in my_set_bis: continue                     # SNP id already treated
        snp_id=allele_id[:-1]
        if my_set[snp_id+'l']>0 and my_set[snp_id+'h']>0:
            my_set_bis[snp_id+'l'] = my_set[snp_id+'l']
            my_set_bis[snp_id+'h'] = my_set[snp_id+'h']
    return my_set_bis
            

def store_abundances(fa_file_name,set_id, RemoveNonVariableSNPS, coverages=None):
    fa_file = open(fa_file_name)
    pos_coverage_determined=False
    pos_coverage=-1
    if not coverages:
        coverages={}                    #snp id (991h) -> coverage in the right read set 
    for oline in fa_file:  #>SNP_higher_path_991|P_1:30_C/G|high|nb_pol_1|C1_38|C2_0|Q1_0|Q2_0|G1_0/0:6,119,764|G2_1/1:664,104,6|rank_1
        if oline[0] != '>': continue
        line=oline.rstrip().split('|')
        id=line[0].split('_')[-1]    #here 991
        id+=line[0].split('_')[1][0] #'h' or 'l'
        if not pos_coverage_determined:
            for pos_coverage in range(len(line)):
                if line[pos_coverage][0]=='C':
                    value=line[pos_coverage][1:].split('_')[0]
                    if value==set_id:
                        pos_coverage_determined=True
                        break
            #if not pos_coverage_determined:
            assert pos_coverage_determined, "Set id "+ str(set_id)+ " not findable in header like "+ oline.rstrip()
        coverages[id]=int(line[pos_coverage].split('_')[1]) # get the right coverage corresponding to the searche read set
    if RemoveNonVariableSNPS: coverages=remove_non_variable_snps_from_set(coverages)
    fa_file.close()
    return coverages
    

def get_upper_sequence_size(sequence):
    upper_sequence_size=0
    for l in sequence:
        if l>='A' and l<='Z':
            upper_sequence_size+=1
    return upper_sequence_size

def store_variant_sizes(fa_file_name,set_id, RemoveNonVariableSNPS, sizes=None):
    fa_file = open(fa_file_name)
    if not sizes:
        sizes={}                    #snp id (991h) -> coverage in the right read set 
    while True:
        header = fa_file.readline()
        if not header: break
        header=header.rstrip().split('|')
        sequence = fa_file.readline()
        id=header[0].split('_')[-1]    #here 991
        id+=header[0].split('_')[1][0] #'h' or 'l'
        sizes[id]=get_upper_sequence_size(sequence)
    if RemoveNonVariableSNPS: sizes=remove_non_variable_snps_from_set(sizes)
    fa_file.close()
    return sizes

    

def store_cc(cc_file):
    cc={}
    for i,oline in enumerate (cc_file): # 852 1891 3484 2641 5758 3247
        oline=oline.rstrip().split()
        for idf in oline: 
            idf=int(idf)
            assert idf not in cc, "ERROR, idf is in more than one connected component"
            cc[idf]=i
    return cc
        


def store_phased_alleles(phased_alleles_file_name): 
    phased_alleles_file = open(phased_alleles_file_name)
    phased_alleles={}
    for oline in phased_alleles_file:               #-129h_0;552l_38;-449h_33; => 2
        oline=oline.lstrip().rstrip()
        if oline[0]=='#': continue
        abundance = int(oline.split(' ')[-1])       # 2
        
        for pair_id in range(len(oline.split(' '))-2):          # Only one loop if data unpaired, two loops else
            ids = oline.split(' ')[pair_id].split(';')[:-1]     # -129h_0 0552l_38 -449h_33                
            # canonical representation: smallest first (removing with the strip function the eventual first '-' sign):
            if int(ids[0].split('_')[0].strip('-')[:-1])>int(ids[-1].split('_')[0].strip('-')[:-1]):
                ids.reverse()
                # change the ortientation = change the sign:
                for i in range(len(ids)):
                    if ids[i][0]=='-': ids[i]=ids[i][1:]
                    else: ids[i]='-'+ids[i]

            
            list_as_string = ""
            for aid in ids:
                list_as_string+=aid.split('_')[0]+';' # concatenation of ids 
            # add the list to the phased_alleles or increase its count if not existing:
            if list_as_string in phased_alleles: 
                phased_alleles[list_as_string]+=abundance
            else: 
                phased_alleles[list_as_string]=abundance

    phased_alleles_file.close()
    return phased_alleles

# def check_phased_alleles_integrity_and_return_cc_only_assert(phased_alleles):
#     snp_ids=set()
#     first_id=abs(int(phased_alleles[0].split('_')[0][:-1]))                                                             # 129
#     snp_ids.add(first_id)
#     assert first_id in cc, "SNP"+str(first_id)+"in facts but not in connected components"
#     this_cc=cc[first_id]                                                                                                # eg 555
#     for j in range(1,len(phased_alleles)):                                                                              # all other allele ids
#         cur_id = abs(int(phased_alleles[j].split('_')[0][:-1]))
#         assert cur_id in cc, "SNP"+str(first_id)+"in facts but not in connected components"
#         assert cc[cur_id] == this_cc, "impossible all variants from "+list_as_string+ "are not in the same CC"
#         assert cur_id not in snp_ids, "impossible, "+str(cur_id)+" exists several times in phased variants "+str(phased_alleles)
#         snp_ids.add(cur_id)
#     return this_cc

def check_phased_alleles_integrity_and_return_cc(phased_alleles,cc):                                                    # ['-129h_0', '552l_38',  '-449h_33']
    snp_ids=set()
    first_id=abs(int(phased_alleles[0].split('_')[0][:-1]))                                                             # 129
    snp_ids.add(first_id)
    assert first_id in cc, "SNP "+str(first_id)+" in facts but not in connected components"
    this_cc=cc[first_id]                                                                                                # eg 555
    for j in range(1,len(phased_alleles)):                                                                              # all other allele ids
        cur_id = abs(int(phased_alleles[j].split('_')[0][:-1]))
        assert cur_id in cc, "SNP"+str(first_id)+"in facts but not in connected components"
        assert cc[cur_id] == this_cc, "impossible all variants from "+list_as_string+ "are not in the same CC"
        if cur_id in snp_ids:
            print ("Warning, "+str(cur_id)+" exists several times in phased variants "+str(phased_alleles)+". We skip this phased variants list", file=sys.stderr)
            return -1
        snp_ids.add(cur_id)
    return this_cc
    
def remove_non_existing_or_non_variable_variants(phased_alleles,coverages):                                                             # ['-129h', '552l',  '-449h']
    """ Often a SNP is detected in one read set and used in a phased fact from an other read set in which it is not variable.
    In this case one of its higher or lower allele as a coverage 0.
    If this happens, we remove this SNP and we update the distance to previous one accordingly
    """
    distance_to_add_to_previous=0
    returned_list=[]
    first=True
    for aid in phased_alleles:
        cud_allele = abs(int(aid[:-1]))
        to_remove = False
        if str(cud_allele)+'h' not in coverages or coverages[str(cud_allele)+'h']==0: to_remove = True
        if str(cud_allele)+'l' not in coverages or coverages[str(cud_allele)+'l']==0: to_remove = True
        
        if not to_remove: 
            returned_list.append(aid)
    return returned_list
            
        
    
def print_distances(phased_alleles_file_name,sizes):
        distances={}                                    #-129 -> 552 -> 38 (for each value: a list of right hand pair with its distance)
        phased_alleles_file = open(phased_alleles_file_name)
        phased_alleles={}
        for oline in phased_alleles_file:               #-129h_0;552l_38;-449h_33; => 2
            oline=oline.lstrip().rstrip()
            if oline[0]=='#': continue
            for pair_id in range(len(oline.split(' '))-2):          # Only one loop if data unpaired, two loops else
                ids=oline.split(' ')[pair_id].split(';')[:-1]       # -129h_0 0552l_38 -449h_33         
                for phasing_pair_id in range(len(ids)-1):
                    
                    variant_R1=ids[phasing_pair_id]                 # -129h_0
                    variant_R2=ids[phasing_pair_id+1]               # 0552l_38
                    id_R1=variant_R1.split('_')[0]                  # -129h
                    id_R2=variant_R2.split('_')[0]                  # 0552l
                    
                    abs_id_R1_allele=variant_R1.split('_')[0]       # 129h
                    abs_id_R2_allele=variant_R2.split('_')[0]       # 0552l
                    if(id_R1[0]=='-'): abs_id_R1_allele=id_R1[1:]
                    else: abs_id_R1_allele=id_R1
                    if(id_R2[0]=='-'): abs_id_R2_allele=id_R2[1:]
                    else: abs_id_R2_allele=id_R2
                    
                    pos_R1=variant_R1.split('_')[1]                 # 0
                    pos_R2=variant_R2.split('_')[1]                 # 38
                    dist_R1_to_R2=int(pos_R2)                       # 38
                        #---------R1-----------<---------x------->
                        #<-l->----------------R2------------------
                        #x=l+|R2|-|R1|
                        # here x is dist_R2_to_R1
                        # here l is dist_R1_to_R2
                    dist_R2_to_R1=dist_R1_to_R2+int(sizes[abs_id_R2_allele]-sizes[abs_id_R1_allele])
                    # store values
                    if id_R1 not in distances: distances[id_R1]=set()
                    distances[id_R1].add((id_R2,dist_R1_to_R2))
                    if id_R1[0]=='-': rc_id_R1=id_R1[1:]
                    else: rc_id_R1='-'+id_R1
                    if id_R2[0]=='-': rc_id_R2=id_R2[1:]
                    else: rc_id_R2='-'+id_R2
                    if rc_id_R2 not in distances: distances[rc_id_R2]=set()
                    distances[rc_id_R2].add((rc_id_R1,dist_R2_to_R1))
                # print (ids)
                # print (distances)

        for id_first in distances: # {'-1003': {('-492', 38)}, '492': {('1003', 38)}}
            if id_first[0]=='-' :
                direction_first   ='n'
                abs_id_first=id_first[1:]
            else:
                direction_first   ='p'
                abs_id_first=id_first
            for values in distances[id_first]:
                id_second=values[0]
                if id_second[0]=='-' :
                    direction_second   ='n'
                    abs_id_second=id_second[1:]
                else:
                    direction_second   ='p'
                    abs_id_second=id_second
                print("dist("+abs_id_first[:-1]+","+direction_first+","+abs_id_first[-1]+","+abs_id_second[:-1]+","+direction_second+","+abs_id_second[-1]+","+str(values[1])+").")
        # sys.exit(0)
        
        
    
def print_formated_phased_variants(coverages,cc,phased_alleles,RemoveNonVariableSNPS):
    """
    Prints paths : 
    a fact is composed of all pairwise links (fact) (one per line) and a count
    """          
    for i,list_as_string in enumerate(phased_alleles):                                                                  #'-129h_0;552l_38;-449h_33;': 2
        ids=list_as_string.split(';')[:-1]                                                                              # ['-129h_0', '552l_38',  '-449h_33']
        this_cc=check_phased_alleles_integrity_and_return_cc(ids,cc)
        if this_cc==-1: continue                                                                                        # There was a problem with this fact
        abundance = phased_alleles[list_as_string]                                                                      # 2

        if RemoveNonVariableSNPS: ids = remove_non_existing_or_non_variable_variants(ids, coverages)        
        for node_order,aid in enumerate(ids):                                                                           # ['-129h', '552l',  '-449h']
            path_id     =aid[:-1]
            # direction (p or m)
            direction   ='p'
            if path_id[0]=='-':
                direction='n'
                path_id=path_id[1:]
            path_hl     =aid[-1]
            # fact(cc_id, fact number, allele number in the fact, allele snp id, allele direction (p/n), allele path (h/l), distance wrt to previous variant in the path)
            print("fact(cc"+str(this_cc)+","+str(i)+","+str(node_order+1)+","+path_id+","+direction+","+path_hl+").")#+","+distance_to_previous+").")
        if len(ids)>0:
            print("count("+str(i)+","+str(abundance)+").")
            
            
def print_variants(coverages,cc):
    """
    Prints all variants that are in connected components 
    """
    for aid in coverages:                                                                                               #snp id (991h) -> coverage
        current_snp_id=int(aid[:-1])                                                                                    #991
        if current_snp_id in cc:                                                                                        #necessary test?
            print("snp(cc"+str(cc[current_snp_id])+","+str(current_snp_id)+","+aid[-1]+","+str(coverages[aid])+").")    #"snp(cc_12,991h,coverage)"
            
def usage():
    usage= """
    #########################################
    format_phased_variants_for_haplotyping.py
    #########################################
    
    -h --help : print this message
    --coherent_file:            <file>.fa:  coherent fa file from discoSnp    [Mandatory]
    --uncoherent_file:          <file>.fa:  uncoherent fa file from discoSnp  [Optional]
    --set_id:                   int:        id of the read set for which one has the phasing information. If set to 'i': corresponds to Ci in the coverage .fa file(s) [Mandatory]
    --connected_components_file <file>.txt: file containing connected components from SNPs. Computed with "sh from_phased_alleles_to_clusters.sh phased_alleles_read_set_id_i.txt"
    --phased_alleles_file       <file>.txt: file from discoSnp (-A option) named phased_alleles_read_set_id_i.txt with 'i' being set_id
    --keep_useless_SNPs         Boolean:    If set: conserve SNPs that are not used in the set_id i (only one of the allele is used). [Optional, default: False]
    
    ----Classical usage example----
    python3 format_phased_variants_for_haplotyping.py --coherent_file disco_specrep_k_31_c_3_D_0_P_10_b_2_coherent.fa --uncoherent_file disco_specrep_k_31_c_3_D_0_P_10_b_2_uncoherent.fa --set_id 2 --connected_components_file connected_components_phased_alleles_read_set_id_2.txt --phased_alleles_file phased_alleles_read_set_id_2.txt  > facts_disco_specrep_k_31_c_3_D_0_P_10_b_2.txt
    """
    print(usage)
    
###OPTIONS 
def main():
    use_uncoherent=False
    RemoveNonVariableSNPS=True      ### With this option se't to True, SNPs having a coverage h or l equal to 0 are removed and they are removed from facts. 

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:u:s:C:p:k",["help","coherent_file=","uncoherent_file=","set_id=","connected_components_file=","phased_alleles_file=","keep_useless_SNPs"])
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
            cc_file = open(arg)
        elif opt in ("-p","--phased_alleles_file"):
            phased_alleles_file_name = arg
        elif opt in ("-k","--keep_useless_SNPs"):
            RemoveNonVariableSNPS=False
        else:
            print("Unkwnown option {} ".format(opt))
            usage()
            sys.exit(2)
            
            
    
    
    coverages=store_abundances(coherent_fa_file_name,set_id,RemoveNonVariableSNPS)                      # Store abundances from coherent snps
    sizes=store_variant_sizes(coherent_fa_file_name,set_id,RemoveNonVariableSNPS)                       # Store sizes from coherent snps
    if use_uncoherent:
        coverages=store_abundances(uncoherent_fa_file_name,set_id,RemoveNonVariableSNPS,coverages)      # Add abundances from uncoherent snps
        sizes=store_variant_sizes(uncoherent_fa_file_name,set_id,RemoveNonVariableSNPS,sizes)           # Add sizes from uncoherent snps
    
    cc=store_cc(cc_file)
    phased_alleles=store_phased_alleles(phased_alleles_file_name)
    print_variants(coverages,cc)
    print_distances(phased_alleles_file_name,sizes)
    print_formated_phased_variants(coverages,cc,phased_alleles,RemoveNonVariableSNPS)
    
if __name__ == "__main__":
    main()
    