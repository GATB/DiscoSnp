import sys
import getopt
#global usage: 
### first create connected components from disco (-A option)
#sh from_phased_alleles_to_clusters.sh phased_alleles_read_set_id_1.txt # creates file connected_components_phased_alleles_read_set_id_1.txt
### them from the .fa file, the id of the set your interested in (e.g. 1 for phased_alleles_read_set_id_1.txt, this will correspond to C1 coverage in the fa file), the file containing the connected components, and the phased_alleles_read_set_id_X.txt file, generate the fact file





def remove_non_variable_snps_from_coverages(coverages):
    coverages_bis={}
    for allele_id in coverages:
        if allele_id in coverages_bis: continue                     # SNP id already treated
        snp_id=allele_id[:-1]
        if coverages[snp_id+'l']>0 and coverages[snp_id+'h']>0:
            coverages_bis[snp_id+'l'] = coverages[snp_id+'l']
            coverages_bis[snp_id+'h'] = coverages[snp_id+'h']
    return coverages_bis
            

def store_abundances(coherent_fa_file,set_id, RemoveNonVariableSNPS, coverages=None):
    pos_coverage_determined=False
    pos_coverage=-1
    if not coverages:
        coverages={}                    #snp id (991h) -> coverage in the right read set 
    for oline in coherent_fa_file:  #>SNP_higher_path_991|P_1:30_C/G|high|nb_pol_1|C1_38|C2_0|Q1_0|Q2_0|G1_0/0:6,119,764|G2_1/1:664,104,6|rank_1
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
    if RemoveNonVariableSNPS: coverages=remove_non_variable_snps_from_coverages(coverages)
    return coverages
    

def store_cc(cc_file):
    cc={}
    for i,oline in enumerate (cc_file): # 852 1891 3484 2641 5758 3247
        oline=oline.rstrip().split()
        for idf in oline: 
            idf=int(idf)
            assert idf not in cc, "ERROR, idf is in more than one connected component"
            cc[idf]=i
    return cc
        


def store_phased_alleles(phased_alleles_file): ## ISSUE: Does not take right part of the pair information. # TODO
    phased_alleles={}
    for oline in phased_alleles_file:               #-129h_0;552l_38;-449h_33; => 2
        oline=oline.lstrip().rstrip()
        if oline[0]=='#': continue
        abundance = int(oline.split(' ')[-1])       # 2

        for pair_id in range(len(oline.split(' '))-2):          # Only one loop if data unpaired, two loops else
            ids = oline.split(' ')[pair_id].split(';')[:-1]     # -129h_0552l_38 -449h_33                
            # canonical representation: smallest first (removing with the strip function the eventual first '-' sign): 
            if int(ids[0].split('_')[0].strip('-')[:-1])>int(ids[-1].split('_')[0].strip('-')[:-1]):
                ids.reverse()
                # change the ortientation = change the sign: 
                for i in range(len(ids)):
                    if ids[i][0]=='-': ids[i]=ids[i][1:]
                    else: ids[i]='-'+ids[i]
                # change the orientation = change the distance to the previous one. 
                # -d_0;-c_5;-b_2;-a_4 is reversed into a_4;b_2_c_5_d_0 which is false. 
                # In fact a and b separated by 4, b and c by 2 and c and d by 5. 
                # d....c.b...a
                # Thus the final order is a_0;b_4;c_2;d_5
                previous="0"
                for i in range(len(ids)):
                    future_previous=ids[i].split('_')[-1]
                    ids[i]=ids[i].split('_')[0]+'_'+previous
                    previous=future_previous
            
            
            list_as_string = ""
            for aid in ids:
                list_as_string+=aid+';'
            # add the list to the phased_alleles or increase its count if not existing:
            if list_as_string in phased_alleles: 
                phased_alleles[list_as_string]+=abundance
            else: 
                phased_alleles[list_as_string]=abundance

                
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
    assert first_id in cc, "SNP"+str(first_id)+"in facts but not in connected components"
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
    
def remove_non_existing_or_non_variable_variants(phased_alleles,coverages):                                                             # ['-129h_0', '552l_38',  '-449h_33']
    """ Often a SNP is detected in one read set and used in a phased fact from an other read set in which it is not variable.
    In this case one of its higher or lower allele as a coverage 0.
    If this happens, we remove this SNP and we update the distance to previous one accordingly
    """
    distance_to_add_to_previous=0
    returned_list=[]
    for aid in phased_alleles:
        cud_allele = abs(int(aid.split('_')[0][:-1]))
        to_remove = False
        if str(cud_allele)+'h' not in coverages or coverages[str(cud_allele)+'h']==0: to_remove = True
        if str(cud_allele)+'l' not in coverages or coverages[str(cud_allele)+'l']==0: to_remove = True
        
        if not to_remove: 
            id_snp=aid.split('_')[0]
            distance_to_revious = int(aid.split('_')[-1])+distance_to_add_to_previous
            returned_list.append(id_snp+'_'+str(distance_to_revious))
            distance_to_add_to_previous=0
        else: 
            distance_to_add_to_previous+=int(aid.split('_')[-1])
            # print ("Warning, SNP "+aid+" does not exist. It was removed from list "+str(phased_alleles), file=sys.stderr)
    return returned_list
            
        
    

    
def print_djack_formated_phased_variants(coverages,cc,phased_alleles,RemoveNonVariableSNPS):
    # for aid in coverages:                                                                                               #snp id (991h) -> coverage
    #     current_snp_id=int(aid[:-1])                                                                                    #991
    #     if current_snp_id in cc:                                                                                        #necessary test?
    #         print("snp(cc"+str(cc[current_snp_id])+","+str(current_snp_id)+","+aid[-1]+","+str(coverages[aid])+").")    #"snp(cc_12,991h,coverage)"
            
    for i,list_as_string in enumerate(phased_alleles):                                                                  #'-129h_0;552l_38;-449h_33;': 2
        ids=list_as_string.split(';')[:-1]                                                                              # ['-129h_0', '552l_38',  '-449h_33']
        this_cc=check_phased_alleles_integrity_and_return_cc(ids,cc)
        if this_cc==-1: continue                                                                                        # There was a problem with this fact
        abundance = phased_alleles[list_as_string]                                                                      # 2

        if RemoveNonVariableSNPS: ids = remove_non_existing_or_non_variable_variants(ids, coverages)        
        for node_order,aid in enumerate(ids):                                                                           # ['-129h_0', '552l_38',  '-449h_33']
            # fact(cc_id, fact number, allele number in the fact, allele snp id, allele path (h/l), distance wrt to previous variant in the path)
            print("fact(cc"+str(this_cc)+","+str(i)+","+str(node_order+1)+","+aid.split('_')[0][:-1]+","+aid.split('_')[0][-1]+","+aid.split('_')[1]+").")
        if len(ids)>0:
            print("count("+str(i)+","+str(abundance)+").")
        
            
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
            coherent_fa_file = open(arg)
        elif opt in ("-u","--uncoherent_file"):
            use_uncoherent=True
            uncoherent_fa_file = open(arg)
        elif opt in ("-s","--set_id"):
            set_id = arg
        elif opt in ("-C","--connected_components_file"):
            cc_file = open(arg)
        elif opt in ("-p","--phased_alleles_file"):
            phased_alleles_file = open(arg)
        elif opt in ("-k","--keep_useless_SNPs"):
            RemoveNonVariableSNPS=False
        else:
            print("Unkwnown option {} ".format(opt))
            usage()
            sys.exit(2)
            
            
    # coherent_fa_file = open(sys.argv[1])
    # uncoherent_fa_file = open(sys.argv[2])
    # set_id = sys.argv[3]
    # cc_file = open(sys.argv[4])
    # phased_alleles_file = open(sys.argv[5])
    
    
    coverages=store_abundances(coherent_fa_file,set_id,RemoveNonVariableSNPS)                   # Store abundances from coherent snps
    if use_uncoherent:
        coverages=store_abundances(uncoherent_fa_file,set_id,RemoveNonVariableSNPS,coverages)   # Add abundances from uncoherent snps
    
    cc=store_cc(cc_file)
    phased_alleles=store_phased_alleles(phased_alleles_file)
    print_djack_formated_phased_variants(coverages,cc,phased_alleles,RemoveNonVariableSNPS)
    
if __name__ == "__main__":
    main()
    