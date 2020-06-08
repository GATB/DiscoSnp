import sys

def get_left_clean_snp(snp):
    return snp.lstrip().lstrip('-')


def set_indexes_from_gfa(gfa_file_name):
    mfile = open(gfa_file_name)
    snp_to_fact_id={}   # each SNP id (without order nor h or l) is linked to some compacted facts id
                        # eg: '1015': {8, 5, 6, 7}
    facts={}            # each compacted fact is accessible by a unique id. We store only the set of SNPs ids (not their orientering or their intra distance)
                        # eg:  5: {' 5: {'1015h', '827h'}', '827h'}
    for line in mfile: 
        
        #                                                    S       0       -5001_l;-8805_h;-12869_h;-25834_l;-47306_l;38133_l;
        if line[0]!="S": continue # not a compacted fact
        line=line.strip().split()
        compactedfact_id = line[1]
        facts[compactedfact_id]=[] # stored as a set to conserve the variants order, needed for computing AC.
        for value in line[2].strip().split(';')[:-1]:
            facts[compactedfact_id].append(value)
            snp_id=get_left_clean_snp(value) # from ' -587h' to '587h'
            if snp_id[:-1] not in snp_to_fact_id:
                snp_to_fact_id[snp_id[:-1]]=set()
            snp_to_fact_id[snp_id[:-1]].add(compactedfact_id)
    mfile.close()
    return facts, snp_to_fact_id
    

query_sign = lambda s : "-" if s[0]=="-" else "+"   

def compatibles(raw_fact,compacted_fact):
    """"given a raw fact, detects if it is totaly included into the compacted_fact[i]
    we already know that the 2 facts share at least a snp
    This methods detects if there exists or not a snp, with distinct alleles in the two facts (eg 1000h in one fact, 1000l in the other). 
    Returns false in this case, True else. 
    When True is returned, a second boolean value is returned. 
    It is set to True of the raw fact and the compacted_fact are in the same orientation, else it is set to false
    """

    ## index query facts:
    query_fact_to_hl = {}
    query_fact_to_pm = {}
    for variant in raw_fact:
        id_only = variant.lstrip('-')[:-1]
        query_fact_to_pm[id_only] = query_sign(variant)
        query_fact_to_hl[id_only] = variant[-1]

    ## checks that all variants from the raw fact are included in the compacted_fact: 
    for variant in raw_fact:
        if variant not in compacted_fact and variant.lstrip('-') not in compacted_fact: return False, None

    ## check compacted fact vs query fact
    same_direction  = False
    direction_set   = False
    for variant in compacted_fact:
        id_only = variant.lstrip('-')[:-1]
        if id_only not in query_fact_to_hl: continue
        sign = query_sign(variant)
        if variant[-1] != query_fact_to_hl[id_only]: return False, None
        # check sign
        if not direction_set: 
            direction_set = True
            same_direction =  (sign == query_fact_to_pm[id_only])
        else:
            if same_direction and sign != query_fact_to_pm[id_only]: return False, None
            if not same_direction and sign == query_fact_to_pm[id_only]: return False, None
    assert (direction_set)
    return True, same_direction
    

def get_compatible_facts(text_raw_fact, compacted_facts, snp_to_fact_id):
    """
    given the text of a raw fact, returns all compacted_facts ids in which this raw fact occurs
    add the sign of each of the compacted_facts id in which this raw fact occurs. Negative if it is not the same orientation as the text_raw_fact, else positive
    Example:
        * text_raw_fact: -10000l_0;92837_h;
        * compacted_facts: ... '29731': {'31938l', '499h'} ...  (factid -> list of non oriented alleles)
        * snp_to_fact_id: ... '31938': {'29731'} ...            (snp id non oriented -> list of fact ids in which the snp occurs)
    """
    subcompacedfacts=set()                      # Stores all facts potentially compatible with any of the snp of the text_raw_fact
    raw_fact_snps = set()                       # Stores the oriented alleles of the text_raw_fact. 
    result = set()                              # Stores the id of the compacted facts that are compatible with the input text_raw_fact
    text_raw_fact=text_raw_fact.rstrip(';')     #Avoids an empty value when splitting with ';'
    for oriented_allele in text_raw_fact.split(";"):
        snp_id_only=get_left_clean_snp(oriented_allele).split("_")[0][:-1]      # get the snp id non oriented
        if snp_id_only in snp_to_fact_id: # the snp may be absent in case it was removed by the sequence concatenation process. 
            subcompacedfacts=subcompacedfacts.union(snp_to_fact_id[snp_id_only])    # fill the subcompacedfacts with all facts in which the snp id non oriented occurs. 
            raw_fact_snps.add(oriented_allele.split("_")[0])
    for j in subcompacedfacts:
        cmpt,same_orientation = compatibles(raw_fact_snps, compacted_facts[j])
        if cmpt:
            if same_orientation:    result.add(int(j))
            else:                   result.add(-int(j))
            
    return result


get_allele_id = lambda x: x.split("_")[-1]+x.split("_")[1][0]   #>SNP_higher_path_9 to "9h"
get_coverage  = lambda x: int(x.split("_")[-1])                 #C1_31 to int(31)
def index_allele_coverage(raw_disco_fa_file_name, read_set_id):
    """
    For each allele name (eg 21112l) provides its coverage in the considered read set)
    Returns a dictionary: allele_id -> coverage
    Used only by "detects_allele_coverage"
    """
    alleles_coverage = {}                    # Fora each allele, store its read coverage
    mfile = open(raw_disco_fa_file_name)
    # first variant rteated appart to recover the position of the coverage in which we are interested
    coverage_field=-1
    while True:
        comment = mfile.readline()
        if not comment: break
        if not comment.startswith(">SNP"): continue # do not deal with indels for now
        s_comment = comment.strip().split("|")
        mfile.readline()        # sequence we don't care
        for i,field_content in enumerate(s_comment):
            if field_content.startswith("C"+str(read_set_id)+"_"):
                coverage_field=i
                break
        assert coverage_field != -1, "Read set id "+str(read_set_id)+" not in "+comment
        allele_id = get_allele_id(s_comment[0])
        coverage  = get_coverage(s_comment[coverage_field])
        alleles_coverage[allele_id]=coverage
    
    
    
    # other variants
    while True:
        # >SNP_higher_path_9|P_1:30_C/T,P_2:35_T/G|high|nb_pol_2|left_unitig_length_31|right_unitig_length_66|C1_31|Q1_63|G1_0/1:398,14,458|rank_0
        # ggtgcagacaacccggcaggtgttgatgataAAGATCTGGTTAAATACGCCGATATTGGCGCGACTTACTATTTCAATAAAAACATGTCCACCTACGttgactataaaatcaacctgttggatgaagatgacagcttctacgctgccaatggcatctctaccg
        # >SNP_lower_path_9|P_1:30_C/T,P_2:35_T/G|high|nb_pol_2|left_unitig_length_31|right_unitig_length_66|C1_28|Q1_63|G1_0/1:398,14,458|rank_0
        # ggtgcagacaacccggcaggtgttgatgataAAGATCTGGTTAAATACGCCGATATTGGCGTGACTGACTATTTCAATAAAAACATGTCCACCTACGttgactataaaatcaacctgttggatgaagatgacagcttctacgctgccaatggcatctctaccg
        comment = mfile.readline()
        if not comment: break
        if not comment.startswith(">SNP"): continue # do not deal with indels for now
        comment = comment.strip().split("|")
        mfile.readline()        # sequence we don't care
        allele_id = get_allele_id(comment[0])
        coverage  = get_coverage(comment[coverage_field])
        alleles_coverage[allele_id]=coverage
        
    mfile.close()
    return alleles_coverage

def detects_allele_coverage(compacted_facts, raw_disco_file_name, read_set_id):
    """
    Given the compacted facts indexed and the raw disco output: 
    for each compacted fact, find all allele that belong to it and store its coverage
    Returns a dictionary: compacted_fact_id -> allele_weight
    """
    alleles_coverage = index_allele_coverage(raw_disco_file_name, read_set_id)
    compacted_fact_allele_weight = {}              # For each compacted fact id, stores its weight
    for fact_id, fact_value in compacted_facts.items():
        compacted_fact_allele_weight[fact_id] = []

        for allele_id in fact_value: #{'10540l', '4734l', '29633h'}
            allele_coverage = alleles_coverage[get_left_clean_snp(allele_id)]
            compacted_fact_allele_weight[fact_id].append(allele_coverage)
    return compacted_fact_allele_weight

def detects_facts_coverage(compacted_facts, snp_to_fact_id, raw_facts_file_name):
    """
    Given the compacted facts indexed and the raw phasing information: for each compacted fact, find all facts that belong to it and compute its estimated coverage
    Returns a dictionary: compacted_fact_id -> weight
    """
    compacted_fact_weight = {}              # For each compacted fact id (int), stores its weight
    mfile = open(raw_facts_file_name)
    for line in mfile.readlines():
        if line[0]=="#" : continue          # comment
        coverage = int(line.strip().split("=>")[-1]) # -10011l_0;13979l_-57;21112l_-22;19270l_-14; => 4
        line=line.strip().split("=>")[0].split()     # remove coverage and split into two facts if needed
        for rawfact in line:
            if len(rawfact)==0: continue
            # detects all compacted_facts in which the rawfact occurs: 
            matching_compacted_fact_ids = get_compatible_facts(rawfact, compacted_facts, snp_to_fact_id)
            for matching_compacted_fact_id in matching_compacted_fact_ids: 
                matching_compacted_fact_id = abs(matching_compacted_fact_id) # remove the useless sign here
                if matching_compacted_fact_id not in compacted_fact_weight: 
                    compacted_fact_weight[matching_compacted_fact_id]=0
                compacted_fact_weight[matching_compacted_fact_id]+=coverage         
    mfile.close()
    return compacted_fact_weight       

def print_facts(phasing_file,compacted_fact_weight, compacted_fact_allele_weight):
    cpt=0
    mfile=open(phasing_file)
    for line in mfile:
        if line[0] != "S" : continue
        #S       0       -5001l;-8805h;-12869h;-25834l;-47306l;38133l;
        line=line.strip()
        compacted_fact_id=int(line.split()[1])
        fact_weight=0
        # assert (compacted_fact_id in compacted_fact_weight)
        fact_weight = compacted_fact_weight[compacted_fact_id]

        # assert(str(compacted_fact_id) in compacted_fact_allele_weight)

        alleles_weight = compacted_fact_allele_weight[str(compacted_fact_id)]
        



        print(f"{str(line)}\tFC:i:{str(fact_weight)}\tmin:{min(alleles_weight)}\tmax:{max(alleles_weight)}\tmean:{sum(alleles_weight)/len(alleles_weight)}\tAC:{';'.join(str(e) for e in alleles_weight)};")
        # print(line+"\tFC:i:"+str(fact_weight)+"\tRC:i:"+str(alleles_weight[1])+"\tmin:i:"+str(alleles_weight[0])+"\tmax:i:"+str(alleles_weight[1])+"\tmean:i:"+str(alleles_weight[2]))  ## RC is max ([1]) as we take the max (17/02/2020)
        cpt+=1
    sys.stderr.write(str(cpt)+" facts written\n")
    mfile.close()

def add_pair_edge (pair_edges, left, right, abundance):
    """
    Add a new pair to all stored pairs. Creates cells and dictionary if needed
    Stores a canonical representation of pairs (smallest id left)
    """
    if left > right: left, right = right, left # swap values if left > right, in order to avoid a -> b and -b -> -a
    if left not in pair_edges: pair_edges[left]={}
    if right not in pair_edges[left]: pair_edges[left][right]=0
    pair_edges[left][right]+=abundance


def detects_pairs_of_linked_compacted_paths(compacted_facts, snp_to_fact_id, raw_facts_file_name, fact_overlaps):
    """
    given the compacted facts indexed and the raw phasing information, detects pairs of facts that are co-mapped by at least one pair of paired non compacted facts
    returns a dictionary compacted_fact_id -> {compacted_fact_id} such that the key is smaller than all values.
    """
    pair_edges = {}                         # For each "left" (arbitrary) compacted fact (key) link to a dictionnary right compacted fact -> number of occurrences 
    mfile = open(raw_facts_file_name)
    for line in mfile.readlines():
        if line[0]=="#" : continue          # comment
        abundance = int(line.strip().split("=>")[1])
        line=line.strip().split("=>")[0]    # remove coverage
        line=line.strip().split()        # two pairs
        if len(line)<2: continue            # we consider only pairs of facts
        # for the first fact, detects all compacted_facts in which it occurs: 
        left_compacted_facts = get_compatible_facts(line[0], compacted_facts, snp_to_fact_id)
        # for the second fact, detects all compacted_facts in which it occurs: 
        right_compacted_facts = get_compatible_facts(line[1], compacted_facts, snp_to_fact_id)
        ### if compacted facts matched, make all pairs
        if len(left_compacted_facts)>0 and len(right_compacted_facts)>0:
            for left_compacted_fact_id in left_compacted_facts:
                for right_compacted_fact_id in right_compacted_facts:
                    if left_compacted_fact_id in fact_overlaps and right_compacted_fact_id in fact_overlaps[left_compacted_fact_id]:
                        continue    # this pair only retreives two facts that overlap
                    if right_compacted_fact_id in fact_overlaps and left_compacted_fact_id in fact_overlaps[right_compacted_fact_id]:
                        continue    # this pair only retreives two facts that overlap

                    add_pair_edge(pair_edges, left_compacted_fact_id, right_compacted_fact_id, abundance)

    mfile.close()
    return pair_edges 
    
def print_pair_edges_gfa_style(pair_edges, occurrence_min=1):
    cpt=0
    for left_fact_id in pair_edges:
        left_sign = ''
        if left_fact_id<0:  left_sign = '-'
        else:               left_sign = '+'
        abs_left_fact_id = abs(left_fact_id)
        for right_fact_id in pair_edges[left_fact_id]:
            # the right sign is reversed as pairend reads map -->  <--. Hence if right is forward on a fact we have to reverse this fact
            right_sign = ''
            if right_fact_id<0:     right_sign = '+'
            else:                   right_sign = '-'
            abs_right_fact_id = abs(right_fact_id)
            
            
            
            if pair_edges[left_fact_id][right_fact_id]>=occurrence_min:
                cpt+=1
                print ("L\t"+str(abs_left_fact_id)+"\t"+left_sign+"\t"+str(abs_right_fact_id)+"\t"+right_sign+"\t0M\tFC:i:"+str(pair_edges[left_fact_id][right_fact_id]))
    sys.stderr.write(str(cpt)+" paired fact written\n")
    
def print_facts_overlaps(phasing_file):
    fact_overlaps={}
    mfile=open(phasing_file)
    cpt=0
    for line in mfile:
        if line[0] != "L" : continue
        cpt+=1
        #L       17012   -       23084   +       2M
        line=line.strip()
        print(line)
        # store pairs: 
        source = int(line.split()[1])
        target = int(line.split()[3])
        if source > target: 
            tmp = source
            source = target
            target = tmp
        if source not in fact_overlaps: 
            fact_overlaps[source]=set()
        fact_overlaps[source].add(target)

    sys.stderr.write(str(cpt)+" facts overlaps written\n")
    return fact_overlaps
    
    
    
    
    

def detects_pairs_of_edges_sharing_snp(compacted_facts, snp_to_fact_id):
    """ 
    detects which facts share at least a snp id with incompatible h/l
    returns a dictionary fact_id -> set(fact_ids) (key is lower than any fact in the value)
    """
    facts_shared_snps = {}
    for key,values in compacted_facts.items():
        for oriented_allele in values:
            snp_id_only=get_left_clean_snp(oriented_allele).split("_")[0][:-1]      # get the snp id non oriented
            horl = get_left_clean_snp(oriented_allele).split("_")[0][-1]            # get the path 'h' or 'l' of the SNP
            # print("snp_id_only",snp_id_only)
            if snp_id_only in snp_to_fact_id: # the snp may be absent in case it was removed by the sequence concatenation process. 
                for fact_id in snp_to_fact_id[snp_id_only]:
                    if int(fact_id)<=int(key): 
                        continue
                    fact = compacted_facts[fact_id]
                    for snp_id in fact: 
                        ### checks that h or l values are disctincts between the two facts 
                        if get_left_clean_snp(snp_id).split("_")[0][:-1] == snp_id_only and get_left_clean_snp(snp_id).split("_")[0][-1]!=horl: # TO VALIDATE 9/9/2019
                            if key not in facts_shared_snps: facts_shared_snps[key] = set()    
                            facts_shared_snps[key].add(fact_id)
    # for key, value in facts_shared_snps.items():
    #     print(key,value)
    return facts_shared_snps
            
def print_pairs_of_edges_sharing_snp(facts_shared_snps):
    cpt=0
    for key, values in facts_shared_snps.items():
        for value in values: 
            print("L\t"+key+"\t+\t"+value+"\t+\t-2M")
            cpt+=1
    sys.stderr.write(str(cpt)+" pairs of facts sharing at least one snp written\n")


def main (phasing_file,raw_facts_file_name, raw_disco_file_name, read_set_id):
    sys.stderr.write("#INDEX FACTS\n")
    compacted_facts, snp_to_fact_id = set_indexes_from_gfa(phasing_file)
    print(f"{compacted_facts['4301']}")
    
    sys.stderr.write("#COMPUTE THE COMPACTED FACT COVERAGES\n")
    compacted_fact_weight = detects_facts_coverage(compacted_facts, snp_to_fact_id, raw_facts_file_name)
    
    sys.stderr.write("#COMPUTE THE COMPACTED FACT ALLELE COVERAGES\n")
    compacted_fact_allele_weight=detects_allele_coverage(compacted_facts, raw_disco_file_name, read_set_id)
    
    sys.stderr.write("#PRINT COMPACTED FACTS \n")
    print_facts(phasing_file,compacted_fact_weight, compacted_fact_allele_weight)
    
    sys.stderr.write("#PRINT COMPACTED FACT OVERLAPS \n")
    fact_overlaps=print_facts_overlaps(phasing_file)

    sys.stderr.write("#COMPUTE PAIRS OF COMPACTED FACT GRAPH\n")
    pair_edges=detects_pairs_of_linked_compacted_paths(compacted_facts, snp_to_fact_id, raw_facts_file_name,fact_overlaps)
    
    sys.stderr.write("#PRINT EDGES OF COMPACTED FACT GRAPH\n")
    print_pair_edges_gfa_style(pair_edges) 
    
    sys.stderr.write("#COMPUTE THE FACTS SHARING AT LEAST ONE SNP\n")
    facts_shared_snps = detects_pairs_of_edges_sharing_snp(compacted_facts, snp_to_fact_id)
    print_pairs_of_edges_sharing_snp(facts_shared_snps)
    
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) # compacted_facts.gfa phased_alleles_read_set_id_1.txt discoRes_k_31_c_2_D_0_P_3_b_2_coherent.fa 1

    
    