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
        facts[compactedfact_id]=set()
        for value in line[2].strip().split(';')[:-1]:
            snp_id=get_left_clean_snp(value) # from ' -587h' to '587h'
            facts[compactedfact_id].add(snp_id)
            if snp_id[:-1] not in snp_to_fact_id:
                snp_to_fact_id[snp_id[:-1]]=set()
            snp_to_fact_id[snp_id[:-1]].add(compactedfact_id)
    mfile.close()
    return facts, snp_to_fact_id
    
    

def compatibles(raw_fact,i,compacted_facts):
    """"given a raw fact, detects if it is (at least partly) included into the compacted_fact[i]
    we already know that the 2 facts share at least a snp
    This methods detects if there exists or not a snp, with distinct alleles in the two facts (eg 1000h in one fact, 1000l in the other). 
    Returns false in this case, True else. 
    """
    # print (compacted_facts[i])
    # print(raw_fact)
    compacted_fact_i_dict = {} # snps id ->  h or l
    for snp in compacted_facts[i]:
        snp=snp.lstrip('-')     # remove the sign
        compacted_fact_i_dict[snp[:-1]]=snp[-1]
    
    nb_shared=0
    for snp in raw_fact:
        snp=snp.lstrip('-')                                 # remove the sign
        if snp[:-1] in compacted_fact_i_dict:
            if compacted_fact_i_dict[snp[:-1]]!=snp[-1]:    # distinct alleles 
                # print ("            INCOMPATIBLE        ", snp, compacted_facts[i])
                # sys.exit(0)
                return False
            else :
                nb_shared+=1
        # else: return False      # uncomment if one needs the raw fact to be fully inlcuded in the compacted_facts[i]
    if nb_shared>1: # at least two shared alleles to link two facts.  # note that if test >1 one may avoid lot of computations upstream.
        return True
    if len(raw_fact)==nb_shared or len(compacted_facts[i])==nb_shared: # the raw fact or the compacted fact is fully mapped
        return True
    return False
    

def get_compatible_facts(text_raw_fact, compacted_facts, snp_to_fact_id):
    """
    given the text of a raw fact, returns all compacted_facts ids in which this raw fact occurs
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
        # print("oriented_allele -"+oriented_allele+"-")
        snp_id_only=get_left_clean_snp(oriented_allele).split("_")[0][:-1]      # get the snp id non oriented
        # print("snp_id_only",snp_id_only)
        if snp_id_only in snp_to_fact_id: # the snp may be absent in case it was removed by the sequence concatenation process. 
            subcompacedfacts=subcompacedfacts.union(snp_to_fact_id[snp_id_only])    # fill the subcompacedfacts with all facts in which the snp id non oriented occurs. 
            # print("subcompacedfacts", subcompacedfacts)
            raw_fact_snps.add(oriented_allele.split("_")[0])
    # print ("raw_fact_snps", raw_fact_snps)
    for j in subcompacedfacts:
        if compatibles(raw_fact_snps, j, compacted_facts):
            result.add(int(j))
    
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
    Given the compacted facts indexed and the raw disco output: for each compacted fact, find all allele that belong to it and compute its estimated coverage (average, min, max)
    Returns a dictionary: compacted_fact_id -> allele_weight
    TODO. In fact we use only the "min" value. Thus, no need to compute and to store the mean and max values. 
    """
    alleles_coverage = index_allele_coverage(raw_disco_file_name, read_set_id)
    compacted_fact_allele_weight = {}              # For each compacted fact id, stores its weight
    for fact_id, fact_value in compacted_facts.items():

        min = sys.maxsize
        max = -1
        nb = 0
        sum=0
        for allele_id in fact_value: #{'10540l', '4734l', '29633h'}
            allele_coverage = alleles_coverage[allele_id]
            if allele_coverage<min: min=allele_coverage
            if allele_coverage>max: max=allele_coverage
            sum+=allele_coverage
            nb+=1
        compacted_fact_allele_weight[fact_id] = [min,max,sum/float(nb)] #min, max, mean
        # print(compacted_fact_allele_weight[fact_id])
    return compacted_fact_allele_weight


def detects_facts_coverage(compacted_facts, snp_to_fact_id, raw_facts_file_name):
    """
    Given the compacted facts indexed and the raw phasing information: for each compacted fact, find all facts that belong to it and compute its estimated coverage
    Returns a dictionary: compacted_fact_id -> weight
    """
    compacted_fact_weight = {}              # For each compacted fact id, stores its weight
    mfile = open(raw_facts_file_name)
    for line in mfile.readlines():
        if line[0]=="#" : continue          # comment
        coverage = int(line.strip().split("=>")[-1]) # -10011l_0;13979l_-57;21112l_-22;19270l_-14; => 4
        line=line.strip().split("=>")[0].split(" ")     # remove coverage and split into two facts if needed
        for rawfact in line:
            if len(rawfact)==0: continue
            # detects all compacted_facts in which the rawfact occurs: 
            matching_compacted_fact_ids = get_compatible_facts(rawfact, compacted_facts, snp_to_fact_id)
            for matching_compacted_fact_id in matching_compacted_fact_ids: 
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
        if compacted_fact_id in compacted_fact_weight: 
            fact_weight = compacted_fact_weight[compacted_fact_id]
        alleles_weight=0
        if str(compacted_fact_id) in compacted_fact_allele_weight: 
            alleles_weight = compacted_fact_allele_weight[str(compacted_fact_id)]
        print(line+"\tFC:i:"+str(fact_weight)+"\tRC:i:"+str(alleles_weight[0]))
        cpt+=1
    sys.stderr.write(str(cpt)+" facts written\n")
    mfile.close()

def detects_pairs_of_linked_compacted_paths(compacted_facts, snp_to_fact_id, raw_facts_file_name,fact_overlaps):
    """
    given the compacted facts indexed and the raw phasing information, detects pairs of facts that are co-mapped by at least one pair of paired non compacted facts
    returns a dictionary compacted_fact_id -> {compacted_fact_id} such that the key is smaller than all values.
    """
    pair_edges = {}                         # For each "left" (arbitrary) compacted fact (key) link to a dictionnary right compacted fact -> number of occurrences 
    mfile = open(raw_facts_file_name)
    for line in mfile.readlines():
        if line[0]=="#" : continue          # comment
        line=line.strip().split("=>")[0]    # remove coverage
        line=line.strip().split(" ")        # two pairs
        if len(line)<2: continue            # we consider only pairs of facts
        # print(line)
        # for the first fact, detects all compacted_facts in which it occurs: 
        left_compacted_facts = get_compatible_facts(line[0], compacted_facts, snp_to_fact_id)
        # for the second fact, detects all compacted_facts in which it occurs: 
        right_compacted_facts = get_compatible_facts(line[1], compacted_facts, snp_to_fact_id)
        ### if compacted facts matched, make all pairs
        if len(left_compacted_facts)>0 and len(right_compacted_facts)>0:
            for left_compacted_fact_id in left_compacted_facts:
                if left_compacted_fact_id not in pair_edges:
                    pair_edges[left_compacted_fact_id]={}
                for right_compacted_fact_id in right_compacted_facts:
                    if left_compacted_fact_id in fact_overlaps and right_compacted_fact_id in fact_overlaps[left_compacted_fact_id]:
                        continue    # this pair only retreives two facts that overlap
                    if right_compacted_fact_id in fact_overlaps and left_compacted_fact_id in fact_overlaps[right_compacted_fact_id]:
                        continue    # this pair only retreives two facts that overlap
                    if right_compacted_fact_id not in pair_edges[left_compacted_fact_id]:
                        pair_edges[left_compacted_fact_id][right_compacted_fact_id]=0
                    pair_edges[left_compacted_fact_id][right_compacted_fact_id]+=1

    mfile.close()
    return pair_edges 
    
def print_pair_edges_gfa_style(pair_edges, occurrence_min=1):
    cpt=0
    for left_fact_id in pair_edges:
        for right_fact_id in pair_edges[left_fact_id]:
            if left_fact_id < right_fact_id and pair_edges[left_fact_id][right_fact_id]>=occurrence_min:
                cpt+=1
                print ("L\t"+str(left_fact_id)+"\t+\t"+str(right_fact_id)+"\t+\t0M\tFC:i:"+str(pair_edges[left_fact_id][right_fact_id]))
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

def main (phasing_file,raw_facts_file_name, raw_disco_file_name, read_set_id):
    sys.stderr.write("#INDEX FACTS\n")
    compacted_facts, snp_to_fact_id = set_indexes_from_gfa(phasing_file)

    # print(facts)
    # print(snp_to_fact_id)
    
    sys.stderr.write("#COMPUTE THE COMPACTED FACT COVERAGES\n")
    compacted_fact_weight=detects_facts_coverage(compacted_facts, snp_to_fact_id, raw_facts_file_name)
    
    sys.stderr.write("#COMPUTE THE COMPACTED FACT ALLELE COVERAGES\n")
    compacted_fact_allele_weight=detects_allele_coverage(compacted_facts, raw_disco_file_name, read_set_id)
    
    sys.stderr.write("#PRINT COMPACTED FACTS \n")
    print_facts(phasing_file,compacted_fact_weight, compacted_fact_allele_weight)
    
    sys.stderr.write("#PRINT COMPACTED FACT OVERLAPS \n")
    fact_overlaps=print_facts_overlaps(phasing_file)

    sys.stderr.write("#COMPUTE PAIRS OF COMPACTED FACT GRAPH\n")
    pair_edges=detects_pairs_of_linked_compacted_paths(compacted_facts, snp_to_fact_id, raw_facts_file_name,fact_overlaps)
    
    sys.stderr.write("#PRINT EDGES OF COMPACTED FACT GRAPH\n")
    print_pair_edges_gfa_style(pair_edges,3) 
    
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) # compacted_facts.gfa phased_alleles_read_set_id_1.txt discoRes_k_31_c_2_D_0_P_3_b_2_coherent.fa 1

    
    