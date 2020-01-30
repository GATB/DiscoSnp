import sys
import K3000_common as kc

def store_fact_extreme_snp_ids(gfa_file_name):
    """ Given a gfa file, store the id of each extreme SNP (without h or l)
    """
    mfile = open(gfa_file_name)
    leftmost_snp_to_fact_id={}   # each SNP id (with order and without h or l) is linked to some compacted facts id
    rightmost_snp_to_fact_id={}   # each SNP id (with order and without h or l) is linked to some compacted facts id
                        # eg: '-1015': {8, 5, 6, 7}

    for line in mfile: 
        
        #                                                    S       10      24617l;10033l;-11833l;  RC:i:55
        if line[0]!="S": continue # not a compacted fact
        line=line.strip().split()
        compactedfact_id = int(line[1])
        
        #left
        leftmost_snp_id=int(line[2].strip().split(';')[0][:-1])
        if leftmost_snp_id not in leftmost_snp_to_fact_id:
            leftmost_snp_to_fact_id[leftmost_snp_id] = set()
        leftmost_snp_to_fact_id[leftmost_snp_id].add(compactedfact_id)
        
        #right
        rightmost_snp_id=int(line[2].strip().split(';')[-2][:-1])# -1 is empty
        if rightmost_snp_id not in rightmost_snp_to_fact_id:
            rightmost_snp_to_fact_id[rightmost_snp_id] = set()
        rightmost_snp_to_fact_id[rightmost_snp_id].add(compactedfact_id)
    
    mfile.close()
    return leftmost_snp_to_fact_id, rightmost_snp_to_fact_id
    

def get_uppercase_sequence(sequence):
    """given a sequence lBr with l and r in acgt and B in ACGT, return B"""
    res=""
    for l in sequence:
        if l>='A' and l<='Z':
            res+=l
    return res
 

def store_remarkable_kmers(fa_file_name, k, leftmost_snp_to_fact_id, rightmost_snp_to_fact_id):
    """
    Given a disco output, store for SNPs the remarkable (k-1)mers.
    
    Remarkable (k-1)mers are
    LO (Left Out) Leftmost (k-1)mer on the left unitig
    LI (Left In) Leftmost (k-1)mer not on the SNP (corresponds to the first k-1 upper case characters)
    RI (Right In) Rightmost (k-1)mer not on the SNP (corresponds to the last k-1 upper case characters)
    RO (Right Out) Rightmost (k-1)mer on the right unitig
    
    We store only LO and LI for SNPs that are a left most SNP of at least a fact
    We store only RO and RI for SNPs that are a right most SNP of at least a fact
    """
    
    LOs = {} # key: (k-1)mer, value: set of fact ids
    LIs = {} # key: (k-1)mer, value: set of fact ids
    RIs = {} # key: (k-1)mer, value: set of fact ids
    ROs = {} # key: (k-1)mer, value: set of fact ids
    
    mfile = open(fa_file_name)
    while True:
        line1 = mfile.readline()
        if not line1: break
        line1 = line1.strip()               #>SNP_higher_path_9995|P_1:30_A/C|high|nb_pol_1|left_unitig_length_1|right_unitig_length_93|C1_39|Q1_63|G1_0/1:494,15,574|rank_0
        line2 = mfile.readline().strip()    #tGCGCGTCTCCGGCCTGAAAAAGCTGTCCGTAAACGGTACAGATAGCAATCCCCAATCGGGAaccgtctactgtagccagcgccggaatataatccgcgactttaccctgaccaatgagcggccgcacttgccgcaagatgttttctaaaattgc
        mfile.readline().strip()            #>SNP_lower_path_9995|P_1:30_A/C|high|nb_pol_1|left_unitig_length_1|right_unitig_length_93|C1_35|Q1_63|G1_0/1:494,15,574|rank_0                                 no need to store
        mfile.readline().strip()            #tGCGCGTCTCCGGCCTGAAAAAGCTGTCCGTCAACGGTACAGATAGCAATCCCCAATCGGGAaccgtctactgtagccagcgccggaatataatccgcgactttaccctgaccaatgagcggccgcacttgccgcaagatgttttctaaaattgc    no need to store
        
        if not line1.startswith(">SNP"): continue #not a SNP
        
        snp_id = int(line1.split("|")[0].split('_')[-1])
        
        # stored = False # DEBUG
        central_sequence = None
        
        ### FORWARD CASES ###
        # This SNP (forward) is the starting of at least one compacted fact. Thus we store its remakable left (k-1)mers and we associate them to the corresponding facts
        if snp_id in leftmost_snp_to_fact_id:
            # stored = True
            LO = line2[:k-1].upper()            # get the first (k-1)mer 
            central_sequence = get_uppercase_sequence(line2)
            LI = central_sequence[:k-1]
            # association (k-1)mer -> facts
            if LO not in LOs: LOs[LO] = set()
            LOs[LO] = LOs[LO].union(leftmost_snp_to_fact_id[snp_id])
            if LI not in LIs: LIs[LI] = set()
            LIs[LI] = LIs[LI].union(leftmost_snp_to_fact_id[snp_id])
            # print("heyL",str(snp_id)," hoL ", leftmost_snp_to_fact_id[snp_id])
            
            
        
        # This SNP (forward) is the ending of at least one compacted fact. Thus we store its remakable right (k-1)mers and we associate them to the corresponding facts
        if snp_id in rightmost_snp_to_fact_id:
            # stored = True
            RO = line2[-k+1:].upper()           # get the last (k-1)mer
            if not central_sequence:
                central_sequence = get_uppercase_sequence(line2)
            RI = central_sequence[-k+1:]
            # association (k-1)mer -> facts
            if RO not in ROs: ROs[RO] = set()
            ROs[RO] = ROs[RO].union(rightmost_snp_to_fact_id[snp_id])
            if RI not in RIs: RIs[RI] = set()
            RIs[RI] = RIs[RI].union(rightmost_snp_to_fact_id[snp_id])
            # print("heyR",str(snp_id)," hoR ", rightmost_snp_to_fact_id[snp_id])
            
        ### REVERSE CASES ###
        # In this cses a facts starts by the reverse of a SNP, eg, -10321. In this case it is sufficient to reverse complement the sequence of the SNP and to have exactly the same treatment as in the forward case. 
        # This SNP (reverse) is the starting of at least one compacted fact. Thus we store its remakable left (k-1)mers and we associate them to the corresponding facts
        if -snp_id in leftmost_snp_to_fact_id:
            # stored = True
            line2 = kc.get_reverse_complement(line2)
            LO = line2[:k-1].upper()            # get the first (k-1)mer 
            central_sequence = get_uppercase_sequence(line2)
            LI = central_sequence[:k-1]
            # association (k-1)mer -> facts
            if LO not in LOs: LOs[LO] = set()
            LOs[LO] = LOs[LO].union(leftmost_snp_to_fact_id[-snp_id])
            if LI not in LIs: LIs[LI] = set()
            LIs[LI] = LIs[LI].union(leftmost_snp_to_fact_id[-snp_id])
            # print("heyL",str(-snp_id)," hoL ", leftmost_snp_to_fact_id[-snp_id])
            
            
        
        # This SNP (reverse) is the ending of at least one compacted fact. Thus we store its remakable right (k-1)mers and we associate them to the corresponding facts
        if -snp_id in rightmost_snp_to_fact_id:
            # stored = True
            line2 = kc.get_reverse_complement(line2)
            RO = line2[-k+1:].upper()           # get the last (k-1)mer
            if not central_sequence:
                central_sequence = get_uppercase_sequence(line2)
            RI = central_sequence[-k+1:]
            # association (k-1)mer -> facts
            if RO not in ROs: ROs[RO] = set()
            ROs[RO] = ROs[RO].union(rightmost_snp_to_fact_id[-snp_id])
            if RI not in RIs: RIs[RI] = set()
            RIs[RI] = RIs[RI].union(rightmost_snp_to_fact_id[-snp_id])
            # print("heyR",str(-snp_id)," hoR ", rightmost_snp_to_fact_id[-snp_id])
        # if stored:
        #     print(line2)
        #     print(LOs)
        #     print(LIs)
        #     print(RIs)
        #     print(ROs)
            # sys.exit(0)
    mfile.close()
    return LOs, LIs, RIs, ROs
    

def print_link_facts(LOs, LIs, RIs, ROs, k, rightmost_snp_to_fact_id,fa_file_name):
    """ given the association (k-1)mers -> leftfacts, we may derive the links between facts
    1. RIA == LOB and ROA == LIB           -> A+ -> B+
    2. RIA == rc(ROB) and ROA == rc(RIB)   -> A+ -> B-
    3. LOA == rc(LIB) and LIA == rc(LOB)   -> A- -> B+
    4. LOA == RIB and LIA == ROB           -> A- -> B- (eq B+ -> A+)
    
    In this function we traverse again all SNPS, 
     for those that are the END of a fact : 
        check if their RI and RO may lead to one of the previous link
     for those whose reverse complement is the END of a fact : 
        Reverse the sequence and then check if their RI and RO may lead to one of the previous link
    """
    mfile = open(fa_file_name)
    while True:
        line1 = mfile.readline()
        if not line1: break
        line1 = line1.strip()               #>SNP_higher_path_9995|P_1:30_A/C|high|nb_pol_1|left_unitig_length_1|right_unitig_length_93|C1_39|Q1_63|G1_0/1:494,15,574|rank_0
        line2 = mfile.readline().strip()    #tGCGCGTCTCCGGCCTGAAAAAGCTGTCCGTAAACGGTACAGATAGCAATCCCCAATCGGGAaccgtctactgtagccagcgccggaatataatccgcgactttaccctgaccaatgagcggccgcacttgccgcaagatgttttctaaaattgc
        mfile.readline().strip()            #>SNP_lower_path_9995|P_1:30_A/C|high|nb_pol_1|left_unitig_length_1|right_unitig_length_93|C1_35|Q1_63|G1_0/1:494,15,574|rank_0                                 no need to store
        mfile.readline().strip()            #tGCGCGTCTCCGGCCTGAAAAAGCTGTCCGTCAACGGTACAGATAGCAATCCCCAATCGGGAaccgtctactgtagccagcgccggaatataatccgcgactttaccctgaccaatgagcggccgcacttgccgcaagatgttttctaaaattgc    no need to store
        
        if not line1.startswith(">SNP"): continue #not a SNP
        
        snp_id = int(line1.split("|")[0].split('_')[-1])
        
        # This SNP (forward) is the ending of at least one compacted fact. Thus we store its remakable right (k-1)mers and we associate them to the corresponding facts
        if snp_id in rightmost_snp_to_fact_id:
            ROA = line2[-k+1:].upper()           # get the last (k-1)mer
            central_sequence = get_uppercase_sequence(line2)
            RIA = central_sequence[-k+1:]
            
            # CASE 1.
            if RIA in LOs and ROA in LIs:
                set_of_left_facts_compatible = LOs[RIA].intersection(LIs[ROA])                                                  # select all facts ids both whose LO == RI and LI == RO                 -> A+ -> B+
                for left_fact_id in set_of_left_facts_compatible:
                    for right_fact_id in rightmost_snp_to_fact_id[snp_id]:
                        print ("L\t"+str(left_fact_id)+"\t+\t"+str(right_fact_id)+"\t+\t"+str(-1)+"M")                          # -1 enables to detect those links

            # CASE 2.
            if kc.get_reverse_complement(RIA) in ROs and kc.get_reverse_complement(ROA) in RIs:
                set_of_left_facts_compatible = ROs[kc.get_reverse_complement(RIA)].intersection(RIs[kc.get_reverse_complement(ROA)])  # select all facts ids both whose RIA == rc(ROB) and ROA == rc(RIB)     -> A+ -> B-
                for left_fact_id in set_of_left_facts_compatible:
                    for right_fact_id in rightmost_snp_to_fact_id[snp_id]:
                        print ("L\t"+str(left_fact_id)+"\t+\t"+str(right_fact_id)+"\t-\t"+str(-1)+"M")                          # -1 enables to detect those links
                    
            # CASE 3. & 4. -> they are symetrical:
            # Cases 1. & 2. : detects A -> B and A -> B_, the other cases are
            # CASE 3. A_ -> B  will be detected when traversing B, detecting then A_ (Case 2.)
            # CASE 4. B -> A   will be detected when traversing Bn detecting then A  (Case 1.)

    mfile.close()

def print_original_gfa(gfa_file_name):
    mfile = open(gfa_file_name)
    print ("H\t#################")
    print ("H\t# GFA of variants")
    print ("H\t#################")
    print ("H\t# Nodes are (compacted) facts with their read mapping coverage. Eg. \"S	34408h;10036h;-3759h; SP:0_450;373_575;480_743; BP:0_86;17_90;35_61;    FC:i:58 RC:i:29\".")
    print ("H\t#  * field SP stands for \"Sequence Position\". A bit useless in this file. It indicates for each allele of the fact its starting and ending position in the ACGT sequence, not shown here")
    print ("H\t#  * field BP stands for \"Bubble Position\". For each allele of the fact it indicates:")
    print ("H\t#     - first: the relative position of first nucleotide of the bubble with respect to the position of the last nucleotide of the bubble of the previous allele. This value is equal to zero for the first allele")
    print ("H\t#     - second: the length of the bubble of the allele") 
    print ("H\t#  * field FC is the coverage of the fact, as provided by the total number of reads that phased at least two alleles of the fact")
    print ("H\t#  * field RC is the coverage of the fact, as provided by the min of the read coverage from all independent allele computed coverages")
    print ("H\t# Four types of edges:")
    print ("H\t#   1. Overlap between facts. These links have an overlap length >0. Eg, \"L	1	-	29384	+	3M\", with:")
    print ("H\t#       \"S	1	10011l;23229h;-21935l;-8929l;-24397l;10011h;	RC:i:24\", and")
    print ("H\t#       \"S	29384	21935l;-23229h;-10011l;24397l;-23229l;-25549h;-10011h;	RC:i:43\".")
    print ("H\t#   2. Facts linked by paired end reads.  Eg \"L	10735	+	29384	+	0M	FC:i:5\".")
    print ("H\t#       These links are non directed and do no validate the facts orientation. The coverage indicates the number of pairend read linking the two facts")
    print ("H\t#       These links have an overlap of length 0.")
    print ("H\t#   3. Facts linked by unitigs. The unitig finishing a fact overlaps the unitig starting another fact. Eg \"L	19946	+	11433	+	-1M\".")
    print ("H\t#       These links are directed and validate the facts orientation. ")
    print ("H\t#       These links have an overlap of length -1.")
    print ("H\t#   4. Facts sharing at least one SNP identifier.")
    print ("H\t#       These links have an overlap of length -2.")
    
    
    
    for l in mfile: 
        print(l,end='')
    mfile.close()
    



def main (gfa_file_name, fa_file_name):
    sys.stderr.write("#Store extreme left and right SNP id for each fact\n")
    leftmost_snp_to_fact_id, rightmost_snp_to_fact_id   = store_fact_extreme_snp_ids(gfa_file_name)
    k                                                   = kc.determine_k(fa_file_name)
    sys.stderr.write("#Store remarkable kmers for each such SNP\n")
    LOs, LIs, RIs, ROs                                  = store_remarkable_kmers(fa_file_name,k,leftmost_snp_to_fact_id, rightmost_snp_to_fact_id)

    sys.stderr.write("#Print original gfa\n")
    print_original_gfa(gfa_file_name)
    sys.stderr.write("#Print links\n")
    print_link_facts(LOs, LIs, RIs, ROs, k, rightmost_snp_to_fact_id, fa_file_name)
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
