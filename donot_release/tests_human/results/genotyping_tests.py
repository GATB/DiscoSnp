import math

import sys

def genotype_simple_model(ch, cl, err, prior_het):
    lik0 = pow(1-err,ch)*pow(err,cl) # homozygous higher (0/0)
    lik1 = pow(1-err,cl)*pow(err,ch) # homozygous higher (1/1)
    lik2 = pow(0.5,ch+cl)            # heterozygous (0/1)
    
    prob0 = lik0*(1-prior_het)/2
    prob1 = lik1*(1-prior_het)/2
    prob2=lik2*prior_het
    
    if prob0>=prob1 and prob0>=prob2: return "1" #0/0
    if prob1>=prob0 and prob1>=prob2: return "1" #1/1
    return "0" #"0/1";


def genotype_threshold_model(ch, cl, ratio):
    if ch==cl: return "0" # avoid the 0 0 case
    f=ch/float(ch+cl)
    if f>=ratio and f<=(1-ratio): return "0"
    return "1"




# Compare mapped positions with reference positions
def analyseDiscoPaths(discores):
    filin = open(discores, 'r')
    while True:
        com1=filin.readline()
        if not com1: break
        seq1=filin.readline()
        com2=filin.readline()
        seq2=filin.readline()

        # >SNP_higher_path_99995|P_1:30_C/T|high|nb_pol_1|C1_0|C2_38|G1_1/1|G2_0/0|rank_1.00000
#         ATTAACCCCACACTGTGATATCTTCTTTATCGAAAGAGGTAATCTGTAGTGCAACTAAAAA
#         >SNP_lower_path_99995|P_1:30_C/T|high|nb_pol_1|C1_42|C2_0|G1_1/1|G2_0/0|rank_1.00000
#         ATTAACCCCACACTGTGATATCTTCTTTATTGAAAGAGGTAATCTGTAGTGCAACTAAAAA
            

        # COVERAGES
        text_coverage1=com1.split("C1")[1].split("|G1")[0] #_0|C2_38
        nb_set=len(text_coverage1.split("|"))
        text_coverage2=com2.split("C1")[1].split("|G1")[0] #_42|C2_0
        
        genotype_text=""
        for set_id in range(nb_set):
            # genotype_text+="_"+genotype_simple_model(int(text_coverage1.split("|")[set_id].split("_")[1]), int(text_coverage2.split("|")[set_id].split("_")[1]), 0.01, 0.33)
            genotype_text+="_"+genotype_threshold_model(int(text_coverage1.split("|")[set_id].split("_")[1]), int(text_coverage2.split("|")[set_id].split("_")[1]), 0.2)
        
        print com1[:-1]+genotype_text
        print seq1,
        print com2[:-1]+genotype_text
        print seq2,
        
        
            
analyseDiscoPaths(sys.argv[1])