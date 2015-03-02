#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt
import parse_gassst_common



####################################################################
#                         REFERENCES
####################################################################
# Store the polymorphism positions
def storeThePolymorphismPositionsGenericList(referencePositionFile, typePolymorphism):
    filin = open(referencePositionFile, 'r')
    list_reference_polymorphism={} # for each chromosome: contains a list of positions
    list_reference_genotypes={} # for each chromosome: contains a list of positions
    while 1:
        line=filin.readline() # SNP 233472 gi C G 0|0 0|1
        if not line: break
        if line.startswith(typePolymorphism):
            pos=int(line.split(" ")[1]) 
            chromosome=line.split(" ")[2].strip()
            if not chromosome in list_reference_polymorphism:                
                list_reference_polymorphism[chromosome]={}
                list_reference_genotypes[chromosome]={}
            if pos in list_reference_polymorphism[chromosome]:
                print "Warning, more than one "+typePolymorphism+" exists position "+str(pos)+" on chromosome "+chromosome
            else:
                list_reference_polymorphism[chromosome][pos]=False # Will be set to true if this SNP is detected. 
                list_reference_genotypes[chromosome][pos]=line.split(" ")[5]+" "+line.split(" ")[6] 

    return list_reference_polymorphism,list_reference_genotypes


####################################################################
#                         REFERENCES
####################################################################


    

# Compare mapped positions with reference positions
def checkMappedDiscoPaths(gassstOutFile,typePolymorphism, list_reference_polymorphism, list_reference_genotypes, list_predicted_polymorphism, threshold,list_predicted_positions, nb_good_genotype_predictions, nb_bad_genotype_predictions, genotypes_predictions):
    filin = open(gassstOutFile, 'r')
    for line in filin.readlines(): #SNP_higher_path_99995|P_1:30_C/T|high|nb_pol_1|C1_0|C2_38|G1_1/1|G2_0/0|rank_1.00000    gi|224384768|gb|CM000663.1|     0       108430783       61M          1       0  
        if line.startswith(typePolymorphism):
            event_id=line.split("|")[0].split("_")[-1]
            this_threshold=1
            if threshold>0:
                this_threshold=float(line.split("\t")[0].split("|")[-1].split("_")[-1])
            if this_threshold<threshold: 
                continue
            matching_chromosome=line.split("\t")[1].split("|")[0].strip()
            matching_position=int(line.split("\t")[3])
            good_prefix=line.split("|")[0] # if this is a discoSnp++ output, this is correct
            predicted_pol_id=int(good_prefix.split("_")[-1])
            
            
            predicted_genotype = genotypes_predictions[event_id]
                
            # avoid to map twice the same variant (both upper and lower path as the mapping is not exact)
            if predicted_pol_id in list_predicted_polymorphism and list_predicted_polymorphism[predicted_pol_id] > 0: continue
            
            span_match=len(line.split("\t")[7])
            forward=False
            if int(line.split("\t")[5])==1: forward=True
             # check if this mapping is in the list of reference stuffs
            localNbMapped=0
            if not matching_chromosome in list_reference_polymorphism: continue
            if predicted_pol_id in list_predicted_polymorphism:
                start=matching_position-1
                stop=matching_position+span_match-1
                if typePolymorphism == "INDEL":
                    localNbMapped = parse_gassst_common.findNbMachedPolymorphism_with_genotype(
                        list_reference_polymorphism[matching_chromosome],
                        list_reference_genotypes[matching_chromosome],
                        start,
                        stop,
                        [i for i in range(stop-start+1)], # fake acceptation of all mapped indel.
                        predicted_genotype,
                        typePolymorphism, 
                        forward) 
                else:
                    localNbMapped = parse_gassst_common.findNbMachedPolymorphism_with_genotype(
                    list_reference_polymorphism[matching_chromosome],
                    list_reference_genotypes[matching_chromosome],
                    start,
                    stop,
                    list_predicted_positions[predicted_pol_id], 
                    predicted_genotype,
                    typePolymorphism, 
                    forward)
            if localNbMapped>0:
                list_predicted_polymorphism[predicted_pol_id] = localNbMapped[0] # We mapped the polymorphism(s) corresponding to this predicted id. `
                nb_good_genotype_predictions+=localNbMapped[1]
                nb_bad_genotype_predictions+=localNbMapped[2]
            
    print "Good/bad genotypes: ", nb_good_genotype_predictions, nb_bad_genotype_predictions
                

def get_genotypes(discoResFile):
    res={}
    filin = open(discoResFile, 'r')
    while True: #>SNP_higher_path_99995|P_1:30_C/T|high|nb_pol_1|C1_0|C2_38|G1_1/1|G2_0/0|rank_1.00000_1_1 
        line = filin.readline ()
        if not line: break
        # GENOTYPE PREDICTION BY KISSREADS
        text_predicted_genotype=line.split("G1")[1].split("|rank")[0]
        nb_set=len(text_predicted_genotype.split("|"))
        event_id=line.split("|")[0].split("_")[-1]
        res[event_id]=[]

    
        # GENOTYPE PREDICTION BY Coverages analysis (python genotyping_tests.py discoRes_k_31_c_4_D_10_P_4_b_1_withlow_coherent.fa )
        for set_id in range(nb_set):
            res[event_id].append(int(line.split("|")[-1].split("_")[set_id+2]))
        filin.readline ()
        filin.readline ()
        filin.readline ()

    return res

####################################################################
#                         MAIN
####################################################################
def doAll(discoResFile, logFile, gassstResFile, polymorphism, threshold, verbose,generic,roc, nb_good_genotype_predictions, nb_bad_genotype_predictions, genotypes_predictions):
    
    res=parse_gassst_common.getPredictedPolymorphism(discoResFile, threshold,  polymorphism)
    list_predicted=res[0]
    nb_predicted=res[1]
    list_predicted_positions=res[2]
    res=storeThePolymorphismPositionsGenericList(logFile, polymorphism)
    list_reference=res[0]
    list_reference_genotypes=res[1]
    
    checkMappedDiscoPaths(gassstResFile,polymorphism, list_reference, list_reference_genotypes, list_predicted,threshold,list_predicted_positions, nb_good_genotype_predictions, nb_bad_genotype_predictions, genotypes_predictions)

    

def usage():
    '''Usage'''
    print "-----------------------------------------------------------------------------"
    print sys.argv[0]," : parse the gassst result on discoMore files"
    print "-----------------------------------------------------------------------------"
    print " A bubble is considered as valid iif: one (not zero, not two) of its paths map with no error the non muted genome"
    print " A reference polymorphism is found if the mapped path of a valid bubble covers it"
    print " TP: reference polymorphism covered by a valid bubble"
    print " FP: reference polymorphism non covered by a valid bubble"
    print " FN: non valid bubble or valid bubble not covering a targeted polymorphism"
    print " ------------"
    print "usage: ",sys.argv[0]," -d discoMore coherent_out_file -l log_file -g gassst_out_file [-t float] [-v] [-r] "
    print "  -d discoSnp outfile. Note that the kissnp output, in which fasta may be on several lines, cannot be read. Need the kissreads output .fa file "
    print "  -l log file (reference to find). Each polymorphism = one line: type[_id] position. e.g. \"SNP_2999 4637891\""
    print "    OR: "
    print "  -L simulator log file. While muting the genomes the simulator generates a log file containing the muted positions that are parsed and used"
    print "  -g gassst output file: file produced by gassst after mapping the disco output on the non-muted reference genome only with 100% identity"
    print "  -t threshold: conserves from the predicted SNPs, only thos with a rank value bigger or equal to this threshold. Default = 0"
    print "  -v verbose mode: outputs also the list of discoMore false positives and the list of false negative"
    print "  -r roc mode: ouputs files roc_SNP and roc_INDEL containing data used for generating roc curves"
    print "  -h: help"
    print "-----------------------------------------------------------------------------"
    sys.exit(2)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hrvd:L:l:g:t:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    generic=True
    threshold=0
    verbose=False
    roc=False
    discoFile=""
    logFile=""
    gassstFile=""
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-r"):
            roc = True
        elif opt in ("-v"):
            verbose = True
        elif opt in ("-t"):
            threshold = float(arg)
        elif opt in ("-d"):
            discoFile = arg
        elif opt in ("-l"):
            logFile = arg
        elif opt in ("-L"):
            logFile = arg
            generic=False
        elif opt in ("-g"):
            gassstFile = arg
        else:
            assert False, "unhandled option"


    if discoFile=="": 
         print "Disco file is missing"
         sys.exit(2)

    if logFile=="": 
         print "log file and reference_positions_file are missing, we need one of them."
         sys.exit(2)
         

    if gassstFile=="": 
         print "gassst file is missing"
         sys.exit(2)
    

    genotypes_predictions=get_genotypes(discoFile)
    
    nb_good_genotype_predictions=0
    nb_bad_genotype_predictions=0
    print "SNPs" , 
    doAll(discoFile, logFile, gassstFile, "SNP", threshold, verbose,generic,roc, nb_good_genotype_predictions, nb_bad_genotype_predictions, genotypes_predictions)
    
    

    nb_good_genotype_predictions=0
    nb_bad_genotype_predictions=0
    print "indels" , 
    doAll(discoFile, logFile, gassstFile, "INDEL", threshold, verbose,generic,roc, nb_good_genotype_predictions, nb_bad_genotype_predictions, genotypes_predictions)
    
    

if __name__ == "__main__":
     main()  
    
# doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "SNP")
#doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "INDEL")
