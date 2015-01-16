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
    while 1:
        line=filin.readline() # SNP_2999 4637891 2
        if not line: break
        if line.startswith(typePolymorphism):
            pos=int(line.split(" ")[1]) 
            chromosome=line.split(" ")[2].strip()
            if not chromosome in list_reference_polymorphism:                
                list_reference_polymorphism[chromosome]={}
            if pos in list_reference_polymorphism[chromosome]:
                print "Warning, more than one "+typePolymorphism+" exists position "+str(pos)+" on chromosome "+chromosome
            else:
                list_reference_polymorphism[chromosome][pos]=False # Will be set to true if this SNP is detected. 

    return list_reference_polymorphism


####################################################################
#                         REFERENCES
####################################################################
# Store the polymorphism positions
def storeThePolymorphismPositionsFromLogFile(simulatorLogFile, typePolymorphism): 
    filin = open(simulatorLogFile, 'r')
    list_simulated_polymorphism={} # for each chromosome: contains a list of positions
    while 1:
        line=filin.readline()
        if not line: break
        if line.startswith(">"+typePolymorphism):
            pos=int(line.split("|")[2]) #>SNP_2999|lower|4637891|chr1|C/T
            chromosome=line.split("|")[3]
            if not chromosome in list_simulated_polymorphism: 
                list_simulated_polymorphism[chromosome]={}
            if pos in list_simulated_polymorphism[chromosome]:
                print "Warning, more than one "+typePolymorphism+" was simulated position "+str(pos)+" on chromosome "+chromosome
            else:
                list_simulated_polymorphism[chromosome][pos]=False # Will be set to true if this SNP is detected. 
                
        filin.readline() # read sequence upper, we don't care
        filin.readline() # read comment lower, we don't care
        filin.readline() # read sequence lower, we don't care
    return list_simulated_polymorphism



 

# Compare mapped positions with reference positions
def checkMappedDiscoPaths(gassstOutFile,typePolymorphism, list_reference_polymorphism, list_predicted_polymorphism, threshold, list_predicted_positions):
    filin = open(gassstOutFile, 'r')
    nbindel=0 #DEB
    for line in filin.readlines(): # SNP_higher_path_1390|high|nb_pol_1|C1_629|C2_12|Q1_54|Q2_69|rank_0.98621        12      1       200256
                                # or INDEL_9_higher_path_995|high|nb_pol_1|C1_32|C2_0|rank_1.00000	gi|545778205|gb|U00096.3|
        if line.startswith(typePolymorphism): 
            this_threshold=1
            if threshold>0:
                this_threshold=float(line.split("\t")[0].split("|")[-1].split("_")[-1])
            if this_threshold<threshold: 
                continue
            matching_chromosome=line.split("\t")[1].split("|")[0].strip()
            matching_position=int(line.split("\t")[3])
            good_prefix=line.split("|")[0] # if this is a discoSnp++ output, this is correct
            if len(good_prefix.split("\t")[0]) < len (good_prefix): # this is a cortex output: "SNP_higher_path_1	gi|545778205|gb|U00096.3|	0	3595700	63M"
                good_prefix=good_prefix.split("\t")[0]
            predicted_pol_id=int(good_prefix.split("_")[-1])
            # avoid to map twice the same variant (both upper and lower path as the mapping is not exact)
            if predicted_pol_id in list_predicted_polymorphism and list_predicted_polymorphism[predicted_pol_id] == True: continue
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
                    nbindel+=1
                    localNbMapped = parse_gassst_common.findNbMachedPolymorphism(
                        list_reference_polymorphism[matching_chromosome],
                        start,
                        stop,
                        [i for i in range(stop-start+1)], # fake acceptation of all mapped indel.
                        typePolymorphism, 
                        forward) 
                else:
                    localNbMapped = parse_gassst_common.findNbMachedPolymorphism(
                    list_reference_polymorphism[matching_chromosome],
                    start,
                    stop,
                    list_predicted_positions[predicted_pol_id], 
                    typePolymorphism, 
                    forward)
            if localNbMapped>0:
                list_predicted_polymorphism[predicted_pol_id] = True # We mapped the polymorphism(s) corresponding to this predicted id. 
                

    print "nbindel",nbindel


####################################################################
#                         MAIN
####################################################################
def doAll(discoResFile, logFile, gassstResFile, polymorphism, threshold, verbose,generic):
    
    res=parse_gassst_common.getPredictedPolymorphism(discoResFile, threshold,  polymorphism)
    list_predicted=res[0]
    nb_predicted=res[1]
    list_predicted_positions=res[2]
    if generic:
        list_reference=storeThePolymorphismPositionsGenericList(logFile, polymorphism)
    else:
        list_reference=storeThePolymorphismPositionsFromLogFile(logFile, polymorphism)
    checkMappedDiscoPaths(gassstResFile,polymorphism, list_reference, list_predicted,threshold,list_predicted_positions)
    parse_gassst_common.print_results(nb_predicted, list_reference, polymorphism, threshold)
    if verbose:
        parse_gassst_common.FP(list_predicted,polymorphism)
        parse_gassst_common.FN(list_reference,polymorphism)
    

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
    print "usage: ",sys.argv[0]," -d discoMore coherent_out_file -l log_file -g gassst_out_file [-t float] [-v] "
    print "  -d discoSnp outfile. Note that the kissnp output, in which fasta may be on several lines, cannot be read. Need the kissreads output .fa file "
    print "  -l log file (reference to find). Each polymorphism = one line: type[_id] position. e.g. \"SNP_2999 4637891\""
    print "    OR: "
    print "  -L simulator log file. While muting the genomes the simulator generates a log file containing the muted positions that are parsed and used"
    print "  -g gassst output file: file produced by gassst after mapping the disco output on the non-muted reference genome only with 100% identity"
    print "  -t threshold: conserves from the predicted SNPs, only thos with a rank value bigger or equal to this threshold. Default = 0"
    print "  -v verbose mode: outputs also the list of discoMore false positives and the list of false negative"
    print "  -h: help"
    print "-----------------------------------------------------------------------------"
    sys.exit(2)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvd:L:l:g:t:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    generic=True
    threshold=0
    verbose=False
    discoFile=""
    logFile=""
    gassstFile=""
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
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
    

    for type_polymorphism in {"SNP","INDEL"}:
         doAll(discoFile, logFile, gassstFile, type_polymorphism, threshold, verbose,generic)

if __name__ == "__main__":
     main()  
    
# doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "SNP")
#doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "INDEL")