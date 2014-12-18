#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt


####################################################################
#                         PREDICTIONS
####################################################################
# Compute Number of predicted polymorphism
def getPredictedPolymorphism (discoResultsFile, threshold, typePolymorphism):
    list_predicted_polymorphism={}
    filin = open(discoResultsFile, 'r') 
    nb_pol=0
    while 1: 
        line=filin.readline() #>SNP_higher_path_95|high|nb_pol_1|C1_0|C2_40|rank_1.00000 OR >INDEL_9_higher_path_891|high|nb_pol_1|C1_0|C2_30|rank_1.00000
        if not line: break
        if line.startswith(">"+typePolymorphism):
            pol_id=int(line.split("|")[0].split("_")[-1])
            rank=float(line.split("|")[-1].split("_")[-1])                        
            if rank>=threshold:
                nb_pol+=int(line.split("|")[2].split("_")[-1])
                list_predicted_polymorphism[pol_id] = False # This polymorphism has not been confirmed yet
            filin.readline() # read sequence upper, we don't care
            filin.readline() # read comment lower, we don't care
            filin.readline() # read sequence lower, we don't care
    return list_predicted_polymorphism,nb_pol
    
   

####################################################################
#                         REFERENCES
####################################################################
# Store the polymorphism positions
def storeThePolymorphismPositions(referencePositionFile, typePolymorphism):
    filin = open(referencePositionFile, 'r')
    list_reference_polymorphism={} # for each chromosome: contains a list of positions
    while 1:
        line=filin.readline() # SNP_2999 4637891 2
        if not line: break
        if line.startswith(typePolymorphism):
            pos=int(line.split(" ")[1]) 
            chromosome=line.split(" ")[2].strip()
            if chromosome in list_reference_polymorphism:                
                if pos in list_reference_polymorphism[chromosome]:
                    print "Warning, more than one "+typePolymorphism+" exists position "+str(pos)
                else:
                    list_reference_polymorphism[chromosome][pos]=False # Will be set to true if this SNP is detected. 
            else:
                list_reference_polymorphism[chromosome]={}
                list_reference_polymorphism[chromosome][pos]=False
    return list_reference_polymorphism


####################################################################
#                         MAPPING
####################################################################
# Stores polymorphisms whose both paths map the reference genome.
def checkDoubleMappedPaths(gassstOutFile,typePolymorphism):
    predicted_polymorphism_lower_mapped={}
    predicted_polymorphism_upper_mapped={}
    predicted_polymorphism_both_mapped={}
    filin = open(gassstOutFile, 'r')
    for line in filin.readlines(): # SNP_higher_path_1390|high|nb_pol_1|C1_629|C2_12|Q1_54|Q2_69|rank_0.98621        12      1       200256
        if line.startswith(typePolymorphism): 
            predicted_pol_id=int(line.split("|")[0].split("_")[-1])
            upper=True
            if line.split("|")[0].split("_")[1]=="lower":
                upper=False
            if upper:
                predicted_polymorphism_upper_mapped[predicted_pol_id] = True
            else:
                predicted_polymorphism_lower_mapped[predicted_pol_id] = True
            if predicted_pol_id in predicted_polymorphism_upper_mapped and predicted_pol_id in predicted_polymorphism_lower_mapped:
                predicted_polymorphism_both_mapped[predicted_pol_id] = True
    return predicted_polymorphism_both_mapped
    
    
def findNbMachedPolymorphism(list_reference_polymorphism, start, stop):
    nbMapped=0 # True positive prediction 
    for pos in range(start, stop+1):
        if pos in list_reference_polymorphism and list_reference_polymorphism[pos]==False:
            list_reference_polymorphism[pos]=True
            nbMapped+=1
            # uncomment if you want to see where FP are.
    # if nbMapped==0:
  #       print start
    return nbMapped
    

# Compare mapped positions with reference positions
def checkMappedPaths(gassstOutFile,typePolymorphism, list_reference_polymorphism, list_predicted_polymorphism,predicted_polymorphism_both_mapped, threshold):
    filin = open(gassstOutFile, 'r')
    for line in filin.readlines(): # SNP_higher_path_1390|high|nb_pol_1|C1_629|C2_12|Q1_54|Q2_69|rank_0.98621        12      1       200256
        if line.startswith(typePolymorphism): 
            this_threshold=1
            if threshold>0:
                this_threshold=float(line.split("|")[7].split("\t")[0].split("_")[-1])
            if this_threshold<threshold: 
                continue
            matching_chromosome=line.split("\t")[1].strip()
            matching_position=int(line.split("\t")[3])
            predicted_pol_id=int(line.split("|")[0].split("_")[-1])
            span_match=int(line.split("\t")[4][:-1])
             # check if this mapping is in the list of reference stuffs
            localNbMapped=0
            if not matching_chromosome in list_reference_polymorphism: continue
            if predicted_pol_id in list_predicted_polymorphism and not predicted_pol_id in predicted_polymorphism_both_mapped:
                localNbMapped = findNbMachedPolymorphism(list_reference_polymorphism[matching_chromosome],matching_position-1,matching_position+span_match-1)
            if localNbMapped>0:
                list_predicted_polymorphism[predicted_pol_id] = True # We mapped the polymorphism(s) corresponding to this predicted id. 


####################################################################
#                         IOs
####################################################################
def print_results(nb_predicted, list_reference, polymorphism, threshold):
    nbTP=0
    for chro in list_reference:
        for i in list_reference[chro]:
            if list_reference[chro][i]:
                nbTP+=1
            
    print "-------------------------------"
    print "             ",polymorphism
    if threshold>0:
        print " from",polymorphism,"predicted with a rank bigger or equal to",threshold,":"
    print len(list_reference), polymorphism,"in the reference.\t Among them", nbTP, "are correctly predicted"
    print nb_predicted, polymorphism,"were predicted.\t Among them", nbTP, "are correctly mapped"
    if nb_predicted>0: print polymorphism,"precision\t %.2f"% float(nbTP/float(nb_predicted)*100)
    if len(list_reference)>0: print polymorphism,"recall\t %.2f"% float(nbTP/float(len(list_reference))*100)
    print "-------------------------------"
    
def FN(list_reference, polymorphism):
    print "-----------------------------------------------"
    print polymorphism,"FALSE NEGATIVE GENOME POSITIONS "
    for chro in list_reference:
        for i in list_reference[chro]:
            if not list_reference[chro][i]:
                print "FN",polymorphism, "chr",chro," position", i
    print "-----------------------------------------------"
    
    
def FP(list_predicted, polymorphism):
    print "-----------------------------------------------"
    print polymorphism,"FALSE POSITIVES PREDICTION IDs "
    for i in list_predicted:
        if not list_predicted[i]:
            print "FP",polymorphism, "id", i
    print "-----------------------------------------------"


####################################################################
#                         MAIN
####################################################################
def doAll(discoResFile, logFile, gassstResFile, polymorphism, threshold, verbose):
    
    res=getPredictedPolymorphism(discoResFile, threshold,  polymorphism)
    list_predicted=res[0]
    nb_predicted=res[1]
    list_reference=storeThePolymorphismPositions(logFile, polymorphism)
    predicted_polymorphism_both_mapped = checkDoubleMappedPaths(gassstResFile,polymorphism)
    checkMappedPaths(gassstResFile,polymorphism, list_reference, list_predicted,predicted_polymorphism_both_mapped,threshold)
    print_results(nb_predicted, list_reference, polymorphism, threshold)
    if verbose:
        FP(list_predicted,polymorphism)
        FN(list_reference,polymorphism)
    

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
    print "  -g gassst output file: file produced by gassst after mapping the disco output on the non-muted reference genome only with 100% identity"
    print "  -t threshold: conserves from the predicted SNPs, only thos with a rank value bigger or equal to this threshold. Default = 0"
    print "  -v verbose mode: outputs also the list of discoMore false positives and the list of false negative"
    print "  -h: help"
    print "-----------------------------------------------------------------------------"
    sys.exit(2)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvd:l:g:t:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
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
        elif opt in ("-g"):
            gassstFile = arg
        else:
            assert False, "unhandled option"


    if discoFile=="": 
         print "Disco file is missing"
         sys.exit(2)

    if logFile=="": 
         print "log file is missing"
         sys.exit(2)
         

    if gassstFile=="": 
         print "gassst file is missing"
         sys.exit(2)
    

    for type_polymorphism in {"SNP", "INDEL"}:
         doAll(discoFile, logFile, gassstFile, type_polymorphism, threshold, verbose)

if __name__ == "__main__":
     main()  
    
# doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "SNP")
#doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "INDEL")