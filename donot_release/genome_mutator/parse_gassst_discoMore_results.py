#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt

def findNbMachedPolymorphism(list_simutated_polymorphism, start, stop):
    nbMapped=0 # True positive prediction 
    for pos in range(start, stop+1):
        if pos in list_simutated_polymorphism and list_simutated_polymorphism[pos]==False:
            list_simutated_polymorphism[pos]=True
            nbMapped+=1
    return nbMapped
    


# Compute Number of predicted polymorphism
def getPredictedPolymorphism (discoResultsFile, threshold, typePolymorphism):
    list_predicted_polymorphism={}
    filin = open(discoResultsFile, 'r') 
    while 1: 
        line=filin.readline() #>SNP_higher_path_95|high|nb_pol_1|C1_0|C2_40|rank_1.00000
        if not line: break
        if line.startswith(">"+typePolymorphism):
            pol_id=int(line.split("|")[0].split("_")[-1])
            rank=float(line.split("|")[-1].split("_")[-1])
            if rank>=threshold:
                list_predicted_polymorphism[pol_id] = False # This polymorphism has not been confirmed yet
            filin.readline() # read sequence upper, we don't care
            filin.readline() # read comment lower, we don't care
            filin.readline() # read sequence lower, we don't care  
    return list_predicted_polymorphism

# Store the polymorphism positions
def storeThePolymorphismPositions(simulatorLogFile, typePolymorphism):
    filin = open(simulatorLogFile, 'r')
    list_simutated_polymorphism={}
    while 1:
        line=filin.readline()
        if not line: break
        if line.startswith(">"+typePolymorphism):
            pos=int(line.split("|")[2]) #>SNP_2999|lower|4637891|C/T
            if pos in list_simutated_polymorphism:
                print "Warning, more than one "+typePolymorphism+" was simulated position "+str(pos)
            else:
                list_simutated_polymorphism[pos]=False # Will be set to true if this SNP is detected. 
        filin.readline() # read sequence upper, we don't care
        filin.readline() # read comment lower, we don't care
        filin.readline() # read sequence lower, we don't care
    return list_simutated_polymorphism




# Compute SNP precision
def checkMappedPaths(gassstOutFile,typePolymorphism, list_simutated_polymorphism, list_predicted_polymorphism):
    filin = open(gassstOutFile, 'r')
    for line in filin.readlines(): # SNP_lower_path_95|high|nb_pol_1|C1_41|C2_0|rank_1.00000 ...
        if line.startswith(typePolymorphism): 
            matching_position=int(line.split("\t")[3])
            predicted_pol_id=int(line.split("|")[0].split("_")[-1])
            span_match=int(line.split("\t")[4][:-1])
             # check if this mapping is in the list of simulated stuffs
            localNbMapped = findNbMachedPolymorphism(list_simutated_polymorphism,matching_position-1,matching_position+span_match-1)
            if localNbMapped>0:
                list_predicted_polymorphism[predicted_pol_id] = True # We mapped the polymorphism(s) corresponding to this predicted id. 


def print_results(list_predicted, list_simulated, polymorphism, threshold):
    nbTP=0
    for i in list_predicted:
        if list_predicted[i]:
            nbTP+=1
    nbSimulatedValidated=0

    for i in list_simulated:
        if list_simulated[i]:
            nbSimulatedValidated+=1
            
    print "-------------------------------"
    print "             ",polymorphism
    if threshold>0:
        print " from",polymorphism,"predicted with a rank bigger or equal to",threshold,":"
    print len(list_simulated), polymorphism,"were simultated.\t Among them", nbSimulatedValidated, "are correctly predicted"
    print len(list_predicted), polymorphism,"were predicted.\t Among them", nbTP, "are correctly mapped"
    if len(list_predicted)>0: print polymorphism,"precision\t %.2f"% float(nbTP/float(len(list_predicted))*100)
    if len(list_simulated)>0: print polymorphism,"recall\t %.2f"% float(nbSimulatedValidated/float(len(list_simulated))*100)
    print "-------------------------------"
    
def FN(list_simulated, polymorphism):
    print "-----------------------------------------------"
    print polymorphism,"FALSE NEGATIVE GENOME POSITIONS "
    for i in list_simulated:
        if not list_simulated[i]:
            print "FN",polymorphism, "position", i
    print "-----------------------------------------------"
    
    
def FP(list_predicted, polymorphism):
    print "-----------------------------------------------"
    print polymorphism,"FALSE POSITIVES PREDICTION IDs "
    for i in list_predicted:
        if not list_predicted[i]:
            print "FP",polymorphism, "id", i
    print "-----------------------------------------------"

def doAll(discoResFile, logFile, gassstResFile, polymorphism, threshold, verbose):
    
    list_predicted=getPredictedPolymorphism(discoResFile, threshold,  polymorphism)
    list_simutated=storeThePolymorphismPositions(logFile, polymorphism)
    checkMappedPaths(gassstResFile,polymorphism, list_simutated, list_predicted)
    print_results(list_predicted, list_simutated, polymorphism, threshold)
    if verbose:
        FP(list_predicted,polymorphism)
        FN(list_simutated,polymorphism)
    

def usage():
    '''Usage'''
    print "-----------------------------------------------------------------------------"
    print sys.argv[0]," : parse the gassst result on discoMore files"
    print "-----------------------------------------------------------------------------"
    print "usage: ",sys.argv[0]," -d discoMore coherent_out_file -l simulator_log_file -g gassst_out_file [-t float] [-v] "
    print "  -d discoSnp outfile. Note that the kissnp output, in which fasta may be on several lines, cannot be read. Need the kissreads output .fa file "
    print "  -l simulator log file. While muting the genomes the simulator generates a log file containing the muted positions that are parsed and used"
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
    

    for type_polymorphism in {"SNP","INDEL"}:
         doAll(discoFile, logFile, gassstFile, type_polymorphism, threshold, verbose)

if __name__ == "__main__":
     main()  
    
# doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "SNP")
#doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "INDEL")