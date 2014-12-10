#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt

def find_mached_polymorphism(list_simutated_polymorphism, start, stop):
    nbMapped=0 # True positive prediction 
    for pos in range(start, stop+1):
        if pos in list_simutated_polymorphism and list_simutated_polymorphism[pos]==False:
            list_simutated_polymorphism[pos]=True
            nbMapped+=1
    return nbMapped
    


# Compute Number of predicted SNPs
def getNumberPredictedPolymorphism (discoResultsFile, typePolymorphism):
    nbPredicted=0
    filin = open(discoResultsFile, 'r') # TODO: user defined
    while 1:
        line=filin.readline()
        if not line: break
        if line.startswith(">"+typePolymorphism):
            nbPredicted+=1
            filin.readline() # read sequence upper, we don't care
            filin.readline() # read comment lower, we don't care
            filin.readline() # read sequence lower, we don't care    
    return nbPredicted

# Store the polymorphism positions
def storeThePolymorphismPositions(simulatorLogFile, typePolymorphism):
    filin = open(simulatorLogFile, 'r')
    list_simutated_polymorphism={}
    nbSimulated=0
    while 1:
        line=filin.readline()
        if not line: break
        if line.startswith(">"+typePolymorphism):
            pos=int(line.split("|")[2]) #>SNP_2999|lower|4637891|C/T
            if pos in list_simutated_polymorphism:
                print "Warning, more than one "+typePolymorphism+" was simulated position "+str(pos)
            else:
                list_simutated_polymorphism[pos]=False # Will be set to true if this SNP is detected. 
                nbSimulated+=1
        filin.readline() # read sequence upper, we don't care
        filin.readline() # read comment lower, we don't care
        filin.readline() # read sequence lower, we don't care
    return nbSimulated,list_simutated_polymorphism




# Compute SNP precision
def getNumberMapped(gassstOutFile,typePolymorphism, list_simutated_polymorphism):
    nbMapped=0
    filin = open(gassstOutFile, 'r')
    for line in filin.readlines():
        if line.startswith(typePolymorphism):
            matching_position=int(line.split("\t")[3])
            span_match=int(line.split("\t")[4][:-1])
            nbMapped+=find_mached_polymorphism(list_simutated_polymorphism,matching_position-1,matching_position+span_match-1)
    return nbMapped


def print_results(simulated, predicted, mapped, polymorphism):
    print "-------------------------------"
    print "             ",polymorphism
    print simulated, polymorphism,"were simultated"
    print predicted, polymorphism,"were predicted"
    print mapped, polymorphism,"were correctly mapped"
    if predicted>0: print polymorphism,"precision %.2f"% float(mapped/float(predicted)*100)
    if simulated>0: print polymorphism,"recall %.2f"% float(mapped/float(simulated)*100)
    print "-------------------------------"
    


def doAll(discoResFile, logFile, gassstResFile, polymorphism):
    nbPredictedSnps=getNumberPredictedPolymorphism(discoResFile,polymorphism)
    res=storeThePolymorphismPositions(logFile, polymorphism)
    nbSnpSimulated=res[0]
    list_simutated_snps=res[1]
    nbSnpMapped=getNumberMapped(gassstResFile,polymorphism, list_simutated_snps)
    print_results(nbSnpSimulated, nbPredictedSnps, nbSnpMapped, polymorphism)
    

def usage():
    '''Usage'''
    print "-----------------------------------------------------------------------------"
    print sys.argv[0]," : parse the gassst result on discoMore files"
    print "-----------------------------------------------------------------------------"
    print "usage: ",sys.argv[0]," -d discoMore coherent_out_file -l simulator_log_file -g gassst_out_file"
    print "  -h: help"
    print "-----------------------------------------------------------------------------"
    sys.exit(2)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:l:g:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    

    discoFile=""
    logFile=""
    gassstFile=""
    for opt, arg in opts:
         if opt in ("-h", "--help"):
             usage()
             sys.exit()
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
         doAll(discoFile, logFile, gassstFile, type_polymorphism)

if __name__ == "__main__":
     main()  
    
# doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "SNP")
#doAll('test_k_31_c_5_D_10_b_1_coherent.fa', 'coli_formated.fasta_log', "afac", "INDEL")