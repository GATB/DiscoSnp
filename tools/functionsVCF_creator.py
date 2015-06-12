#!/bin/python
# -*- coding: utf-8 -*-
###############################################
#v10022015Disco++ : Check des résultats - "compilation" séparée - voir pour mieux trier les fonctions
import os
import sys
import subprocess
import re
import time

from sets import Set # Added by pierre (12/06/2015)
#FUNCTIONS_________________________________________________________________________________________________________
#INDEX_____________________________________________________________________________________________________________
#       Counting : """Function that counts the number of SNPs and the number of genotype"""
#       ParsingDiscoSNP
#       GetCouple :"""Retrieves for each path alignment information in a list ; retrieves a dictionary with all the positions of a path and the number of associated mismatch ; retrieves a Boolean indicating whether the SNP is multimapped or not"""
#       GetCoverage :"""Takes as input lists list  """
#       AddPosition :"""Add a position to a set of positions only if the difference between the two is greater than delta"""
#       ValidationSNP :"""Main function of the snp validation : check if the couple is validated with only one mapping position at the number of mismatch tested"""
#       ValidationSNPBestHits :"""Main function of the snp validation : check if the couple is validated with only one mapping position """
#       ReverseComplement 
#       CigarCodeChecker :"""Function which allows to recover the position of the SNPs according to the position of the deletions or insertions relative to the reference sequence (CigarCode Parsing : checks for insertion deletion or soft clipping"""
#       ReferenceChecker : """Function which allows to get the MD tag parsing; checks if path nucleotide is identical to the reference nucleotide"""
#       GetSequence :"""Gets and reverses if it is needed sequence"""
#       RecupPosSNP :"""Determines : if the paths are mapped in reverse or forward sense ; which paths have the variant identical to the genome (information from the MD Tag) ; the position of mapping of the variant for each path by taking into account the shift with the genome ( soft clipping insertion, deletion  and so on : information from the cigar code) ; the nucleotide variant for the upper path and the lower path"""
#       MismatchChecker :"""In case of divergent main position (case snp whose two paths are mapped ) to define the reference = > check the number of mismatch( If the number of mismatch is the same in both cases it is the lower lexicographical SNP which is selected for reference .The Boolean allows to know the reference SNP ) """
#       GetPolymorphism :'''Gets from the dicoHeader all the positions, and the nucleotides (R means that it's on the reverse strand)SNP :  one variant correspond to listPos[0], listPosR[0], listnucleoUp[0], listnucleoLow[0],listnucleoUpR[0],listnucleoLowR[0]''' 
#       fillVCFSimpleSnp :""" Fills the different fields of vcf based on boolean ( on whether the SNP is identical or not the reference) , if neither is identical to the reference = > we take the SNP comes first in the lexicographical order"""
#       printOneline : """Prints the line of the current SNP in the VCF file."""
#       printVCFSNPclose : """Prints the line of the close SNPs in the VCF file. The first variant will determine which path will be used as reference (checks with the boolean boolRef which path is identical to the reference) """ 
#       GetGenotype : """Gets the genotype, the coverage and the likelihood by sample and print it in the correspondand fields. The genotype is determined by DiscoSnp++ (which considered the upper path as reference). If the “REF” corresponds the upper path, the genotype in the VCF is identical to the genotype in DiscoSnp++, else  it's the opposite ( 1/1 becomes 0/0 and so on)."""
#       PrintVCFGhost :'''Without samfile some fields will have default value : CHROM table[0],POS table[1], QUAL table[5], FILTER [6]''' 
#       FillVCF:"""Take all necessary input variables to fill the vcf;  Fills the fields of the table which will be printed in the vcf ; return the table"""
#       ReverseSeq :"""Take a sequence and reverse it""" 
#       CheckBitwiseFlag :"""Checks if the BitwiseFlag contains the tested value such as : read reverse strand, read unmmaped and so on."""
#       CheckContigUnitig : """Checks if there is an extension in the form of contig or unitig to take into account in the mapping position"""
##############################################################
##############################################################
def Counting(fichier):
    """Function that counts the number of SNPs and the number of genotype"""
    fileName, fileExtension = os.path.splitext(fichier)
    samfile=open(fichier,'r')
    nbSnp=0
    nbGeno=0
    nbclose=0
    nbsimplesnp=0
    nbINDEL=0
    if ".sam"==fileExtension:
            for line in samfile :
                matchInt=re.match(r'^@',line)
                if matchInt:#Removes the header lines from the loop
                    continue
                else:
                    nbSnp+=1
                    if isinstance(line,str): #Test the samline ti know if we have to parse it
                        line=line.rstrip('\n')
                        line=line.split('\t')#Splits the line into a list
                    for i in line:
                        if 'SNP' or 'INDEL' in i:
                            nomDisco=i.split('|')
                            nb_pol=nomDisco[3].split("_")
                            nbSnp+=int(nb_pol[2])# Takes into account the number of variant => close snps
                            break
                    if nbGeno==0: #Count the number of genotypes in the header
                        for k in nomDisco:
                            match=re.match(r'^G',k)
                            if match:
                                nbGeno+=1
    else:
         for line in samfile :
                if ">" in line:
                        listLine=line.split("|")
                        if nbGeno==0:#We just count once the number of genotype
                                for nom in listLine:
                                    match=re.match(r'^G\d+',nom)
                                    if match:
                                        nbGeno+=1   
                        if "SNP_" in line:
                                listLine=line.split("|")
                                listnbpol=listLine[3].split("_")
                                if int(listnbpol[2])>1:
                                        nbclose+=int(listnbpol[2])
                                if int(listnbpol[2])==1:
                                        nbsimplesnp+=1
                        elif "INDEL_" in line:
                                nbINDEL+=1
         nbSnp=nbsimplesnp+nbclose+nbINDEL
    samfile.close()
    return(int(nbSnp)/2,nbGeno)

##############################################################
#Take the new discoSnp++ header : example SNP_higher_path_13736|P_1:30_A/G|high|nb_pol_1|C1_53425|C2_1|C3_30|C4_3|C5_19|G1_0/0:1947,160066,1067676|G2_0/1:40,9,20|G3_1/1:1895477,284398,3575|G4_0/1:34,9,54|G5_1/1:1306661,196101,2493|Q1_66|Q2_61|Q3_36|Q4_55|Q5_33|rank_0.99909
#INDEL_lower_path_2355|P_1:30_1_2|high|nb_pol_1|C1_56|C2_23|C3_1721540|C4_41|C5_89|G1_0/1:776,17,955|G2_0/1:318,16,457|G3_1/1:34427840,5179573,73399|G4_0/1:570,16,710|G5_0/1:1510,84,313|Q1_66|Q2_66|Q3_69|Q4_63|Q5_66|rank_0.60362
##############################################################
def ParsingDiscoSNP(snp,boolNum):
    if isinstance(snp,str):#Converts the line of the samfile into a list
        snp=snp.rstrip('\r')
        snp=snp.rstrip('\n')
        snp=snp.split('\t')
    listCoverage=[]
    listC=[]
    j=0
    numSNP=None
    unitigLeft=None
    unitigRight=None
    contigLeft=None
    contigRight=None
    rank=None
    i=0
    nb_pol=0
    valRank=None
    nomDisco=None
    pos=None
    geno=[]
    posD=0
    ntUp=''
    ntLow=''
    dicoHeader={}
    dicoGeno={}
    SNP=0
    INDEL=False
    ID=None
    matchInt=None
    key=None
    chaine=None
    ind=None
    chaine1=None
    chaine2=None
    amb=None
    rank=None
    coverage=None
    nb_pol=None
    gen=None
    listegeno=None
    listlikelihood=None
    for i in snp:
        if 'SNP' or 'INDEL' in i:
            nomDisco=i.split('|')
            break
    for i in nomDisco:# Loops on the header of discoSnp++
        if "SNP" in i: #Important to know if we have a snp for the creation of the dictionary with the information with nucleotide and position
            discoName=i
            SNP=1 
            ID= i.split('_') 
            numSNP=ID[3] #Gets the ID of the snp
        elif 'P_' in i: ## Esssential to have the position on the path of the allele, and the nucleotide of the upper and the lower path (forward direction) => information given by DiscoSnp++
            pos=i.split(',')#Parsing if there is multiple snps on the paths
            ntUp=''
            ntLow=''
            key=''
            chaine=[]
            ind=0
            dicoHeader={}
            for j in pos:#Loops on the list if there are close snps
                posD=0
                if SNP==1: #Case of snps
                    if ":" in j:#P_1:30_A/G or P_1:30_A/G,P_2:31_G/A
                        chaine=j.split(":")
                        key=chaine[0]
                        chaine1=chaine[1].split('_')
                        posD=chaine1[0]
                        chaine2=chaine1[1].split('/')
                        ntUp=chaine2[0]
                        ntLow=chaine2[1]
                        dicoHeader[key]=[posD,ntUp,ntLow] ##### !!! Essential in case of snp
                else:#Case of indel 
                    if ":" in j:# P_1:30_3_2
                        ind=0
                        posD=0
                        chaine=j.split(":")
                        amb=0
                        key=chaine[0]
                        chaine1=chaine[1].split('_')
                        posD=int(chaine1[0])
                        ind=int(chaine1[1])
                        amb=int(chaine1[2])
                        dicoHeader[key]=[posD,ind,amb] ##### !!! Essential in case of indel 
        elif "INDEL" in i: #Gets the ID of an indel
            INDEL=True
            ID= i.split('_')
            discoName=i
            numSNP=ID[3]
        elif "unitig" in i: #Checks the presence and the size of unitig 
            if "left" in i:
                unitig=i.split('_')
                unitigLeft=unitig[3]
            elif "right" in i :
                unitig=i.split('_')
                unitigRight=unitig[3]
        elif "contig" in i:#Checks the presence and the size of contig
            if "left" in i:
                contig=i.split('_')
                contigLeft=contig[3]
            elif "right" in i :
                contig=i.split('_')
                contigRight=contig[3]
        elif "rank" in i: #Gets the rank value
            rank=i.split('_')
            valRank=rank[1]
        elif "C" in i: ###Gets the coverage for the path
            coverage=i.split('_')
            j=0
            for j in range(len(coverage)):
                matchInt=re.match(r'^\d+$',coverage[j])
                if matchInt:
                    listCoverage.append(str(coverage[j]))
                    listC.append(str(coverage[j-1]))
        elif "nb_pol" in i : ###Gets the number of variant for the paths
            nb_pol=i.split('_')
            nb_pol=nb_pol[2]
        elif "G" in i: ##Gets the genotype and likelihood by samples
            gen=i.replace("_",":")
            listgeno=gen.split(":")
            if len(listgeno)>2:
                listlikelihood=listgeno[2].split(",")
            else:
                listlikelihood=0
            dicoGeno[listgeno[0]]=[listgeno[1],listlikelihood] ##dictionary with the genotype by sample and a list with the likelihood
    
    if boolNum==0:
        return(discoName,snp,numSNP,unitigLeft,unitigRight,contigLeft,contigRight,valRank,listCoverage,listC,nb_pol,posD,ntUp,ntLow,dicoGeno,dicoHeader)
    else:
        return(numSNP)

##############################################################
# snpUp: sam line 1
# snpLow: sam line 2
##############################################################
def GetCouple(snpUp,snpLow):
    """Retrieves for each path alignment information in a list ; retrieves a dictionary with all the positions of a path and the number of associated mismatch ; retrieves a Boolean indicating whether the SNP is multimapped or not"""
    posUp = {} #dictionary of all the positions (keys) associated with the number of mismatches (values) 
    posLow = {}
    i=0
    j=0
    #Boolean snps multimapped
    boolXAUp=False
    boolXALow=False
    listXA=None
    strXA=None
    position=None
    garbage=None
    nbMismatchUp=None
    nbMismatchLow=None
    #Error list with mapping positions very close to the firt position give by bwa
    listerreurUp=set([(int(snpUp[3])-1),(int(snpUp[3])+1),(int(snpUp[3])+2),(int(snpUp[3])+3),(int(snpUp[3])-3),(int(snpUp[3])-2),int(snpUp[3])])
    listerreurLow=set([(int(snpLow[3])-1),(int(snpLow[3])+1),(int(snpLow[3])+2),(int(snpLow[3])+3),(int(snpLow[3])-3),(int(snpLow[3])-2),int(snpLow[3])])
    #Creation of a dict with mapping position associate with number of mismatch
    if 'XA:Z' in ''.join(snpUp): # XA: tag for multiple mapping : Check if the upper path is multiple mapped
        i=0
        #Parsing XA tag
        listXA=snpUp[19].split(';')
        strXA = ','.join(listXA)
        listXA = strXA.split(',')
        listXA.pop()
        position=listXA[1:] #position=[pos1,cigarcode1,number of mismatch1 , pos2,cigarcode2,number of mismatch2,...]
        while i<len(position): #Runs through the list 4 by 4 to get all the positions 
            if abs(int(position[i])) not in listerreurUp : #Checks if the position is not too close to the main one
                boolXAUp=True
                posUp[abs(int(position[i]))]=int(position[i+2]) #the position is associated to the number of mismatch in a dictionary
            i+=4
    if 'XA:Z' in ''.join(snpLow):#XA: tag for multiple mapping : Checks if the lower path is multiple mapped 
        i=0
        listXA=snpLow[19].split(';')
        strXA = ','.join(listXA)
        listXA = strXA.split(',')
        listXA.pop()
        position=listXA[1:]
        while i<len(position):
            if abs(int(position[i])) not in listerreurLow :
                boolXALow=True
                posLow[abs(int(position[i]))]=int(position[i+2])
            i+=4
    #Adds first position give by BWA to the dictionary if the path is mapped
    if abs(int(snpUp[3]))>0:
        posMutUp,nbMismatchUp=GetTag(snpUp)
        posUp[abs(int(snpUp[3]))]=int(nbMismatchUp)
    if abs(int(snpLow[3]))>0:
        posMutLow,nbMismatchLow=GetTag(snpLow)    
        posLow[abs(int(snpLow[3]))]=int(nbMismatchLow)
    return(posUp,posLow,boolXAUp,boolXALow)
    ##posUp/posLow : dictionary with all the position associated with their mismatch number
    ##boolXAUp/boolXALow : boolean if TRUE : the path is multiple mapped for BWA  and the max mapping distance

##############################################################
#listCup=["C1","C2"]
#listCoverageUp=[23,45]
##############################################################
def GetCoverage( listCUp, listCLow, listCoverageUp, listCoverageLow):
    """Takes as input lists list  """
    i=0
    covUp=''
    listCovGeno=[]
    covLow=''
    ##Creates a string if the number of coverage is superior to 1
    if len( listCoverageUp)>1 and len( listCoverageLow)>1:
        while i<len(listCoverageUp): #Goes through the list of coverage for the upper path => Usefull if the reference path is the upper path
            listCovGeno.append(int( listCoverageUp[i])+int( listCoverageLow[i]))#Sums the two coverages
            if covUp!='': #Adds the following coverage to the string (upper coverage in the first position ; lower in second)
                covUp=str(covUp)+';'+str( listCUp[i])+'='+str( listCoverageUp[i])+","+str( listCoverageLow[i]) # creates a string with the coverage
            else: #Inits with the first coverages
                covUp=str( listCUp[i])+'='+str( listCoverageUp[i])+","+str( listCoverageLow[i])
            i+=1
        i=0
        covLow=''
        while i<len( listCoverageLow): #Goes through the list of coverage for the lower path => Usefull if the reference path is the lower path
            if covLow!='':#Adds the following coverage to the string (lower coverage in the first position ; upper in second)
                covLow=str(covLow)+';'+str( listCLow[i])+'='+str( listCoverageLow[i])+","+str( listCoverageUp[i])
            else:#Inits with the first coverages
                covLow=str( listCLow[i])+'='+str( listCoverageLow[i])+","+str( listCoverageUp[i])
            i+=1
    else: #In case of one dataset
        listCovGeno.append(int(listCoverageUp[0])+int(listCoverageLow[0]))
        covUp=str( listCUp[0])+'='+str( listCoverageUp[0])+","+str( listCoverageLow[0])
        covLow=str( listCLow[0])+'='+str( listCoverageLow[0])+","+str( listCoverageUp[0])
    return(covUp,covLow,listCovGeno) #String covUp "C1=5,23;C2=35,1" listCovGeno=[28,36]
##############################################################
#position : current position to add at the ensemble
#delta : minimum number of difference allowed between two positions
##############################################################
def AddPosition(position,ensemble,delta):
    """Add a position to a set of positions only if the difference between the two is greater than delta"""
    if len(ensemble)==0: #Case : it's the first position add to the set
        ensemble.append(position)
        return(ensemble)
    if abs(position-ensemble[0])>delta: #Checks in the postions to test have a difference greatest than delta with the first position of the set
        ensemble.append(position)
        return(ensemble)
    return(ensemble)
    #ensemble : set of positions for a path
##############################################################
# snpUp/snpLow : sam line
# posUp/posLow : dictionary with all the position associated with their mismatch number example : posUp[5687884684]=2
# NM : number of mismatch to test : in the main test every distance 0 to distance max 
##############################################################
def ValidationSNP(snpLow,posLow,snpUp,posUp,NM):
    """Main function of the snp validation : check if the couple is validated with only one mapping position at the number of mismatch tested"""
    couple=None
    listPos= []
    ensemble=None
    delta=3
    position=None
    nbMismatch=None
    if abs(int(snpUp[3]))<=0 and abs(int(snpLow[3]))<=0:
        couple="unmapped"
        return(couple)
    else:
        for position,nbMismatch in posUp.items(): # Checks for the upper path
            if nbMismatch==int(NM): #Checks if the distance with the reference for the tested position is identical to the tested distance 
                listPos=AddPosition(position,listPos,delta) #If the position is mapped at NM (number of mismatch tested) test if its not too close with an other position
                ensemble=set(listPos)
                if abs(len(ensemble))>1: #Case of many position of mapping for NM => multiple mapped  Exit function
                    couple="multiple"
                    return(couple)
        for position,nbMismatch in posLow.items():#Checks for the lower path
            if nbMismatch==int(NM):
                listPos=AddPosition(position,listPos,delta)
                ensemble=set(listPos)
                if abs(len(ensemble))>1:#Case of many position of mapping for NM => multiple mapped  Exit function
                    couple="multiple"
                    return(couple)
        if ensemble!=None: #In case of an ensemble==1 uniq genomic mapping position 
            couple="ok"
            return(couple)
    return(couple)
##############################################################
##############################################################
##############################################################
# posUp/posLow : dictionary with all the position associated with their mismatch number example : posUp[5687884684]=2
# returns 
#  * "." if the union of the positions of the best hit of the two sequences is of size 0 (if none of the two sequences is mapped)
#  * "PASS" if the union of the positions of the best hit of the two sequences is of size 1 (or zero, meaning that the best hit of one of the two sequences is higher than the authorized limit)
#  * "MULTIPLE" else (if the union of the positions of the best hit of the two sequences is of size >1)
##############################################################
def ValidationSNPBestHits (posUp,posLow):
    """Prediction validation : check if the couple is validated with only one mapping position """
  

    # get the best mapping distance for upper path 
    best_up=1024
    for position,nbMismatch in posUp.items(): 
        if nbMismatch<best_up:
            best_up=nbMismatch

    # get the best mapping distance for lower path 
    best_low=1024
    for position,nbMismatch in posLow.items(): 
        if nbMismatch<best_Low:
            best_Low=nbMismatch
    
    # get the union of the mapping position at the best mapping positions
    position_set = Set([])
    for position,nbMismatch in posUp.items():
        if nbMismatch == best_up:
            position_set.add(position)
            if len(position_set) > 1: 
                return "MULTIPLE"

    for position,nbMismatch in posLow.items():
        if nbMismatch == best_up:
            position_set.add(position)
            if len(position_set) > 1: 
                return "MULTIPLE"
    
    if len(position_set) == 1: 
        return "PASS"
    
    return "."
    
##############################################################
##############################################################
def ReverseComplement(nucleotide):
    if nucleotide=="A": return "T"
    if nucleotide=="T": return "A"
    if nucleotide=="C": return "G"
    return "C"
##############################################################
# cigarcode : in the samfile column 5 example 2S3M1I25M
# listpol : For close snps example listpol=[31,36,45] (with nb_pol=3)
# indel : boolean True if the current variant is an indel
##############################################################
def CigarCodeChecker(cigarcode,listpol):
    """Function which allows to recover the position of the SNPs according to the position of the deletions or insertions relative to the reference sequence (CigarCode Parsing : checks for insertion deletion or soft clipping"""
    parsingCigarCode=re.findall('(\d+|[A-Za-z])',cigarcode) #ParsingCigarCode=['2', 'S', '3', 'M', '1', 'I', '25', 'M']
    listPosRef=[]
    listShift=[]
    somme=0
    shift=0
    pos=0
    i=1
    j=0
    posCentraleRef=None
    #Reference (REF) =Genome
    #SEQ =query = discosnp path
    #Close snps
    if len(listpol)>1:#Checks : in case of close snps 
        while i<len(parsingCigarCode): #Goes through the list by twos to get all the letters and to take them into account
            #Soft clipping
            if parsingCigarCode[i]=="S":
                shift-=int(parsingCigarCode[i-1]) #It is the shift in the alignment between the reference and the sequence of the variant 
                pos+=int(parsingCigarCode[i-1]) #pos corresponds to the current position in the cigarcode
            #Match or Mismatch  
            elif parsingCigarCode[i]=="M":
                pos+=int(parsingCigarCode[i-1])
            #Deletion
            elif parsingCigarCode[i]=="D":
                shift+=int(parsingCigarCode[i-1])
            #Insertion
            elif parsingCigarCode[i]=="I":
                shift-=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            #Hard clipping (clipped sequences NOT present in SEQ)
            elif parsingCigarCode[i]=="H":
                shift-=int(parsingCigarCode[i-1]) #It is the shift in the alignment between the reference and the sequence of the variant 
                pos+=int(parsingCigarCode[i-1])
            #Padding (silent deletion from padded reference)
            elif parsingCigarCode[i]=="P":
                shift+=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="=":
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="X":
                pos+=int(parsingCigarCode[i-1])
            while int(pos)>=int(listpol[j]):#Goes through the list of position (close snps) and see if the position is affected by the shift (means shift before the position)
                posRef=int(listpol[j])+shift #Add the shift to the position (to get the real position of the snp on the reference)
                listPosRef.append(posRef) #Add the position to the list by taking into account the shift only if the current position 
                listShift.append(shift)
                if j<(len(listpol)-1):
                    j+=1
                else:
                    return(listPosRef,listShift)
            i+=2
    #Simple snp and Indel
    else:
        while i<len(parsingCigarCode):#Goes through the list by twos to get all the letters and to take them into account
            if parsingCigarCode[i]=="S":
                shift-=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="M":
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="D":
                shift+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="I":
                shift-=int(parsingCigarCode[i-1])#There is a nucleotide of shift compared to the reference
                pos+=int(parsingCigarCode[i-1]) #We advance in the query SEQ
            #Hard clipping (clipped sequences NOT present in SEQ)
            elif parsingCigarCode[i]=="H":
                shift-=int(parsingCigarCode[i-1]) # It's the shift in the alignment between the reference and the sequence of the variant 
                pos+=int(parsingCigarCode[i-1])
            #Padding (silent deletion from padded reference)
            elif parsingCigarCode[i]=="P":
                shift+=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="=":
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="X":
                pos+=int(parsingCigarCode[i-1])
            if len(listpol)==1:
                if pos>=int(listpol[0]):
                    posCentraleRef=int(listpol[0])+shift#takes into account the shift to add the after the mapping position and the variant position in the sequence
                    return(posCentraleRef,shift)
            i+=2
##############################################################
#shift : integer shift in the alignment between the reference and the sequence of the variant 
#posMut : column 19 in the samfile example 30T30 means 30 matches between the ref and the seq of the variant, one mismatch (T on the reference) and 30 matches between the ref and the seq of the variant
#posCentraleRef : position of the variant on the reference (taking the shift into account)
##############################################################
def ReferenceChecker(shift,posMut,posCentraleRef):
        """Function which allows to get the MD tag parsing; checks if path nucleotide is identical to the reference nucleotide"""
        i=0
        parsingPosMut=None
        boolDel=True
        pos=shift
        nucleoRef=None
        boolEgalRef=None
        matchInt=None
        parsingPosMut=re.findall('(\d+|[A-Za-z]|\^)',posMut)
        while i<len(parsingPosMut):#Converts the number into integer 
                matchInt=re.match(r'\d+',parsingPosMut[i]) #Integer motif
                if matchInt:
                        parsingPosMut[i]=int(parsingPosMut[i])
                i+=1
        i=0
        while i<len(parsingPosMut): #Goes through the list to know if the variant is identical to the reference
                if isinstance(parsingPosMut[i],int): #Checks if it's a number of nucleotide or a letter
                        pos+=parsingPosMut[i] #Adds to pos the current number of the list
                elif parsingPosMut[i]=="^":#In case of deletion from the reference we have to substract the deletion from the current position  
                        i+=1
                        while boolDel: #Checks if it is still a letter in the deletion to substract it from the current position
                                if isinstance(parsingPosMut[i],int): #If it is an integer we achieved to take into account the deletion : adds the integer to the current position
                                        boolDel=False
                                        pos+=parsingPosMut[i]
                                else:
                                        pos-=1 #Substracts the nucleotide from the current position
                                        i+=1
                else:
                        pos+=1
                if pos==posCentraleRef: # Checks if the current position pos is identical to the position of the variant 
                        if isinstance(parsingPosMut[i],str): #=> it means that the nucleotide is different in the variant and in the reference
                                boolEgalRef=False
                                nucleoRef=parsingPosMut[i]
                        else: #If the last item of the list of the MD tag is an intger => it means that the nucleotide of the allele is identical to the reference
                                boolEgalRef=True
                        break
                if pos>posCentraleRef: #If the current position is bigger than the variant position it means that the nucleotide of the variant is identical to the reference
                        boolEgalRef=True
                        break
                i+=1
        return(boolEgalRef,nucleoRef)   
    # boolEgalRef : boolean if TRUE the nucleotide in the variant is equal to the reference
    # nucleoRef : if the nucleotide is different from the reference return the nucleotide of reference
##############################################################
#snpUp/snpLow: sam line
##############################################################
def GetSequence(snpUp,snpLow):
    """Gets and reverses if it is needed sequence"""
    seqLow=snpLow[9]
    seqUp=snpUp[9]
    listSeqUp=[]
    listSeqLow=[]
    i=0
    if CheckBitwiseFlag(snpUp[1],4): #Case of reverse mapping 
        i=0
        listSeqUp=list(seqUp)
        seqUp=''
        #Sequence of the upper path
        while i<len(listSeqUp): #Reverse the sequence : usefull for future comparison and to get the position of an eventual indel
            if seqUp!='':#If the sequence is not empty adds the following reverse nucleotide
                seqUp=str(ReverseComplement(listSeqUp[i]))+seqUp
            else : #Starts the reverse sequence with the first reverse nucleotide 
                seqUp=str(ReverseComplement(listSeqUp[i]))
            i+=1
    if CheckBitwiseFlag(snpLow[1],4):#Case of mapping reverse
        i=0
        listSeqLow=list(str(seqLow))
        seqLow=''
        #Sequence of the lower path
        while i<len(listSeqLow):
            if seqLow!='':
                seqLow=str(ReverseComplement(listSeqLow[i]))+seqLow
            else :
                seqLow=str(ReverseComplement(listSeqLow[i]))
            i+=1
    return(seqUp,seqLow)
    #Return the two sequences in direction forward 
##############################################################
#dicoHeader snps : dicoHeader[key]=[posD,ntUp,ntLow] 
#dicoHeader indel : dicoHeader[key]=[posD,ind,amb]
#snpUp/snpLow : line of the sam file
#posUp/posLow : dictionary with all the mapping position of every snps associated with their number of mismatch
#nb_polUp/nb_polLow : number of variants for the current path
#indel : boolean True if the current variant is an indel
##############################################################
def RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel):
    """Determines : if the paths are mapped in reverse or forward sense ; which paths have the variant identical to the genome (information from the MD Tag) ; the position of mapping of the variant for each path by taking into account the shift with the genome ( soft clipping insertion, deletion  and so on : information from the cigar code) ; the nucleotide variant for the upper path and the lower path"""
    #INIT VALUES
    posSNPUp=dicoHeaderUp["P_1"][0] 
    posSNPLow=dicoHeaderUp["P_1"][0]
    boolRefUp=None
    boolRefLow=None
    nucleoUp=None
    nucleoLow=None
    nucleoRefUp=None
    nucleoRefLow=None
    listPolymorphismPosUp=None
    listPolymorphismPosLow=None
    shiftUp=0
    shiftLow=0
    posCentraleUp=None
    posCentraleLow=None
    listPos=None
    listnucleoUp=None
    listnucleoLow=None
    listPosR=None
    listnucleoUpR=None
    listnucleoLowR=None
    seqUp=None
    seqLow=None
    reverseUp=None
    reverseLow=None
    MD=None
    Z=None
    posMutUp=None
    posMutLow=None
    ambiguityPos=None
    insert=None
    ntStart=None
    seqUp=list(snpUp[9])
    seqLow=list(snpLow[9])
    reverseUp=1
    reverseLow=1
    dicopolUp={}
    dicopolLow={}
    boolSmallestUp=False
    boolSmallestLow=False
    i=0
    #Defines the first value of the reference position (first given by bwa) : will be usefull when the mapping position of each path is the same
    if int(snpUp[3])>0:
        posRef=int(snpUp[3])
    else:
        posRef=int(snpLow[3])
#---------------------------------------------------------------------------------------------------------------------------
    #List of snps : Gets the positions of the variant in reverse and forward position ; gets the nucleotide for each allele for the forward and reverse strand
    if indel==False:
        listPos,listnucleoUp,listnucleoLow,listPosR,listnucleoUpR,listnucleoLowR=GetPolymorphism(dicoHeaderUp,snpUp[9],indel,False)
    else: #In case of Indel we will get the insertion in the longest sequence 
        seqUp=snpUp[9]
        seqLow=snpLow[9]          
        if len(seqUp)<len(seqLow): #Determines the longest sequence to get the insertion
            boolSmallestUp=True
            boolSmallestLow=False
        else:
            boolSmallestUp=False
            boolSmallestLow=True
        listPos,listPosRUp,insert,ntStart,ambiguityPos=GetPolymorphism(dicoHeaderUp,seqUp,indel,boolSmallestUp) #For indel get the insert, the list of position and the possible ambiguity for the position
        listPos,listPosRLow,insert,ntStart,ambiguityPos=GetPolymorphism(dicoHeaderUp,seqLow,indel,boolSmallestLow)
#---------------------------------------------------------------------------------------------------------------------------
#Gets the shift by positions (insertion,deletion,sofclipping) and update of the position on the path
#---------------------------------------------------------------------------------------------------------------------------
###Case mapped variant Up : Gets the shift by positions (insertion,deletion,sofclipping) 
    if int(snpUp[3])>0: #Mapped path
        #Check cigarCode : Presence of insertion, softclipping, deletion
        if not CheckBitwiseFlag(snpUp[1],4) : #Forward strand : check the bitwise flag of the samfile
            posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],listPos) #Gets the positions of the variants with an eventual shift from the reference (insertion,deletion,soft clipping) ; in case of close snps return a list 
            listPolymorphismPosUp=listPos
            if len(listPos)==1 and indel==False: #Simple snp
                nucleoUp=listnucleoUp[0] #Gets the nucleotide of the allele
        elif CheckBitwiseFlag(snpUp[1],4):#Reverse Strand
            reverseUp=-1
            if indel==True:
                posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],listPosRUp)
                listPolymorphismPosUp=listPosRUp #List of all the reverse positions
            else:
                posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],listPosR)
                listPolymorphismPosUp=listPosR #List of all the reverse positions
                if len(listPos)==1 and indel==False:
                        nucleoUp=listnucleoUpR[0]
        else:#Default case we can not determine the bitwise flag.
            posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],listPos)
            listPolymorphismPosUp=listPos
            if len(listPos)==1 and indel==False: #Simple snp
                nucleoUp=listnucleoUp[0] #Gets the nucleotide 
    else: ###Case unmapped path
        listPolymorphismPosUp=listPos
        posCentraleUp=listPos
        shiftUp=None
        reverseUp="."
#---------------------------------------------------------------------------------------------------------------------------
###Case mapped variant Low : Gets the shift by positions(insertion,deletion,sofclipping)
    if int(snpLow[3])>0:
        #Check cigarCode : Presence of insertion, softclipping, deletion
        if not CheckBitwiseFlag(snpLow[1],4):#Forward Strand
            posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],listPos)
            listPolymorphismPosLow=listPos
            if len(listPos)==1 and indel==False:#Defines the nucleotide of the allele in case of forward strand and simple snp
                nucleoLow=listnucleoLow[0]
        elif CheckBitwiseFlag(snpLow[1],4):#Reverse strand
            reverseLow=-1
            if indel==True:
                posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],listPosRLow)
                listPolymorphismPosLow=listPosRLow
            else:
                posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],listPosR)
                listPolymorphismPosLow=listPosR
            if len(listPos)==1 and indel==False:#Defines the nucleotide of the allele in case of reverse strand and simple snp
                nucleoLow=listnucleoLowR[0]
        else:
            posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],listPos)
            listPolymorphismPosLow=listPos
            if len(listPos)==1 and indel==False:
                nucleoLow=listnucleoLow[0]             
    else:###Case unmapped variant Low
        listPolymorphismPosLow=listPos
        posCentraleLow=listPos
        shiftLow=None
        reverseLow="."
#---------------------------------------------------------------------------------------------------------------------------
#Checks for all the variants if there are identical to the reference by parsing the MD tag ; and defines the nucleotide on the reference
#---------------------------------------------------------------------------------------------------------------------------
    ####Upper Case : Checks for all the variants if there are identical to the reference 
    if int(snpUp[3])>0:
        posMutUp,nbMismatchUp=GetTag(snpUp) #MD tag parsing
        if len(listPos)>1: #CLOSE SNPS
            if reverseUp==1:
                i=0
                for i in range(len(listPos)):#In case of close snps : goes through the list of allele position
                    boolRefUp,nucleoRefUp=ReferenceChecker(shiftUp[i],posMutUp,posCentraleUp[i]) #Checks if the variant is identical to the reference ; returns a boolean and the nucleotide of the reference
                    if nucleoRefUp==None: #If there is no reference nucleotide given by ReferenceChecker, it means that the variant is equal to the reference so we defined it !
                        nucleoRefUp=listnucleoUp[i] 
                    dicopolUp[listPos[i]]=[boolRefUp,nucleoRefUp,posCentraleUp[i],listnucleoUp[i],reverseUp,(int(snpUp[3])+posCentraleUp[i])] #dictionary for close snps to keep all the results given by the different functions
            elif reverseUp==-1:
                i=0
                for i in range(len(listPosR)):#In case of close snps with reverse mapping : goes through the list of reverse allele position
                    boolRefUp,nucleoRefUp=ReferenceChecker(shiftUp[i],posMutUp,posCentraleUp[i])
                    if nucleoRefUp==None: 
                        nucleoRefUp=listnucleoUpR[i]
                    dicopolUp[listPosR[i]]=[boolRefUp,nucleoRefUp,posCentraleUp[i],listnucleoUpR[i],reverseUp,(int(snpUp[3])+posCentraleUp[i])]                    
        else: #SIMPLE SNP
            boolRefUp,nucleoRefUp=ReferenceChecker(shiftUp,posMutUp,posCentraleUp)
            if nucleoRefUp==None and reverseUp==1:#If there is no reference nucleotide given by ReferenceChecker, it means that the variant is equal to the reference so we defined it !
                nucleoRefUp=dicoHeaderUp["P_1"][1]
            elif nucleoRefUp==None and reverseUp==-1:
                nucleoRefUp=ReverseComplement(dicoHeaderUp["P_1"][1])
	    
    
    elif int(snpUp[3])<=0 :
        nucleoUp = dicoHeaderUp["P_1"][1]
        if reverseLow==-1:
               nucleoUp = ReverseComplement(dicoHeaderUp["P_1"][1])
#---------------------------------------------------------------------------------------------------------------------------
    ####Lower Case 
    if int(snpLow[3])>0:
        #Extraction tag MD
        posMutLow,nbMismatchLow=GetTag(snpLow)
        if len(listPos)>1:
            if reverseLow==1:
                i=0
                for i in range(len(listPos)):
                    boolRefLow,nucleoRefLow=ReferenceChecker(shiftLow[i],posMutLow,posCentraleLow[i])
                    if nucleoRefLow==None:
                        nucleoRefLow=listnucleoLow[i]
                    dicopolLow[listPos[i]]=[boolRefLow,nucleoRefLow,posCentraleLow[i],listnucleoLow[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
            elif reverseLow==-1:
                i=0
                for i in range(len(listPosR)):
                    boolRefLow,nucleoRefLow=ReferenceChecker(shiftLow[i],posMutLow,posCentraleLow[i])
                    if nucleoRefLow==None:
                        nucleoRefLow=listnucleoLowR[i]
                    dicopolLow[listPosR[i]]=[boolRefLow,nucleoRefLow,posCentraleLow[i],listnucleoLowR[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
        else:
            boolRefLow,nucleoRefLow=ReferenceChecker(shiftLow,posMutLow,posCentraleLow)
            if nucleoRefLow==None and reverseLow==1:#If there is no reference nucleotide given by ReferenceChecker, it means that the variant is equal to the reference so we defined it !
                nucleoRefLow=dicoHeaderUp["P_1"][2]
            elif nucleoRefLow==None and reverseLow==-1:
                nucleoRefLow=ReverseComplement(dicoHeaderUp["P_1"][2])
    
    elif int(snpLow[3])<=0:
        nucleoLow = dicoHeaderUp["P_1"][2]
        if reverseUp==-1:
               nucleoLow = ReverseComplement(dicoHeaderUp["P_1"][2])

#---------------------------------------------------------------------------------------------------------------------------
#Checks the exception : different mapping position or both paths identical to the reference
#---------------------------------------------------------------------------------------------------------------------------
####Two paths mapped at different main position (first given by bwa) : test which path will be the reference with the function mismatchChecker
    if ((int(snpUp[3])>0 and int(snpLow[3])>0) and (int(snpUp[3])!=int(snpLow[3]))) or (int(snpUp[3])>0 and int(snpLow[3])>0 and boolRefUp==True and boolRefLow==True):
        if len(listPos)>1:
            i=0
            for i in range(len(listPos)):
                boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,dicopolUp[listPolymorphismPosUp[i]][1],dicopolLow[listPolymorphismPosLow[i]][1],dicopolUp[listPolymorphismPosUp[0]][3],dicopolLow[listPolymorphismPosLow[0]][3],dicopolUp[listPolymorphismPosUp[i]][0],dicopolLow[listPolymorphismPosLow[i]][0],indel)
                dicopolUp[listPolymorphismPosUp[i]][0]=boolRefUp
                dicopolUp[listPolymorphismPosUp[i]][1]=nucleoRefUp
                dicopolLow[listPolymorphismPosLow[i]][0]=boolRefLow
                dicopolLow[listPolymorphismPosLow[i]][1]=nucleoRefLow
                dicopolLow[listPolymorphismPosLow[i]][5]=int(dicopolLow[listPolymorphismPosLow[i]][2])+int(posRef)
                dicopolUp[listPolymorphismPosUp[i]][5]=int(dicopolUp[listPolymorphismPosUp[i]][2])+int(posRef)
        else:
                boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel)
#---------------------------------------------------------------------------------------------------------------------------
#Defines the mapping position 
#---------------------------------------------------------------------------------------------------------------------------
    #Variants (INDEL and SIMPLE SNPS) Positions : adds the position of mapping given by bwa to the variant position and for indel add ambiguity position
    if int(nb_polUp)==1: #Checks if there is only one allele in the path
        if indel==True: #Case of indel 
            nucleoLow="." #It will be define later
            nucleoUp="."
            nucleoRefUp="." #No nucleotide for the reference
            nucleoRefLow="."
        if (int(snpUp[3])>0 and int(snpLow[3])>0): #Checks if both paths are mapped
            if indel==True: #In case of indel : defines which path will be considered as the reference in function of the mapping position
                posSNPUp=posCentraleUp+posRef-int(ambiguityPos) #taking into account the possible ambiguity for the indel position
                posSNPLow=posCentraleLow+posRef-int(ambiguityPos)
                if snpUp[9]<snpLow[9]: #Checks with path has the lefmost position to determine the reference
                    boolRefUp=True
                    boolRefLow=False
                else:
                    boolRefUp=False
                    boolRefLow=True
            else: #In case of simple snps
                posSNPUp=posCentraleUp+posRef #Position of the variant on the path (+shift from the reference)+ mapping position given by bwa
                posSNPLow=posCentraleLow+posRef
        elif int(snpUp[3])<=0 and int(snpLow[3])>0:
            if indel==True:
                boolRefUp=False
                boolRefLow=True
                posSNPUp=None
                posSNPLow=posCentraleLow+posRef-int(ambiguityPos)
            else:
                posSNPUp=None
                posSNPLow=posCentraleLow+posRef
        elif int(snpUp[3])>0 and int(snpLow[3])<=0: #Upper mapped path  ; lower unmapped path 
            if indel==True:
                boolRefUp=True
                boolRefLow=False
                posSNPLow=None
                posSNPUp=posCentraleUp+posRef-int(ambiguityPos)
            else:
                posSNPLow=None
                posSNPUp=posCentraleUp+posRef
        elif (int(snpUp[3])<=0 and int(snpLow[3])<=0): #Both unmapped paths
                posSNPUp=dicoHeaderUp["P_1"][0]
                posSNPLow=dicoHeaderUp["P_1"][0]
#---------------------------------------------------------------------------------------------------------------------------        
        ##Checks the strand : if the alternative nucleotide is not on the same strand as the reference : reverse it
        if boolRefLow==True and reverseLow==-1 and reverseUp==1:
               nucleoUp=ReverseComplement(nucleoUp)
        elif boolRefLow==True and reverseLow==1 and reverseUp==-1:
               nucleoUp=ReverseComplement(nucleoUp)
        if boolRefUp==True and reverseUp==-1 and reverseLow==1:
               nucleoLow=ReverseComplement(nucleoLow)
        elif boolRefUp==True and reverseUp==1 and reverseLow==-1:
               nucleoLow=ReverseComplement(nucleoLow)
#---------------------------------------------------------------------------------------------------------------------------
        ##Checks which path is chosen : if there are both different from the reference=> use the function MismatchChecker
        if boolRefUp==False and boolRefLow==False:
               boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel)
        return(nucleoLow,posSNPLow,nucleoUp,posSNPUp,boolRefLow,boolRefUp,reverseUp,reverseLow,nucleoRefUp,nucleoRefLow)
        #nucleoLow/nucleoUp : nucleotide of each path corresponding to the variant
        #posSNPLow/posSNPUp : position of the variant taking into account the shift
        #boolRefUp/boolRefLow : boolean TRUE the variant (of the path Up or Low) is the reference
    else:
        return(dicopolUp,dicopolLow,listPolymorphismPosUp,listPolymorphismPosLow)
##############################################################
##############################################################
def MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel):
    """In case of divergent main position (case snp whose two paths are mapped ) to define the reference = > check the number of mismatch
( If the number of mismatch is the same in both cases it is the lower lexicographical SNP which is selected for reference .
The Boolean allows to know the reference SNP ) """
    posRef=0
    nmUp=None
    nmLow=None
    #Two paths mapped
    if int(snpUp[3])>0 and int(snpLow[3])>0:#Checks if both paths are mapped 
        nmUp=posUp[int(snpUp[3])] #Distance with the reference for the snpUp
        nmLow=posLow[int(snpLow[3])] #Distance with the reference for the snpLow
        if nmUp<nmLow: #Checks if the upper path has a distance with the reference smaller than the lower path
            boolRefLow=False
            boolRefUp=True 
            nucleoRefLow=nucleoRefUp
            posRef=int(snpUp[3])
        elif nmUp>nmLow :#Checks if the lower path has a distance with the reference smaller than the upper path
            boolRefLow=True
            boolRefUp=False
            nucleoRefUp=nucleoRefLow
            posRef=int(snpLow[3])
        elif nmUp==nmLow: #Checks if both path have the same number of difference
            if indel==False: #In case of simple snp
                if nucleoUp<nucleoLow: #Checks the lexicographical order
                    boolRefLow=False
                    boolRefUp=True
                    nucleoRefLow=nucleoRefUp
                    posRef=int(snpUp[3])
                elif nucleoUp>nucleoLow: #Checks the lexicographical order
                    boolRefLow=True
                    boolRefUp=False
                    nucleoRefUp=nucleoRefLow
                    posRef=int(snpLow[3])
                else : #If none of the alleles is lexicographically less : checks the mapping position and keeps the lefmost position
                    if int(snpUp[3])<int(snpLow[3]):
                        boolRefLow=False
                        boolRefUp=True
                        nucleoRefLow=nucleoRefUp
                        posRef=int(snpUp[3])
                    else:
                        boolRefLow=True
                        boolRefUp=False
                        nucleoRefUp=nucleoRefLow
                        posRef=int(snpLow[3])
            else: #In case of indel
                if int(snpUp[3])<int(snpLow[3]): #Checks the mapping position and keeps the lefmost position
                    boolRefLow=False
                    boolRefUp=True
                    nucleoRefLow=nucleoRefUp
                    posRef=int(snpUp[3])
                else:
                    boolRefLow=True
                    boolRefUp=False
                    nucleoRefUp=nucleoRefLow
                    posRef=int(snpLow[3])
    return(boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef)
##############################################################
#Simple Variant : dicoHeader[key]=[posD,ntUp,ntLow]
#Indel : dicoHeader[key]=[posD,ind,amb]
##############################################################
def GetPolymorphism(dicoHeader,seq,indel,boolSmallest):
    '''Gets from the dicoHeader all the positions, and the nucleotides (R means that it's on the reverse strand)
      SNP :  one variant correspond to listPos[0], listPosR[0], listnucleoUp[0], listnucleoLow[0],listnucleoUpR[0],listnucleoLowR[0]'''
    #Forward variables
    listPos=[]
    listnucleoUp=[]
    listnucleoLow=[]
    #Reverse variables
    listPosR=[]
    listnucleoUpR=[]
    listnucleoLowR=[]
    tailleSeq=len(seq)
    insert=None
    ntStart=None
    if indel==False:##Case of simple snp
        for key,(posD,ntUp,ntLow) in dicoHeader.items(): #Goes through the dictionary of parsed header
            listPos.append(int(posD)) #Adds the position of the variant
            listPosR.append(tailleSeq-int(posD)) #Adds the reverse position of the variant
            listnucleoUp.append(ntUp) #Adds the allele of the upper path
            listnucleoUpR.append(ReverseComplement(ntUp))#Adds the reverse allele of the upper path
            listnucleoLow.append(ntLow)#Adds the allele of the lower path
            listnucleoLowR.append(ReverseComplement(ntLow))#Adds the reverse allele of the lower path
        return(listPos,listnucleoUp,listnucleoLow,listPosR,listnucleoUpR,listnucleoLowR)
    else:##Case of indel
        for key,(posD,ind,amb) in dicoHeader.items():#Goes through the dictionary of parsed header
            listPos.append(posD+1)#Adds the position of the variant
            listPosR.append(tailleSeq-int(posD)+1)#Adds the reverse position of the variant
            if boolSmallest==False: #If we have teh sequence of the longest path : gets the insert
                insert=seq[(int(posD-1)-1):(int(posD-1)+int(ind))] #Gets the insert with the position on the variant (just on forward sequence)
                ntStart=seq[(int(posD-1)-1)] #Get the nucleotide just before the insertion
            ambiguityPos=amb #Gets the possible ambiguity for the insertion
        return(listPos,listPosR,insert,ntStart,ambiguityPos)
##############################################################
##############################################################
def fillVCFSimpleSnp(snpUp,snpLow,nucleoLow,positionSnpLow,nucleoUp,positionSnpUp,boolRefLow,boolRefUp,table,nbSnp,filterfield,ok,tp,phased,listCovGeno,nucleoRefUp,nucleoRefLow,reverseUp,reverseLow,geno,nbGeno,covUp,covLow):
    """ Fills the different fields of vcf based on boolean ( on whether the SNP is identical or not the reference) , if neither is identical to the reference = > we take the SNP comes first in the lexicographical order"""
#---------------------------------------------------------------------------------------------------------------------------    
    ##Gets the variable of the header of disco snps
    ##Upper path
    numSNPUp=None
    unitigLeftUp=None
    unitigRightUp=None
    contigLeftUp=None
    contigRightUp=None
    valRankUp=None
    listCoverageUp=None
    listCUp=None
    nb_polUp=None
    lnUp=None
    posDUp=None
    genoUp=None
    dicoHeaderUp=None
    ntUp=None
    posUnmappedUp=None
    ##Lower path
    numSNPLow=None
    unitigLeftLow=None
    unitigRightLow=None
    contigLeftLow=None
    contigRightLow=None
    valRankLow=None
    listCoverageLow=None
    listClow=None
    nb_polLow=None
    lnlow=None
    posDLow=None
    ntLow=None
    genoLow=None
    dicoHeaderLow=None
#---------------------------------------------------------------------------------------------------------------------------
##Parsing of discosnp++ header to fill all the vcf fields        
    discoNameUp,snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp, listCoverageUp, listCUp,nb_polUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
    discoNameLow,snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow, listCoverageLow,listClow,nb_polLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
    posUnmappedUp=CheckContigUnitig(unitigLeftUp,contigLeftUp) #Takes into account the lenght of the unitig/contig for the position of unmapped allele (position of the allele on the upper path)
    posUnmappedLow=CheckContigUnitig(unitigLeftLow,contigLeftLow) #Takes into account the lenght of the unitig/contig for the position of unmapped allele (position of the allele on the lower path)
#---------------------------------------------------------------------------------------------------------------------------
##Case : two mapped paths
    if int(snpUp[3])>0 and int(snpLow[3])>0:
        ##The path identical to the reference is the lower path 
        if boolRefLow==True and boolRefUp==False:
            table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,".",filterfield,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
            
        ##The path identical to the reference is the upper path 
        elif boolRefUp==True and boolRefLow==False:
            table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,".",filterfield,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
        ##No path is identical to the reference => lexicographique choice
        elif boolRefUp==False and boolRefLow==False:
            if nucleoUp<nucleoLow:
                table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,".",filterfield,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
            elif nucleoUp>nucleoLow:
              table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,".",filterfield,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
            
              
#---------------------------------------------------------------------------------------------------------------------------
##Case : Lower path mapped and upper path unmapped      
    elif int(snpUp[3])<=0 and int(snpLow[3])>0:
        table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,".",filterfield,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
        if "."=="*":
                table[5]="."
        else:
                table[5]="."  
#---------------------------------------------------------------------------------------------------------------------------
##Case : Upper path mapped and lower path unmapped        
    elif int(snpUp[3])>0 and int(snpLow[3])<=0:
        table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,".",filterfield,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
        if "."=="*":
                table[5]="."
        else:
                table[5]="."
#---------------------------------------------------------------------------------------------------------------------------
##Case : Both paths are unmapped       
    elif int(snpUp[3])<=0 and int(snpLow[3])<=0:#FillVCF(table,numSNP,chrom,pos,ref,alt,qual,filterfield,tp,valRank,ok,unitigLeft,unitigRight,contigLeft,contigRight,cov,nucleoRef,reverse,geno,nbGeno,phased,listCovGeno,boolRefLow)
        if nucleoLow<nucleoUp:
            table=FillVCF(table,numSNPLow,discoNameUp,(int(positionSnpUp)+int(posUnmappedUp)),nucleoLow,nucleoUp,".",filterfield,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,".",".",geno,nbGeno,phased,listCovGeno,boolRefLow)
        else:
            table=FillVCF(table,numSNPUp,discoNameLow,(int(positionSnpLow)+int(posUnmappedLow)),nucleoUp,nucleoLow,".",filterfield,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,".",".",geno,nbGeno,phased,listCovGeno,boolRefLow)
        
    return(table)
##############################################################
##############################################################   
def printOneline(table,VCF):
    """Prints the line of the current SNP in the VCF file."""
    if table==[0, 0, 0, 0, 0, 0, 0, 0, 0, 0] or (table[0]==0 and table[1]==0 and table[3]==0 and table[4]==0):
        return
    for i in range(len(table)):
        element=table[i]
        VCF.write((str(element)).strip())
        if i<len(table)-1: VCF.write("\t")
    VCF.write('\n')    

##############################################################
#dicopolLow[listPos[i]]=[boolRefLow,nucleoRefLow,posCentraleLow[i],listnucleoLowR[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
##############################################################
def printVCFSNPclose(dicoUp,dicoLow,table,filterField,snpUp,snpLow,listPolymorphismPosUp,listPolymorphismPosLow,ok,covUp,covLow,listnucleoUp,listnucleoLow,geno,nbGeno,listCovGeno, VCF,posUp,posLow):
    """Prints the line of the close SNPs in the VCF file. The first variant will determine which path will be used as reference (checks with the boolean boolRef which path is identical to the reference) """ 
    info=''
    champAlt=0
    comptPol=0
    seqUp=snpUp[9]
    seqLow=snpLow[9]
    tp="SNP"
    phased=True
    table = [0] * 10 #Create a 10 cols array
    tablebis = []
    k=0
    listSortedPosUp=None
    indexSmallestPosUp=None
    listSortedPosLow=None
    indexSmallestPosLow=None
    indexSmallestPos=None
    numSNPUp=None
    unitigLeftUp=None
    unitigRightUp=None
    contigLeftUp=None
    contigRightUp=None
    valRankUp=None
    listCoverageUp=None
    listCUp=None
    nb_polUp=None
    lnUp=None
    posDUp=None
    ntUp=None
    ntLow=None
    genoUp=None
    dicoHeaderUp=None
    numSNPLow=None
    unitigLeftLow=None
    unitigRightLow=None
    contigLeftLow=None
    contigRightLow=None
    valRankLow=None
    listCoverageLow=None
    listClow=None
    nb_polLow=None
    lnlow=None
    posDLow=None
    ntUp=None
    ntLow=None
    genoLow=None
    dicoHeaderLow=None
    positionSnpUp1=None
    nucleoUp1=None
    nucleoLow1=None
    positionSnpLow1=None
    reverseLow=None
    reverseUp=None
    boolRefUp=None
    boolRefLow=None
    positionSnpUp=None
    nucleoUp=None
    nucleoLow=None
    positionSnpLow=None
    ID=0
    comptPol=0
    posUnmappedUp=None
    #Variables
    discoNameUp,snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp, listCoverageUp, listCUp,nb_polUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
    discoNameLow,snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow, listCoverageLow, listCLow,nb_polLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
    posUnmappedUp=CheckContigUnitig(unitigLeftUp,contigLeftUp)
    posUnmappedLow=CheckContigUnitig(unitigLeftLow,contigLeftLow)
#---------------------------------------------------------------------------------------------------------------------------
##Case : two mapped paths
    if int(snpUp[3])>0 and int(snpLow[3])>0:
        ##Sorts the list of position to get the smallest one and its position in the list unsorted!
        ##Keeps the close snps to sort them : indeed all the lists and dictionnaries :listnucleoUp,listPolymorphismPosUp,listPolymorphismPosLow,listnucleoLow dicoUp,dicoLow are classified according to dicoHeader so if we start by sorting we lose the correspondence between data
        listSortedPosUp=list(listPolymorphismPosUp)
        listSortedPosUp.sort()
        indexSmallestPosUp=listPolymorphismPosUp.index(listSortedPosUp[0])
        listSortedPosLow=list(listPolymorphismPosLow)
        listSortedPosLow.sort()
        indexSmallestPosLow=listPolymorphismPosLow.index(listSortedPosLow[0])
        boolRefUp=dicoUp[listPolymorphismPosUp[indexSmallestPosUp]][0]
        boolRefLow=dicoLow[listPolymorphismPosLow[indexSmallestPosLow]][0]
        #Decides what is the smallest position according to the reference path (useful if the paths are not aligned on the same strand)
        if boolRefUp==True and boolRefLow==True:
            #dicopolLow[listPos[i]]=[boolRefLow,nucleoRefLow,posCentraleLow[i],listnucleoLowR[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
            #boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel)
            boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,dicoUp[listPolymorphismPosUp[indexSmallestPosUp]][1],dicoLow[listPolymorphismPosLow[indexSmallestPosLow]][1],dicoUp[listPolymorphismPosUp[indexSmallestPosUp]][3],dicoLow[listPolymorphismPosLow[indexSmallestPosLow]][3],boolRefUp,boolRefLow,False)    
        if boolRefUp==True and boolRefLow==False: #The smallest position identical to the reference is on the upper path
            indexSmallestPos=indexSmallestPosUp
            listPolymorphismPos=listPolymorphismPosUp
        elif boolRefLow==True and boolRefUp==False:#The smallest position identical to the reference is on the lower path
            indexSmallestPos=indexSmallestPosLow
            listPolymorphismPos=listPolymorphismPosLow
        elif boolRefUp==False and boolRefLow==False: #Both paths are different from the reference => choice with the lexicographical order
            nucleoUp1=dicoUp[listPolymorphismPosUp[indexSmallestPosUp]][3]
            nucleoLow1=dicoLow[listPolymorphismPosLow[indexSmallestPosLow]][3]
            if nucleoUp1<nucleoLow1:
                indexSmallestPos=indexSmallestPosUp
                listPolymorphismPos=listPolymorphismPosUp
            else:
                indexSmallestPos=indexSmallestPosLow
                listPolymorphismPos=listPolymorphismPosLow
        #Remembers the values for the first snps ==> the others snps will depend on these parameters
        boolRefUp=dicoUp[listPolymorphismPosUp[indexSmallestPos]][0]
        boolRefLow=dicoLow[listPolymorphismPosLow[indexSmallestPos]][0]
        #Checks if both allele are different of the reference => decides which variant to choose as reference 
        if boolRefUp==True and boolRefLow==True:
            #dicopolLow[listPos[i]]=[boolRefLow,nucleoRefLow,posCentraleLow[i],listnucleoLowR[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
            #boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel)
            boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,dicoUp[listPolymorphismPosUp[indexSmallestPos]][1],dicoLow[listPolymorphismPosLow[indexSmallestPos]][1],dicoUp[listPolymorphismPosUp[indexSmallestPos]][3],dicoLow[listPolymorphismPosLow[indexSmallestPos]][3],boolRefUp,boolRefLow,False)  
        nucleoUp1=dicoUp[listPolymorphismPosUp[indexSmallestPos]][3]
        nucleoLow1=dicoLow[listPolymorphismPosLow[indexSmallestPos]][3]
        positionSnpUp1=dicoUp[listPolymorphismPosUp[indexSmallestPos]][5]
        positionSnpLow1=dicoLow[listPolymorphismPosLow[indexSmallestPos]][5]
        reverseLow=dicoLow[listPolymorphismPosLow[indexSmallestPos]][4]
        reverseUp=dicoUp[listPolymorphismPosUp[indexSmallestPos]][4]
        for comptPol in range(len(listPolymorphismPos)): #Goes through the list of the variant position starting with the smallest
            positionSnpUp=dicoUp[listPolymorphismPosUp[comptPol]][5]
            positionSnpLow=dicoLow[listPolymorphismPosLow[comptPol]][5]
            nucleoUp=dicoUp[listPolymorphismPosUp[comptPol]][3]
            nucleoRefUp=dicoUp[listPolymorphismPosUp[comptPol]][1]
            nucleoLow=dicoLow[listPolymorphismPosLow[comptPol]][3]
            nucleoRefLow=dicoLow[listPolymorphismPosLow[comptPol]][1]
            #Fills the variable table with the vcf fields ; Checks the "REF" path to fill the vcf
            if boolRefLow==True and boolRefUp==False: #The lower path is defined as REF
                 table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,".",filterField,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
                 table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoUp)
            elif boolRefUp==True and boolRefLow==False:#The upper path is defined as REF
                table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,".",filterField,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
                table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoLow)
            elif boolRefUp==False and boolRefLow==False: #The two paths are different from the reference => defines which one will be the reference (thanks to the first allele of the path)
                if nucleoUp1<nucleoLow1:
                    table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,".",filterField,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
                    table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoLow)
                elif  nucleoUp1>nucleoLow1:
                      table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,".",filterField,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
                      table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoUp)
                elif positionSnpUp1<positionSnpLow1:
                    table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,".",filterField,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
                    table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoLow)
                elif positionSnpUp1>positionSnpLow1:
                    table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,".",filterField,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
                    table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoUp)
            tablebis.append(list(table))#Stocks the variable with all the vcf fields for each close snp to sort it and print it in the vcf
        tablebis=sorted(tablebis, key=lambda colonnes: colonnes[1])
        l=0
        ID=1
        for l in range(len(tablebis)):
                tablebis[l][2]=str(table[2])+"_"+str(ID)
                printOneline(tablebis[l],VCF)
                ID+=1
#---------------------------------------------------------------------------------------------------------------------------
##Case : Both paths are unmapped     
    elif int(snpUp[3])<=0 and int(snpLow[3])<=0:
        reverseUp="."
        i=0
        for i in range(len(listPolymorphismPosUp)):
            positionSnpUp=int(listPolymorphismPosUp[i])+int(posUnmappedUp)
            positionSnpLow=int(listPolymorphismPosUp[i])+int(posUnmappedLow)
            nucleoUp=listnucleoUp[i]
            nucleoRefUp="."
            nucleoLow=listnucleoLow[i]
            nucleoRefLow="."
            table=FillVCF(table,numSNPUp,discoNameUp,positionSnpUp,nucleoUp,nucleoLow,".",filterField,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,0)
            tablebis.append(list(table))
        tablebis=sorted(tablebis, key=lambda colonnes: colonnes[1])
        l=0
        ID=1
        for l in range(len(tablebis)):
                tablebis[l][2]=str(table[2])+"_"+str(ID)
                printOneline(tablebis[l],VCF)
                ID+=1
#---------------------------------------------------------------------------------------------------------------------------
##Case : Upper path mapped and lower path unmapped  
    elif int(snpUp[3])>0 and int(snpLow[3])<=0:
        comptPol=0
        reverseUp=dicoUp[listPolymorphismPosUp[0]][4]
        for comptPol in range(len(listPolymorphismPosUp)):
            if (int(reverseUp)==-1):
                nucleoLow=ReverseComplement(listnucleoLow[comptPol])
            elif int(reverseUp)==1:
                nucleoLow=listnucleoLow[comptPol]
            nucleoRefLow="."
            positionSnpUp=dicoUp[listPolymorphismPosUp[comptPol]][5]
            nucleoUp=dicoUp[listPolymorphismPosUp[comptPol]][3]
            nucleoRefUp=dicoUp[listPolymorphismPosUp[comptPol]][1]
            table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,".",filterField,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,0)
            tablebis.append(list(table))
        tablebis=sorted(tablebis, key=lambda colonnes: colonnes[1])
        l=0
        ID=1
        for l in range(len(tablebis)):
                tablebis[l][2]=str(table[2])+"_"+str(ID)
                printOneline(tablebis[l],VCF)
                ID+=1
#---------------------------------------------------------------------------------------------------------------------------
##Case : Lower path mapped and upper path unmapped            
    elif int(snpUp[3])<=0 and int(snpLow[3])>0:
        reverseLow=dicoLow[listPolymorphismPosLow[0]][4]
        for comptPol in range(len(listPolymorphismPosLow)):
            if (int(reverseLow)==-1):
                nucleoUp=ReverseComplement(listnucleoUp[comptPol])
            elif int(reverseLow)==1:
                nucleoUp=listnucleoUp[comptPol]
            nucleoRefUp="."
            positionSnpLow=dicoLow[listPolymorphismPosLow[comptPol]][5]
            nucleoLow=dicoLow[listPolymorphismPosLow[comptPol]][3]
            nucleoRefLow=dicoLow[listPolymorphismPosLow[comptPol]][1]
            table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,".",filterField,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,0)
            tablebis.append(list(table))
        tablebis=sorted(tablebis, key=lambda colonnes: colonnes[1])
        l=0
        ID=1
        for l in range(len(tablebis)):
                tablebis[l][2]=str(table[2])+"_"+str(ID)
                printOneline(tablebis[l],VCF)
                ID+=1  
            
##############################################################
# geno: [key][0] : genotype (0/0, 0/1 or 1/1)
# geno: [key][1] : likelihood for the 3 possible genotypes ("x,y,z")
##############################################################
def GetGenotype(geno,boolRefLow,table,nbGeno,phased,listCovGeno,cov):
    """Gets the genotype, the coverage and the likelihood by sample and print it in the correspondand fields. The genotype is determined by DiscoSnp++ (which considered the upper path as reference). If the “REF” corresponds the upper path, the genotype in the VCF is identical to the genotype in DiscoSnp++, else  it's the opposite ( 1/1 becomes 0/0 and so on)."""
    j=0
    genotypes=""
    key=None
    current_genotype=None
    likelihood=None
    coverage=None
    #cov="C1=53425,22;C2=1,2;C3=30,94832;C4=3,2;C5=19,65370"
    listcov=cov.split(";") #['C1=53425,22', 'C2=1,2', 'C3=30,94832', 'C4=3,2', 'C5=19,65370']
    if int(nbGeno)==0:
        return table
    else:
            for i in range(0,nbGeno):
                key="G"+str(i+1) # Creates the dictionary key
                current_genotype = geno[key]
                likelihood=current_genotype[1]
                if boolRefLow==True: #Checks if the mapped path is the lower (in this case exchange 0/0 to 1/1 and 1/1 to 0/0 exchanges the likelihood to have the good one for each genotypes)
                    likelihoodStart=likelihood[2]
                    likelihoodEnd=likelihood[0]
                    likelihood[0]=likelihoodStart
                    likelihood[2]=likelihoodEnd
                    if "1/1" in current_genotype[0]:
                        current_genotype[0]=current_genotype[0].replace("1/1","0/0")
                    elif "0/0" in current_genotype[0]:
                        current_genotype[0]=current_genotype[0].replace("0/0","1/1")
                
                if phased==True: #In case of phasing we change the "/" symbol
                    current_genotype[0]=current_genotype[0].replace("/","|")
                #TODO REMOVE THIS TEST WHEN WE HAVE ALL FILES WITH RIGHT HEADER
                if isinstance(likelihood,list):
                    likelihood=str(','.join(current_genotype[1]))
                else:
                    likelihood=str(likelihood)
                coverage=listcov[i].split("=") 
                genotypes+=str(current_genotype[0])+":"+str(listCovGeno[i])+":"+likelihood+":"+str(coverage[1])
                #genotypes+=str(current_genotype[0])+":"+str(listCovGeno[i])+":"+str(','.join(current_genotype[1])) # Add the current genotype
                
                if i<nbGeno-1 :
                    genotypes+="\t" #Adds a \t except if this is the last genotype

            table[8]="GT:DP:PL:AD"
            table[9]=genotypes
    return table
##############################################################
##############################################################
def PrintVCFGhost(table,numSNPUp,chrom,position,tp,valRankUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,cov,ntUp,ntLow,geno,nbGeno,phased,listCovGeno,VCF):
    '''Without samfile some fields will have default value : CHROM table[0],POS table[1], QUAL table[5], FILTER [6]''' 
    table[0]=chrom
    table[1]=position
    table[5]="."
    table[6]="."
    
    table[2]=numSNPUp
    table[3]=ntUp
    table[4]=ntLow
    table[7]="Ty="+str(tp)+";"+"Rk="+str(valRankUp)+";"+"UL="+str(unitigLeftUp)+";"+"UR="+str(unitigRightUp)+";"+"CL="+str(contigLeftUp)+";"+"CR="+str(contigRightUp)
    if "higher" in chrom:
        table=GetGenotype(geno,0,table,nbGeno,phased,listCovGeno,cov)
    elif "lower" in chrom:
        table=GetGenotype(geno,1,table,nbGeno,phased,listCovGeno,cov)
    table[7]=table[7].replace("None",".")
    table[7]=table[7].replace("none",".")
    printOneline(table,VCF)
    
##############################################################
##############################################################
def FillVCF(table,numSNP,chrom,pos,ref,alt,qual,filterfield,tp,valRank,ok,unitigLeft,unitigRight,contigLeft,contigRight,cov,nucleoRef,reverse,geno,nbGeno,phased,listCovGeno,boolRefLow):
    """Take all necessary input variables to fill the vcf;  Fills the fields of the table which will be printed in the vcf ; return the table"""
    if chrom=="*":
        table[0]="."
    else:
        table[0]=chrom
    table[1]=pos
    table[2]=numSNP
    table[3]=ref
    table[4]=alt
    if qual=="*": 
        table[5]="."
    else:
        table[5]=qual
    table[6]=filterfield
    table[7]="Ty="+str(tp)+";"+"Rk="+str(valRank)+";"+"DT="+str(ok)+";"+"UL="+str(unitigLeft)+";"+"UR="+str(unitigRight)+";"+"CL="+str(contigLeft)+";"+"CR="+str(contigRight)+";"+"Genome="+str(nucleoRef)+";"+"Sd="+str(reverse)
    table[7]=table[7].replace("None",".")
    table[7]=table[7].replace("none",".")
    table=GetGenotype(geno,boolRefLow,table,nbGeno,phased,listCovGeno,cov)
    return table  
##############################################################
##############################################################
def ReverseSeq(seq):
    """Take a sequence and reverse it"""  
    i=0
    listSeq=list(seq)
    seq=''
    while i<len(listSeq):
        if seq!='':
                seq=str(ReverseComplement(listSeq[i]))+seq
        else :
                seq=str(ReverseComplement(listSeq[i]))
        i+=1
    return(seq)


##############################################################
##############################################################
def ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleo):
    """Checks strand of the reference (forward or reverse) and reverse or not the nucleotide of the alternative path """
    if (int(reverseUp)==-1 and int(reverseLow)==1) or (int(reverseUp)==1 and int(reverseLow)==-1):
        nucleo=ReverseComplement(nucleo)
    return(nucleo)

##############################################################
#data=str(bin(int(FLAG SAM)))[::-1]
#bit = field to test
#0   read paired
#1   read mapped in proper pair
#2   read unmapped
#3   mate unmapped
#4   read reverse strand
#5   mate reverse strand
#6   first in pair
#7   second in pair
#8   not primary alignment
#9   read fails platform/vendor quality checks
#10  read is PCR or optical duplicate
#11  supplementary alignment
#
##############################################################

def CheckBitwiseFlag(FLAG,bit):
        """Checks if the BitwiseFlag contains the tested value such as : read reverse strand, read unmmaped and so on."""
        try:
                data=str(bin(int(FLAG)))[::-1]
        except ValueError:
                return False
        try:
                return(data[bit]=='1')
        except IndexError:
                return False
##############################################################
##############################################################
def CheckContigUnitig(unitig,contig):
        """Checks if there is an extension in the form of contig or unitig to take into account in the mapping position"""
        if contig:
                return(int(contig))
        elif unitig:
                return(int(unitig))
        else:
                return 0

##############################################################
##############################################################
def GetTag(snp):
        if isinstance(snp,str):#Converts the line of the samfile into a list
                snp=snp.rstrip('\r')
                snp=snp.rstrip('\n')
                snp=snp.split('\t')        
        for field in snp:
                if "MD" in field:
                     MD,Z,posMut = field.split(":") #MD tag parsing
                if "NM" in field:
                     garbage,garbage,nbMismatch=field.split(":")
        return(posMut,nbMismatch)




