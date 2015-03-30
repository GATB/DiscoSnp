#!/bin/python
# -*- coding: utf-8 -*-
###############################################
#v10022015Disco++ : Check des résultats - "compilation" séparée - voir pour mieux trier les fonctions
import os
import sys
import subprocess
import re
import time
#FUNCTIONS_________________________________________________________________________________________________________
##############################################################
##############################################################
def Comptage(fichier):
    """Function that counts the number of SNPs"""
    samfile=open(fichier,'r')
    nbSnp=0
    i=0
    nbGeno=0
    for snp in samfile :
        matchInt=re.match(r'^@',snp)
        if matchInt:
            continue
        else:
            nbSnp+=1
            if isinstance(snp,str):
                snp=snp.rstrip('\n')
                snp=snp.split('\t')
            for i in snp:
                if 'SNP' or 'INDEL' in i:
                    nomDisco=i.split('|')
                    break
            i=0
            for i in nomDisco: # Takes into account the number of variant 
                if "nb_pol" in i :
                    nb_pol=i.split('_')
                    j=0
                    for j in nb_pol:
                        matchInt=re.match(r'^\d+$',j)
                        if matchInt:
                            nbSnp=nbSnp+int(j)-1
            if nbGeno==0: #Count the number of genotypes in the header
                for k in nomDisco:
                    match=re.match(r'^G',k)
                    if match:
                        nbGeno+=1
    return(int(nbSnp)/2,nbGeno)

##############################################################
#Take the new discoSnp++ header : example SNP_higher_path_13736|P_1:30_A/G|high|nb_pol_1|C1_53425|C2_1|C3_30|C4_3|C5_19|G1_0/0:1947,160066,1067676|G2_0/1:40,9,20|G3_1/1:1895477,284398,3575|G4_0/1:34,9,54|G5_1/1:1306661,196101,2493|Q1_66|Q2_61|Q3_36|Q4_55|Q5_33|rank_0.99909
##############################################################
def ParsingDiscoSNP(snp,boolNum):
    if isinstance(snp,str):
        snp=snp.rstrip('\r')
        snp=snp.rstrip('\n')
        snp=snp.split('\t')
    listCoverture=[]
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
    ln=0
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
    i=0
    for i in nomDisco:
        if "SNP" in i:
            discoName=i
            SNP=1
            ID= i.split('_')
            for j in ID:
                matchInt=re.match(r'^\d+$',j)
                if matchInt:
                    numSNP=j
        elif 'P_' in i: ## Esssential
            pos=i.split(',')
            j=0
            ln=0
            ntUp=''
            ntLow=''
            key=''
            chaine=[]
            ind=0
            dicoHeader={}
            for j in pos:
                posD=0
                if SNP==1:
                    if ":" in j:
                        chaine=j.split(":")
                        k=0
                        for k in chaine:
                            if "P_" in k:
                                key=k
                            elif "/" in k:
                                chaine1=[]
                                chaine1=k.split('_')
                                l=0
                                for l in chaine1:
                                    if "/" in l:
                                        chaine2=[]
                                        chaine2=l.split('/')
                                        ntUp=chaine2[0]
                                        ntLow=chaine2[1]
                                    else:
                                        posD=int(l)+1
                    dicoHeader[key]=[posD,ntUp,ntLow] ##### !!! Essential : In case of snps 
                else:
                    if ":" in j:
                        ind=0
                        posD=0
                        chaine=j.split(":")
                        k=0
                        amb=0
                        for k in chaine:
                            if "P_" in k:
                                key=k
                            else:
                                chaine=k.split('_')
                                l=0
                                for l in chaine:
                                    if posD==0:
                                        posD=int(l)+1
                                    elif ind==0:
                                        ind=int(l)
                                    elif amb==0:
                                        amb=int(l)
                    dicoHeader[key]=[posD,ind,amb] ##### !!! Essential In case of indel 
        elif "INDEL" in i:
            INDEL=True
            ID= i.split('_')
            j=0
            discoName=i
            for j in ID:
                matchInt=re.match(r'^\d+|$',j)
                if matchInt:
                    numSNP=j
        
        elif "unitig" in i:
            if "left" in i:
                unitig=i.split('_')
                j=0
                for j in unitig:
                    matchInt=re.match(r'^\d+$',j)
                    if matchInt:
                        unitigLeft=j
            elif "right" in i :
                unitig=i.split('_')
                j=0
                for j in unitig:
                    matchInt=re.match(r'^\d+$',j)
                    if matchInt:
                        unitigRight=j
        elif "contig" in i:
            if "left" in i:
                contig=i.split('_')
                j=0
                for j in contig:
                    matchInt=re.match(r'^\d+$',j)
                    if matchInt:
                        contigLeft=j
            elif "right" in i :
                contig=i.split('_')
                j=0
                for j in contig:
                    matchInt=re.match(r'^\d+$',j)
                    if matchInt:
                        contigRight=j
        elif "rank" in i: 
            rank=i.split('_')
            j=0
            for j in rank:
                matchInt=re.match(r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?',j)
                if matchInt:
                    valRank=j
        elif "C" in i: ### Essential
            coverage=i.split('_')
            j=0
            for j in range(len(coverage)):
                matchInt=re.match(r'^\d+$',coverage[j])
                if matchInt:
                    listCoverture.append(str(coverage[j]))
                    listC.append(str(coverage[j-1]))
        elif "nb_pol" in i : ###Esssential
            nb_pol=i.split('_')
            j=0
            for j in nb_pol:
                matchInt=re.match(r'^\d+$',j)
                if matchInt:
                    nb_pol=j
        elif "G" in i:
            gen=i.replace("_",":")
            listgeno=gen.split(":")
            if len(listgeno)>2:
                listlikelihood=listgeno[2].split(",")
            else:
                listlikelihood=0
            dicoGeno[listgeno[0]]=[listgeno[1],listlikelihood]
    
    if boolNum==0:
        return(discoName,snp,numSNP,unitigLeft,unitigRight,contigLeft,contigRight,valRank,listCoverture,listC,nb_pol,ln,posD,ntUp,ntLow,dicoGeno,dicoHeader)
    else:
        return(numSNP)

##############################################################
# snpUp: sam line 1
# snpLow: sam line 2
# 
##############################################################
def GetCouple(snpUp,snpLow):
    """Retrieves for each path alignment information in a list ; retrieves a dictionary with all the positions of a path and the number of associated mismatch ; retrieves a Boolean indicating whether the SNP is multimapped or not"""
    posUp = {} #dictionnary of all the position (keys) associated with the number of mismatches (values) 
    posLow = {}
    i=0
    j=0
    #Boolean snps mulimapped
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
    if 'XA' in ''.join(snpUp): # XA: tag for multiple mapping : Check if the upper path is multiple mapped 
        i=0
        #parsing XA tag
        listXA=snpUp[19].split(';')
        strXA = ','.join(listXA)
        listXA = strXA.split(',')
        listXA.pop()
        position=listXA[1:] #position=[pos,cigarcode,number of mismatch , pos,cigarcode,number of mismatch]
        while i<len(position): # runs through the list 4 by 4 to get all the positions 
            if abs(int(position[i])) not in listerreurUp : #Checks if the position is not too close to the main one
                boolXAUp=True
                posUp[abs(int(position[i]))]=int(position[i+2]) #the position is associated to the number of mismatch in a dictionnary
            i+=4
    if 'XA' in ''.join(snpLow):# XA: tag for multiple mapping : Check if the lower path is multiple mapped 
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
    #Add first position give by BWA to the dictionnary if the path is mapped
    if abs(int(snpUp[3]))>0:
        garbage,garbage,nbMismatchUp=snpUp[12].split(":")
        posUp[abs(int(snpUp[3]))]=int(nbMismatchUp)
    if abs(int(snpLow[3]))>0:
        garbage,garbage,nbMismatchLow=snpLow[12].split(":")
        posLow[abs(int(snpLow[3]))]=int(nbMismatchLow)
    return(posUp,posLow,boolXAUp,boolXALow)
    ## posUp/posLow : dictionnary with all the position associated with their mismatch number
    ## boolXAUp/boolXALow : boolean if TRUE : the path is multiple mapped for BWA  and the max mapping distance

##############################################################
#listCup=["C1","C2"]
#listCovertureUp=[23,45]
##############################################################
def GetCoverage( listCUp, listCLow, listCoverageUp, listCoverageLow):
    """Takes as input lists list  """
    i=0
    covUp=''
    listCovGeno=[]
    covLow=''
    ##Create a string if the number of coverage is superior to 1
    if len( listCoverageUp)>1 and len( listCoverageLow)>1:
        while i<len( listCoverageUp):
            listCovGeno.append(int( listCoverageUp[i])+int( listCoverageLow[i]))
            if covUp!='':
                covUp=str(covUp)+';'+str( listCUp[i])+'='+str( listCoverageUp[i])+","+str( listCoverageLow[i])
            else:
                covUp=str( listCUp[i])+'='+str( listCoverageUp[i])+","+str( listCoverageLow[i])
            i+=1
        i=0
        covLow=''
        while i<len( listCoverageLow):
            if covLow!='':
                covLow=str(covLow)+';'+str( listCLow[i])+'='+str( listCoverageLow[i])+","+str( listCoverageUp[i])
            else:
                covLow=str( listCLow[i])+'='+str( listCoverageLow[i])+","+str( listCoverageUp[i])
            i+=1
    else:
        listCovGeno.append(int(listCoverageUp[0])+int(listCoverageLow[0]))
        covUp=str( listCUp[0])+'='+str( listCoverageUp[0])+","+str( listCoverageLow[0])
        covLow=str( listCLow[0])+'='+str( listCoverageLow[0])+","+str( listCoverageUp[0])
    return(covUp,covLow,listCovGeno) #string covUp C1=5,23;C2=35,1 listCovGeno=[28,36]
##############################################################
#position : current position to add at the ensemble
#delta : minimum number of difference allowed between two positions
##############################################################
def addPosition(position,ensemble,delta):
    """Add a position to a set of positions only if the difference between the two is greater than delta"""
    if len(ensemble)==0: # Case : it's the first position add to the set
        ensemble.append(position)
        return(ensemble)
    if abs(position-ensemble[0])>delta:
        ensemble.append(position)
        return(ensemble)
    return(ensemble)
    #ensemble : set of positions for a path
##############################################################
# snpUp/snpLow : sam line
# posUp/posLow : dictionnary with all the position associated with their mismatch number example : posUp[5687884684]=2
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
            if nbMismatch==int(NM):
                listPos=addPosition(position,listPos,delta) #If the position is mapped at NM (number of mismatch tested) test if its not too close with an other position
                ensemble=set(listPos)
                if abs(len(ensemble))>1: #case of many position of mapping for NM => multiple mapped  Exit function
                    couple="multiple"
                    return(couple)
        for position,nbMismatch in posLow.items():# Checks for the lower path
            if nbMismatch==int(NM):
                listPos=addPosition(position,listPos,delta)
                ensemble=set(listPos)
                if abs(len(ensemble))>1:#case of many position of mapping for NM => multiple mapped  Exit function
                    couple="multiple"
                    return(couple)
        if ensemble!=None: # In case of an ensemble==1 uniq genomic mapping position 
            couple="ok"
            return(couple)
    return(couple)
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
# posModif usefull in case of simple snps or indel : give the position of the variant on the sequence
# indel : boolean True if the current variant is an indel
##############################################################
def CigarCodeChecker(cigarcode,listpol,posModif,indel):
    """Function which allows to recover the position of the SNPs according to the position of the deletions or insertions relative to the reference sequence (CigarCode Parsing : checks for insertion deletion or soft clipping"""
    parsingCigarCode=re.findall('(\d+|[A-Za-z])',cigarcode) #parsingCigarCode=['2', 'S', '3', 'M', '1', 'I', '25', 'M']
    listPosRef=[]
    listShift=[]
    somme=0
    shift=0
    pos=0
    i=1
    j=0
    lenDemiSeq=posModif
    posCentraleRef=None
    #Close snps
    if len(listpol)>1:
        while i<len(parsingCigarCode): # Goes through the list by twos to get all the letters and to take them into account
            #Soft clipping
            if parsingCigarCode[i]=="S":
                shift-=int(parsingCigarCode[i-1]) # It's the shift in the alignment between the reference and the sequence of the variant 
                pos+=int(parsingCigarCode[i-1]) # pos corresponds to the current position in the cigarcode
            #Match or Mismatch  
            elif parsingCigarCode[i]=="M":
                pos+=int(parsingCigarCode[i-1])
            #Deletion
            elif parsingCigarCode[i]=="D":
                shift+=int(parsingCigarCode[i-1])
            #Insertion
            elif parsingCigarCode[i]=="I":
                shift-=int(parsingCigarCode[i-1])
            # hard clipping (clipped sequences NOT present in SEQ)
            elif parsingCigarCode[i]=="H":
                shift-=int(parsingCigarCode[i-1]) # It's the shift in the alignment between the reference and the sequence of the variant 
                pos+=int(parsingCigarCode[i-1])
            # padding (silent deletion from padded reference)
            elif parsingCigarCode[i]=="P":
                shift+=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="=":
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="X":
                pos+=int(parsingCigarCode[i-1])
            while int(pos)>=int(listpol[j]):# Goes through the list of position (close snps) and see if the position is affected by the shift (means shift before the position)
                posRef=int(listpol[j])+shift # Add the shift to the position (to get the real position of the snp on the reference)
                listPosRef.append(posRef) #Add the position to the list by taking into account the shift only if the current position 
                listShift.append(shift)
                if j<(len(listpol)-1):
                    j+=1
                else:
                    #if listShift==[]:
                    #    listShift=[0]*len(listpol)
                    return(listPosRef,listShift)
            i+=2
    #Simple snp and Indel
    else:
        while i<len(parsingCigarCode):# Goes through the list by twos to get all the letters and to take them into account
            if parsingCigarCode[i]=="S":
                shift-=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="M":
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="D":
                shift+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="I":
                shift-=int(parsingCigarCode[i-1])
            # hard clipping (clipped sequences NOT present in SEQ)
            elif parsingCigarCode[i]=="H":
                shift-=int(parsingCigarCode[i-1]) # It's the shift in the alignment between the reference and the sequence of the variant 
                pos+=int(parsingCigarCode[i-1])
            # padding (silent deletion from padded reference)
            elif parsingCigarCode[i]=="P":
                shift+=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="=":
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="X":
                pos+=int(parsingCigarCode[i-1])
            if len(listpol)==1:
                if pos>=int(lenDemiSeq):
                    posCentraleRef=int(lenDemiSeq)+shift#takes into account the shift to add the after the mapping position and the variant position in the sequence
                    return(posCentraleRef,shift)
            i+=2
##############################################################
#shift : integer shift in the alignment between the reference and the sequence of the variant 
#posMut : column 19 in the samfile example 30T30 means 30 matches between the ref and the seq of the variant, one mismatch (T on the reference) and 30 matches between the ref and the seq of the variant
#posCentraleRef : position of the variant on the reference (taking the shift into account)
##############################################################
def ReferenceChecker(shift,posMut,posCentraleRef):
    """MD tag parsing; checks if path nucleotide is identical to the reference nucleotide """
    pos=0
    i=0
    matchInt=None
    nucleoRef=None
    boolEgalRef=None
    motifDel=None
    dictDel={}
    deletion=None
    numeroDel=None
    parsingPosMut=None
    strNumerodel=None
    if '^' in posMut: #Means there is a deletion relative to the reference
        pos=shift
        motifDel=re.compile("[0-9]*\^[A-Za-z]*") #motif of the deletion example  : 31^CGC
        dictDel={}
        deletion=re.findall('\d+\^[A-Za-z]*',posMut) # Parsing with the deletion
        deletion=''.join(deletion)
        numeroDel=re.findall('\d+',deletion)
        parsingPosMut=re.findall('(\d+|[A-Za-z]|\^)',posMut)
        strNumerodel=''.join(numeroDel)
        for m in motifDel.finditer(posMut):
            dictDel[m.group()]=parsingPosMut.index('^')-1 #Get the position of the deletion in posMut
        posMut=posMut.replace(deletion,"")
        parsingPosMut=re.findall('(\d+|[A-Za-z])',posMut)
        parsingPosMut.insert(dictDel[deletion],numeroDel[0])
    else:
        parsingPosMut=re.findall('(\d+|[A-Za-z])',posMut)
    while i<len(parsingPosMut):#Convert the number into integer 
        matchInt=re.match(r'\d+',parsingPosMut[i])
        if matchInt:
            parsingPosMut[i]=int(parsingPosMut[i])
        i+=1
    i=0
    while i<len(parsingPosMut): # Goes to the list to know if the variant is identical to the reference
        if isinstance(parsingPosMut[i],int): #Checks if it's a number of nucleotide or a letter
            pos+=parsingPosMut[i] # Add to pos the current number of the list
        else :
            pos+=1 #else add the letter to the pos
        if pos==posCentraleRef: # Checks if the current position pos is identical to the position of the variant 
            if isinstance(parsingPosMut[i],str): #=> it means that the nucleotide is differente in the variant and in the reference
                boolEgalRef=False
                nucleoRef=parsingPosMut[i]
            else:
                boolEgalRef=True
            break
        if pos>posCentraleRef:
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
    """Get reverse sequence"""
    seqLow=snpLow[9]
    seqUp=snpUp[9]
    listSeqUp=[]
    listSeqLow=[]
    i=0
    if CheckBitwiseFlag(snpUp[1],4): #Case of mapping reverse
        i=0
        listSeqUp=list(seqUp)
        seqUp=''
        while i<len(listSeqUp): #Reverse the sequence usefull for future comparison
            if seqUp!='':
                seqUp=str(ReverseComplement(listSeqUp[i]))+seqUp
            else :
                seqUp=str(ReverseComplement(listSeqUp[i]))
            i+=1
    if CheckBitwiseFlag(snpLow[1],4):#Case of mapping reverse
        i=0
        listSeqLow=list(str(seqLow))
        seqLow=''
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
#posUp/posLow : dictionnary with all the mapping position of every snps associated with their number of mismatch
#nb_polUp/nb_polLow : number of variants for the current path
#indel : boolean True if the current variant is an indel
##############################################################
def RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel):
    #INIT VALUES
    posSNPUp="0"+str(dicoHeaderUp["P_1"][0]) # in case of unmapped variant example position : 031
    posSNPLow="0"+str(dicoHeaderUp["P_1"][0])
    boolRefUp=None
    boolRefLow=None
    nucleoUp=None
    nucleoLow=None
    nucleoRefUp=None
    nucleoRefLow=None
    listPolymorphismePosUp=None
    listPolymorphismePosLow=None
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
    
    posModif=dicoHeaderUp["P_1"][0] #Position of the first variant gived by discosnp
    seqUp=list(snpUp[9])
    seqLow=list(snpLow[9])
    if posModif==0:
        posModif=len(seqUp)/2
    reverseUp=1
    reverseLow=1
    dicopolUp={}
    dicopolLow={}
    i=0
    if int(snpUp[3])>0:
        posRef=int(snpUp[3])
    else:
        posRef=int(snpLow[3])
#---------------------------------------------------------------------------------------------------------------------------
    #List of snps
    if indel==False:
        listPos,listnucleoUp,listnucleoLow,listPosR,listnucleoUpR,listnucleoLowR=GetPolymorphisme(dicoHeaderUp,snpUp[9],indel)
    else:
        seqUp=snpUp[9]
        seqLow=snpLow[9]
        if len(seqUp)<len(seqLow): #Determines the longest sequence to get the insertion
            seq=seqLow
        else:
            seq=seqUp
        listPos,listPosR,insert,ntStart,ambiguityPos=GetPolymorphisme(dicoHeaderUp,seq,indel) # For indel get the insert, the list of position and the possible ambiguity for the position 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
###Case mapped variant Up : Shift by positions (insertion,deletion,sofclipping) and update of the position in alignment
    if int(snpUp[3])>0:
        #Check cigarCode : Presence of insertion, softclipping, deletion
        if not CheckBitwiseFlag(snpUp[1],4) : #Forward Strand : check the bitwise flag of the samfile
            posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],listPos,posModif,indel) # Gets the positions of the variants with an eventual shift from the reference (insertion,deletion,soft clipping) ; in case of close snps return a list 
            listPolymorphismePosUp=listPos
            if len(listPos)==1 and indel==False: #simple snp
                nucleoUp=listnucleoUp[0] # Gets the nucleotide 
        elif CheckBitwiseFlag(snpUp[1],4):# Reverse Strand
            reverseUp=-1
            posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],listPosR,posModif,indel)
            listPolymorphismePosUp=listPosR #List of all the reverse position
            if len(listPos)==1 and indel==False:
                nucleoUp=listnucleoUpR[0]
    else: ###Case unmapped variant Up
        listPolymorphismePosUp=listPos
        posCentraleUp=listPos
        shiftUp=None
        reverseUp="."
###Case mapped variant Low :Shift by positions (insertion,deletion,sofclipping) and update of the position in alignment
    if int(snpLow[3])>0:
        #Check cigarCode : Presence of insertion, softclipping, deletion
        if not CheckBitwiseFlag(snpLow[1],4):#Forward Strand
            posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],listPos,posModif,indel)
            listPolymorphismePosLow=listPos
            if len(listPos)==1 and indel==False:
                nucleoLow=listnucleoLow[0]
        elif CheckBitwiseFlag(snpLow[1],4):# Reverse Strand
            reverseLow=-1
            posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],listPosR,posModif,indel)
            listPolymorphismePosLow=listPosR
            if len(listPos)==1 and indel==False:
                nucleoLow=listnucleoLowR[0]
    else:###Case unmapped variant Low
        listPolymorphismePosLow=listPos
        posCentraleLow=listPos
        shiftLow=None
        reverseLow="."
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
    ####Upper Case : Check for all the variants if there are identical to the reference and proc 
    if int(snpUp[3])>0:
        MD,Z,posMutUp = snpUp[18].split(":") # MD tag parsing
        if len(listPos)>1: # CLOSE SNPS
            if reverseUp==1:
                i=0
                for i in range(len(listPos)):
                    boolRefUp,nucleoRefUp=ReferenceChecker(shiftUp[i],posMutUp,posCentraleUp[i]) # Check if the variant is identical to the reference ; return a boolean and the nucleotide of the reference
                    if nucleoRefUp==None: #If there is no reference nucleotide gived by ReferenceChecker, it means that the variant is equal to the reference so we defined it !
                        nucleoRefUp=listnucleoUp[i] 
                    dicopolUp[listPos[i]]=[boolRefUp,nucleoRefUp,posCentraleUp[i],listnucleoUp[i],reverseUp,(int(snpUp[3])+posCentraleUp[i])] # Dictionnary for close snps to keep all the results gived by the different functions
            elif reverseUp==-1:
                i=0
                for i in range(len(listPosR)):
                    boolRefUp,nucleoRefUp=ReferenceChecker(shiftUp[i],posMutUp,posCentraleUp[i])
                    if nucleoRefUp==None: 
                        nucleoRefUp=listnucleoUpR[i]
                    dicopolUp[listPosR[i]]=[boolRefUp,nucleoRefUp,posCentraleUp[i],listnucleoUpR[i],reverseUp,(int(snpUp[3])+posCentraleUp[i])]
        else: #SIMPLE SNPS
            boolRefUp,nucleoRefUp=ReferenceChecker(shiftUp,posMutUp,posCentraleUp)
            if nucleoRefUp==None and reverseUp==1:#If there is no reference nucleotide gived by ReferenceChecker, it means that the variant is equal to the reference so we defined it !
                nucleoRefUp=dicoHeaderUp["P_1"][1]
            elif nucleoRefUp==None and reverseUp==-1:
                nucleoRefUp=ReverseComplement(dicoHeaderUp["P_1"][1])
	    
    
    elif int(snpUp[3])<=0 :
        nucleoUp = dicoHeaderUp["P_1"][1]
        if reverseLow==-1:
               nucleoUp = ReverseComplement(dicoHeaderUp["P_1"][1])
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
    ####Lower Case 
    if int(snpLow[3])>0:
        #Extraction tag MD
        MD,Z,posMutLow = snpLow[18].split(":")
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
            if nucleoRefLow==None and reverseLow==1:#If there is no reference nucleotide gived by ReferenceChecker, it means that the variant is equal to the reference so we defined it !
                nucleoRefLow=dicoHeaderUp["P_1"][2]
            elif nucleoRefLow==None and reverseLow==-1:
                nucleoRefLow=ReverseComplement(dicoHeaderUp["P_1"][2])
    
    elif int(snpLow[3])<=0:
        nucleoLow = dicoHeaderUp["P_1"][2]
        if reverseUp==-1:
               nucleoLow = ReverseComplement(dicoHeaderUp["P_1"][2])

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
####Two paths mapped at different main position (first given by bwa) : test which path will be the reference with the function mismatchChecker
    if (int(snpUp[3])>0 and int(snpLow[3])>0) and (int(snpUp[3])!=int(snpLow[3])):
        if len(listPos)>1:
            i=0
            for i in range(len(listPos)):
                boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,dicopolUp[listPolymorphismePosUp[i]][1],dicopolLow[listPolymorphismePosLow[i]][1],dicopolUp[listPolymorphismePosUp[0]][3],dicopolLow[listPolymorphismePosLow[0]][3],dicopolUp[listPolymorphismePosUp[i]][0],dicopolLow[listPolymorphismePosLow[i]][0],indel)
                dicopolUp[listPolymorphismePosUp[i]][0]=boolRefUp
                dicopolUp[listPolymorphismePosUp[i]][1]=nucleoRefUp
                dicopolLow[listPolymorphismePosLow[i]][0]=boolRefLow
                dicopolLow[listPolymorphismePosLow[i]][1]=nucleoRefLow
                dicopolLow[listPolymorphismePosLow[i]][5]=int(dicopolLow[listPolymorphismePosLow[i]][2])+int(posRef)
                dicopolUp[listPolymorphismePosUp[i]][5]=int(dicopolUp[listPolymorphismePosUp[i]][2])+int(posRef)
        else:
                boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel)
#---------------------------------------------------------------------------------------------------------------------------    #---------------------------------------------------------------------------------------------------------------------------
    #Variants (INDEL and SIMPLE SNPS) Positions : add the position of mapping gived by bwa to the variant position and for indel add ambiguity position
    if int(nb_polUp)==1:
        if indel==True:
            nucleoLow="."
            nucleoUp="."
            nucleoRefUp="."
            nucleoRefLow="."
        if (int(snpUp[3])>0 and int(snpLow[3])>0):
            if indel==True:
                posSNPUp=posCentraleUp+posRef-int(ambiguityPos)
                posSNPLow=posCentraleLow+posRef-int(ambiguityPos)
                if snpUp[9]<snpLow[9]:
                    boolRefUp=True
                    boolRefLow=False
                else:
                    boolRefUp=False
                    boolRefLow=True
            else:
                posSNPUp=posCentraleUp+posRef
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
        elif int(snpUp[3])>0 and int(snpLow[3])<=0:
            if indel==True:
                boolRefUp=True
                boolRefLow=False
                posSNPLow=None
                posSNPUp=posCentraleUp+posRef-int(ambiguityPos)
            else:
                posSNPLow=None
                posSNPUp=posCentraleUp+posRef
        elif (int(snpUp[3])<=0 and int(snpLow[3])<=0):
                posSNPUp=dicoHeaderUp["P_1"][0]
                posSNPLow=dicoHeaderUp["P_1"][0]
        ##Check the strand : if the alternative nucleotide is not on the same strand as the reference : reverse it
        if boolRefLow==True and reverseLow==-1 and reverseUp==1:
               nucleoUp=ReverseComplement(nucleoUp)
        elif boolRefLow==True and reverseLow==1 and reverseUp==-1:
               nucleoUp=ReverseComplement(nucleoUp)
        if boolRefUp==True and reverseUp==-1 and reverseLow==1:
               nucleoLow=ReverseComplement(nucleoLow)
        elif boolRefUp==True and reverseUp==1 and reverseLow==-1:
               nucleoLow=ReverseComplement(nucleoLow)
        ##Check which path chooses if there are both different from the reference
        if boolRefUp==False and boolRefLow==False:
               boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel)
        return(nucleoLow,posSNPLow,nucleoUp,posSNPUp,boolRefLow,boolRefUp,reverseUp,reverseLow,nucleoRefUp,nucleoRefLow)
        #nucleoLow/nucleoUp : nucleotide of each path corresponding to the variant
        #posSNPLow/posSNPUp : position of the variant taking into account the shift
        #boolRefUp/boolRefLow : boolean TRUE the variant (of the path Up or Low) is the reference
    else:
        return(dicopolUp,dicopolLow,listPolymorphismePosUp,listPolymorphismePosLow)
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
    if int(snpUp[3])>0 and int(snpLow[3])>0:
        nmUp=posUp[int(snpUp[3])]
        nmLow=posLow[int(snpLow[3])]
        if nmUp<nmLow:
            boolRefLow=False
            boolRefUp=True
            nucleoRefLow=nucleoRefUp
            posRef=int(snpUp[3])
        elif nmUp>nmLow :
            boolRefLow=True
            boolRefUp=False
            nucleoRefUp=nucleoRefLow
            posRef=int(snpLow[3])
        elif  nmUp==nmLow:
            if indel==False:
                if nucleoUp<nucleoLow:
                    boolRefLow=False
                    boolRefUp=True
                    nucleoRefLow=nucleoRefUp
                    posRef=int(snpUp[3])
                elif nucleoUp>nucleoLow:
                    boolRefLow=True
                    boolRefUp=False
                    nucleoRefUp=nucleoRefLow
                    posRef=int(snpLow[3])
                else :
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
            
            else:
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
    return(boolRefLow,boolRefUp,nucleoRefUp,nucleoRefLow,posRef)
##############################################################
#Simple Variant : dicoHeader[key]=[posD,ntUp,ntLow]
#Indel : dicoHeader[key]=[posD,ind,amb]
##############################################################
def GetPolymorphisme(dicoHeader,seq,indel):
    '''Gets from the dicoHeader all the positions, and the nucleotides (R means that it's on the reverse strand)
      SNP :  one variant correspond to listPos[0], listPosR[0], listnucleoUp[0], listnucleoLow[0],listnucleoUpR[0],listnucleoLowR[0]'''
    #Forward
    listPos=[]
    listnucleoUp=[]
    listnucleoLow=[]
    #Reverse
    listPosR=[]
    listnucleoUpR=[]
    listnucleoLowR=[]
    tailleSeq=len(seq)
    if indel==False:##Case of simple snp
        for key,(posD,ntUp,ntLow) in dicoHeader.items():
            listPos.append(posD)
            listPosR.append(tailleSeq-int(posD)+1)
            listnucleoUp.append(ntUp)
            listnucleoUpR.append(ReverseComplement(ntUp))
            listnucleoLow.append(ntLow)
            listnucleoLowR.append(ReverseComplement(ntLow))
        return(listPos,listnucleoUp,listnucleoLow,listPosR,listnucleoUpR,listnucleoLowR)
    else:##Case of indel
        for key,(posD,ind,amb) in dicoHeader.items():
            listPos.append(posD)
            listPosR.append(tailleSeq-int(posD)+1)
            insert=seq[(int(posD-1)-1):(int(posD-1)+int(ind))]
            ntStart=seq[(int(posD-1)-1)]
            ambiguityPos=amb
        return(listPos,listPosR,insert,ntStart,ambiguityPos)
##############################################################
##############################################################
def fillVCFSimpleSnp(snpUp,snpLow,nucleoLow,positionSnpLow,nucleoUp,positionSnpUp,boolRefLow,boolRefUp,table,nbSnp,dmax,filterfield,multi,ok,tp,phased,listCovGeno,nucleoRefUp,nucleoRefLow,reverseUp,reverseLow,geno,nbGeno,covUp,covLow):
    """ Fills the different fields of vcf based on boolean ( on whether the SNP is identical or not the reference) , if neither is identical to the reference = > we take the SNP comes first in the lexicographical order"""
    ##Gets the variable of the header of disco snps
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
    posUnmappedUp=None
    discoNameUp,snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp, listCoverageUp, listCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
    discoNameLow,snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow, listCoverageLow,listClow,nb_polLow,lnlow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
    posUnmappedUp=CheckContigUnitig(unitigLeftUp,contigLeftUp)
    posUnmappedLow=CheckContigUnitig(unitigLeftLow,contigLeftLow)
#---------------------------------------------------------------------------------------------------------------------------
##Case : two mapped paths
    if int(snpUp[3])>0 and int(snpLow[3])>0:
        ##The path identical to the reference is the lower path 
        if boolRefLow==True and boolRefUp==False:
            table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterfield,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
            
         ##The path identical to the reference is the upper path 
        elif boolRefUp==True and boolRefLow==False:
            table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterfield,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
        ##No path is identical to the reference => lexicographique choice
        elif boolRefUp==False and boolRefLow==False:
            if nucleoUp<nucleoLow:
                table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterfield,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
            elif nucleoUp>nucleoLow:
              table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterfield,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
              
#---------------------------------------------------------------------------------------------------------------------------
##Case : Lower path mapped and upper path unmapped      
    elif int(snpUp[3])<=0 and int(snpLow[3])>0:
        table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterfield,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
        if snpLow[10]=="*":
                table[5]="."
        else:
                table[5]=snpLow[10]  
#---------------------------------------------------------------------------------------------------------------------------
##Case : Upper path mapped and lower path unmapped        
    elif int(snpUp[3])>0 and int(snpLow[3])<=0:
        table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterfield,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
        if snpUp[10]=="*":
                table[5]="."
        else:
                table[5]=snpUp[10]
#---------------------------------------------------------------------------------------------------------------------------
##Case : Both paths are unmapped       
    elif int(snpUp[3])<=0 and int(snpLow[3])<=0:#FillVCF(table,numSNP,chrom,pos,ref,alt,qual,filterfield,tp,valRank,multi,ok,unitigLeft,unitigRight,contigLeft,contigRight,cov,nucleoRef,reverse,geno,nbGeno,phased,listCovGeno,boolRefLow)
        if nucleoLow<nucleoUp:
            table=FillVCF(table,numSNPLow,discoNameUp,(int(positionSnpUp)+int(posUnmappedUp)),nucleoLow,nucleoUp,snpLow[10],filterfield,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,".",".",geno,nbGeno,phased,listCovGeno,boolRefLow)
        else:
            table=FillVCF(table,numSNPUp,discoNameLow,(int(positionSnpLow)+int(posUnmappedLow)),nucleoUp,nucleoLow,snpUp[10],filterfield,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,".",".",geno,nbGeno,phased,listCovGeno,boolRefLow)
    return(table)
##############################################################
##############################################################   
def printOneline(table,VCF):
    if table==[0, 0, 0, 0, 0, 0, 0, 0, 0, 0] or (table[0]==0 and table[1]==0 and table[3]==0 and table[4]==0):
        return
    for i in range(len(table)):
        element=table[i]
        VCF.write(str(element))
        if i<len(table)-1: VCF.write("\t")
    VCF.write('\n')    

##############################################################
#dicopolLow[listPos[i]]=[boolRefLow,nucleoRefLow,posCentraleLow[i],listnucleoLowR[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
##############################################################
def printVCFSNPclose(dicoUp,dicoLow,table,filterField,dmax,snpUp,snpLow,listPolymorphismePosUp,listPolymorphismePosLow,listPolymorphismePos,multi,ok,covUp,covLow,listnucleoUp,listnucleoLow,geno,nbGeno,listCovGeno, VCF):
    info=''
    champAlt=0
    comptPol=0
    seqUp=snpUp[9]
    seqLow=snpLow[9]
    tp="SNP"
    phased=True
    table = [0] * 10 # create a 10 cols array
    ##Keep the close snps to sort them : indeed all the lists and dictionnaries :listnucleoUp,listPolymorphismePosUp,listPolymorphismePosLow,listnucleoLow dicoUp,dicoLow are classified according to dicoHeader so if we start by sorting we lose the correspondence between data
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
    discoNameUp,snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp, listCoverageUp, listCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
    discoNameLow,snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow, listCoverageLow, listCLow,nb_polLow,lnLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
    posUnmappedUp=CheckContigUnitig(unitigLeftUp,contigLeftUp)
    posUnmappedLow=CheckContigUnitig(unitigLeftLow,contigLeftLow)
#---------------------------------------------------------------------------------------------------------------------------
##Case : two mapped paths
    if int(snpUp[3])>0 and int(snpLow[3])>0:
        ##Sorts the list of position to get the smallest one and its position in the list unsorted!
        listSortedPosUp=list(listPolymorphismePosUp)
        listSortedPosUp.sort()
        indexSmallestPosUp=listPolymorphismePosUp.index(listSortedPosUp[0])
        listSortedPosLow=list(listPolymorphismePosLow)
        listSortedPosLow.sort()
        indexSmallestPosLow=listPolymorphismePosLow.index(listSortedPosLow[0])
        boolRefUp=dicoUp[listPolymorphismePosUp[indexSmallestPosUp]][0]
        boolRefLow=dicoLow[listPolymorphismePosLow[indexSmallestPosLow]][0]
        #Decide what is the smallest position according to the reference path (useful if the paths are not aligned on the same strand)
        if boolRefUp==True:
            indexSmallestPos=indexSmallestPosUp
        elif boolRefLow==True:
            indexSmallestPos=indexSmallestPosLow
        elif boolRefUp==False and boolRefLow==False:
            nucleoUp1=dicoUp[listPolymorphismePosUp[indexSmallestPosUp]][3]
            nucleoLow1=dicoLow[listPolymorphismePosLow[indexSmallestPosLow]][3]
            if nucleoUp1<nucleoLow1:
                indexSmallestPos=indexSmallestPosUp
            else:
                indexSmallestPos=indexSmallestPosLow
        #Remembers the values for the first snps ==> the others snps will dependent on these parameters
        boolRefUp=dicoUp[listPolymorphismePosUp[indexSmallestPos]][0]
        boolRefLow=dicoLow[listPolymorphismePosLow[indexSmallestPos]][0]
        nucleoUp1=dicoUp[listPolymorphismePosUp[indexSmallestPos]][3]
        nucleoLow1=dicoLow[listPolymorphismePosLow[indexSmallestPos]][3]
        positionSnpUp1=dicoUp[listPolymorphismePosUp[indexSmallestPos]][5]
        positionSnpLow1=dicoLow[listPolymorphismePosLow[indexSmallestPos]][5]
        reverseLow=dicoLow[listPolymorphismePosLow[indexSmallestPos]][4]
        reverseUp=dicoUp[listPolymorphismePosUp[indexSmallestPos]][4]
        for comptPol in range(len(listPolymorphismePos)):
            positionSnpUp=dicoUp[listPolymorphismePosUp[comptPol]][5]
            positionSnpLow=dicoLow[listPolymorphismePosLow[comptPol]][5]
            nucleoUp=dicoUp[listPolymorphismePosUp[comptPol]][3]
            nucleoRefUp=dicoUp[listPolymorphismePosUp[comptPol]][1]
            nucleoLow=dicoLow[listPolymorphismePosLow[comptPol]][3]
            nucleoRefLow=dicoLow[listPolymorphismePosLow[comptPol]][1]
            #Remplissage VCF
            if boolRefLow==True and boolRefUp==False:
                 table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterField,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
                 table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoUp)
            elif boolRefUp==True and boolRefLow==False:
                table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
                table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoLow)
            elif boolRefUp==False and boolRefLow==False:
                if nucleoUp1<nucleoLow1:
                    table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
                    table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoLow)
                elif  nucleoUp1>nucleoLow1:
                      table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterField,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
                      table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoUp)
                elif positionSnpUp1<positionSnpLow1:
                    table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
                    table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoLow)
                elif positionSnpUp1>positionSnpLow1:
                    table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterField,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
                    table[4]=ReverseCheckerCloseSNP(reverseUp,reverseLow,nucleoUp)
            tablebis.append(list(table))
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
        for i in range(len(listPolymorphismePos)):
            positionSnpUp=int(listPolymorphismePos[i])+int(posUnmappedUp)
            positionSnpLow=int(listPolymorphismePos[i])+int(posUnmappedLow)
            nucleoUp=listnucleoUp[i]
            nucleoRefUp="."
            nucleoLow=listnucleoLow[i]
            nucleoRefLow="."
            table=FillVCF(table,numSNPUp,discoNameUp,positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,0)
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
        reverseUp=dicoUp[listPolymorphismePosUp[0]][4]
        for comptPol in range(len(listPolymorphismePos)):
            if (int(reverseUp)==-1):
                nucleoLow=ReverseComplement(listnucleoLow[comptPol])
            elif int(reverseUp)==1:
                nucleoLow=listnucleoLow[comptPol]
            nucleoRefLow="."
            positionSnpUp=dicoUp[listPolymorphismePosUp[comptPol]][5]
            nucleoUp=dicoUp[listPolymorphismePosUp[comptPol]][3]
            nucleoRefUp=dicoUp[listPolymorphismePosUp[comptPol]][1]
            table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,0)
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
        reverseLow=dicoLow[listPolymorphismePosLow[0]][4]
        for comptPol in range(len(listPolymorphismePos)):
            if (int(reverseLow)==-1):
                nucleoUp=ReverseComplement(listnucleoUp[comptPol])
            elif int(reverseLow)==1:
                nucleoUp=listnucleoUp[comptPol]
            nucleoRefUp="."
            positionSnpLow=dicoLow[listPolymorphismePosLow[comptPol]][5]
            nucleoLow=dicoLow[listPolymorphismePosLow[comptPol]][3]
            nucleoRefLow=dicoLow[listPolymorphismePosLow[comptPol]][1]
            table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterField,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,0)
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
def GetGenotype(geno,boolRefLow,table,nbGeno,phased,listCovGeno):
    j=0
    genotypes=""
    key=None
    current_genotype=None
    likelihood=None
    if int(nbGeno)==0:
        return table
    else:
            for i in range(0,nbGeno):
                key="G"+str(i+1) # Create the dictionnary key
                current_genotype = geno[key]
                likelihood=current_genotype[1]
                if boolRefLow==True: # check if the mapped path is the lower (in this case exchange 0/0 to 1/1 and 1/1 to 0/0
                    if "1/1" in current_genotype[0]:
                        current_genotype[0]=current_genotype[0].replace("1/1","0/0")
                    elif "0/0" in current_genotype[0]:
                        current_genotype[0]=current_genotype[0].replace("0/0","1/1")
                
                if phased==True: # in case of phasing we change the / symbol
                    current_genotype[0]=current_genotype[0].replace("/","|")
                #TODO REMOVE THIS TEST WHEN WE HAVE ALL FILES WITH RIGHT HEADER
                if isinstance(likelihood,list):
                    likelihood=str(','.join(current_genotype[1]))
                else:
                    likelihood=str(likelihood)
                genotypes+=str(current_genotype[0])+":"+str(listCovGeno[i])+":"+likelihood
                #genotypes+=str(current_genotype[0])+":"+str(listCovGeno[i])+":"+str(','.join(current_genotype[1])) # Add the current genotype
                
                if i<nbGeno-1 :
                    genotypes+="\t" # Add a \t except if this is the last genotype

            table[8]="GT:DP:PL"
            table[9]=genotypes
    return table
##############################################################
##############################################################
def PrintVCFGhost(table,numSNPUp,chrom,position,tp,valRankUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,ntUp,ntLow,geno,nbGeno,phased,listCovGeno,VCF):
    '''Without samfile some fields will have default value : CHROM table[0],POS table[1], QUAL table[5], FILTER [6]''' 
    table[0]=chrom
    table[1]=position
    table[5]="."
    table[6]="."
    
    table[2]=numSNPUp
    table[3]=ntUp
    table[4]=ntLow
    table[7]="Ty="+str(tp)+";"+"Rk="+str(valRankUp)+";"+"UL="+str(unitigLeftUp)+";"+"UR="+str(unitigRightUp)+";"+"CL="+str(contigLeftUp)+";"+"CR="+str(contigRightUp)+";"+str(covUp)
    table=GetGenotype(geno,0,table,nbGeno,phased,listCovGeno)
    table[7]=table[7].replace("None",".")
    table[7]=table[7].replace("none",".")
    printOneline(table,VCF)
    
##############################################################
##############################################################
def FillVCF(table,numSNP,chrom,pos,ref,alt,qual,filterfield,tp,valRank,multi,ok,unitigLeft,unitigRight,contigLeft,contigRight,cov,nucleoRef,reverse,geno,nbGeno,phased,listCovGeno,boolRefLow):
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
    table[7]="Ty="+str(tp)+";"+"Rk="+str(valRank)+";"+"MULTI="+str(multi)+";"+"DT="+str(ok)+";"+"UL="+str(unitigLeft)+";"+"UR="+str(unitigRight)+";"+"CL="+str(contigLeft)+";"+"CR="+str(contigRight)+";"+str(cov)+";"+"Genome="+str(nucleoRef)+";"+"Sd="+str(reverse)
    table[7]=table[7].replace("None",".")
    table[7]=table[7].replace("none",".")
    table=GetGenotype(geno,boolRefLow,table,nbGeno,phased,listCovGeno)
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
        data=str(bin(int(FLAG)))[::-1]
        try:
                return(data[bit]=='1')
        except IndexError:
                return False
##############################################################
##############################################################
def CheckContigUnitig(unitig,contig):
        if contig:
                return(int(contig))
        elif unitig:
                return(int(unitig))
        else:
                return 0







