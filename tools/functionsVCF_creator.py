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
            for i in nomDisco:
                if "nb_pol" in i :
                    nb_pol=i.split('_')
                    j=0
                    for j in nb_pol:
                        matchInt=re.match(r'^\d+$',j)
                        if matchInt:
                            nbSnp=nbSnp+int(j)-1
            if nbGeno==0:
                for k in nomDisco:
                    match=re.match(r'^G',k)
                    if match:
                        nbGeno+=1
    return(int(nbSnp)/2,nbGeno)

##############################################################
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
    for i in snp:
        if 'SNP' or 'INDEL' in i:
            nomDisco=i.split('|')
            break
    i=0
    for i in nomDisco:
        if "SNP" in i:
            SNP=1
            ID= i.split('_')
            for j in ID:
                matchInt=re.match(r'^\d+$',j)
                if matchInt:
                    numSNP=j
        elif 'P_' in i:
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
                    dicoHeader[key]=[posD,ntUp,ntLow]
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
                    dicoHeader[key]=[posD,ind,amb]
        elif "INDEL" in i:
            INDEL=True
            ID= i.split('_')
            j=0
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
        elif "C" in i:
            coverage=i.split('_')
            j=0
            for j in range(len(coverage)):
                matchInt=re.match(r'^\d+$',coverage[j])
                if matchInt:
                    listCoverture.append(str(coverage[j]))
                    listC.append(str(coverage[j-1]))
        elif "nb_pol" in i :
            nb_pol=i.split('_')
            j=0
            for j in nb_pol:
                matchInt=re.match(r'^\d+$',j)
                if matchInt:
                    nb_pol=j
        elif "G" in i:
            i=i.replace("_",":")
            listgeno=i.split(":")
            if len(listgeno)>2:
                listlikelihood=listgeno[2].split(",")
            else:
                listlikelihood=0
            dicoGeno[listgeno[0]]=[listgeno[1],listlikelihood]
    
    if boolNum==0:
        return(snp,numSNP,unitigLeft,unitigRight,contigLeft,contigRight,valRank,listCoverture,listC,nb_pol,ln,posD,ntUp,ntLow,dicoGeno,dicoHeader)
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
    return(posUp,posLow,snpLow,snpUp,boolXAUp,boolXALow)


##############################################################
#listCup=["C1","C2"]
#listCovertureUp=[23,45]
##############################################################
def GetCoverage( listCUp, listCLow, listCoverageUp, listCoverageLow):
    """Takes as input lists list  """
    i=0
    covUp=''
    listCovGeno=[]
    ##Create a string if the number of coverage is superior to 1
    if len( listCoverageUp)>1 and len( listCoverageLow)>1:
        while i<len( listCoverageUp):
            listCovGeno.append(int( listCoverageUp[i])+int( listCoverageLow[i]))
            if covUp!='':
                covUp=str(covUp)+';'+str( listCUp[i])+':'+str( listCoverageUp[i])+"|"+str( listCoverageLow[i])
            else:
                covUp=str( listCUp[i])+':'+str( listCoverageUp[i])+"|"+str( listCoverageLow[i])
            i+=1
        i=0
        covLow=''
        while i<len( listCoverageLow):
            if covLow!='':
                covLow=str(covLow)+';'+str( listCLow[i])+':'+str( listCoverageLow[i])+"|"+str( listCoverageUp[i])
            else:
                covLow=str( listCLow[i])+':'+str( listCoverageLow[i])+"|"+str( listCoverageUp[i])
            i+=1
    else:
        covUp=str( listCUp[0])+':'+str( listCoverageUp[0])
        covLow=str( listCLow[0])+':'+str( listCoverageLow[0])
    return(covUp,covLow,listCovGeno) #string covUp C1:5|23;C2:35|1 listCovGeno=[28,36]
##############################################################
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

##############################################################
##############################################################
def ValidationSNP(snpLow,posLow,snpUp,posUp,NM):
    """Main function of the snp validation : check if the couple is validated with only one mapping position at the number of mismatch tested"""
    couple=None
    listPos= []
    ensemble=None
    delta=3
    if abs(int(snpUp[3]))<=0 and abs(int(snpLow[3]))<=0:
        couple="unmapped"
        return(couple)
    else:
        for position,nbMismatch in posUp.items(): # Checks for the upper path
            if nbMismatch==int(NM):
                listPos=addPosition(position,listPos,delta) #If the position is mapped at NM (number of mismatch tested) test if its not too close with an other position
                ensemble=set(listPos)
                if abs(len(ensemble))>1: #case of many position of mapping for NM => multiple mapped
                    couple="multiple"
                    return(couple)
        for position,nbMismatch in posLow.items():# Checks for the lower path
            if nbMismatch==int(NM):
                listPos=addPosition(position,listPos,delta)
                ensemble=set(listPos)
                if abs(len(ensemble))>1:
                    couple="multiple"
                    return(couple)
        if ensemble!=None:
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
##############################################################
def CigarCodeChecker(cigarcode,seq,listpol,posModif,indel):
    """Function which allows to recover the position of the SNPs according to the position of the deletions or insertions relative to the reference sequence (CigarCode Parsing : checks for insertion deletion or soft clipping"""
    parsingCigarCode=re.findall('(\d+|[A-Za-z])',cigarcode)
    listPosRef=[]
    listShift=[]
    somme=0
    shift=0
    pos=0
    i=1
    j=0
    #Close snps
    if len(listpol)>1:
        while i<len(parsingCigarCode):
            #Soft clipping
            if parsingCigarCode[i]=="S":
                shift-=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            #Match or Mismatch
            elif parsingCigarCode[i]=="M":
                pos+=int(parsingCigarCode[i-1])
            #Deletion
            elif parsingCigarCode[i]=="D":
                shift+=int(parsingCigarCode[i-1])
            #Insertion
            elif parsingCigarCode[i]=="I":
                shift-=int(parsingCigarCode[i-1])
            while int(pos)>=int(listpol[j]):
                posRef=int(listpol[j])+shift
                listPosRef.append(posRef)
                listShift.append(shift)
                if j<(len(listpol)-1):
                    j+=1
                else:
                    return(listPosRef,listShift)
            i+=2
    #Simple snp and Indel
    else:
        lenDemiSeq=posModif
        while i<len(parsingCigarCode):
            if parsingCigarCode[i]=="S":
                shift-=int(parsingCigarCode[i-1])
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="M":
                pos+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="D":
                shift+=int(parsingCigarCode[i-1])
            elif parsingCigarCode[i]=="I":
                shift-=int(parsingCigarCode[i-1])
            if len(listpol)==1:
                if pos>=int(lenDemiSeq):
                    posCentraleRef=int(lenDemiSeq)+shift
                    return(posCentraleRef,shift)
            i+=2
    return(listPosRef,listShift)
##############################################################
##############################################################
def ReferenceChecker(shift,posMut,posCentraleRef):
    """MD tag parsing; checks if path nucleotide is identical to the reference nucleotide """
    pos=0
    i=0
    matchInt=None
    nucleoRef=None
    boolEgalRef=None
    if '^' in posMut:
        pos=shift
        motifDel=re.compile("[0-9]*\^[A-Za-z]*")
        dictDel={}
        deletion=re.findall('\d+\^[A-Za-z]*',posMut)
        deletion=''.join(deletion)
        numeroDel=re.findall('\d+',deletion)
        parsingPosMut=re.findall('(\d+|[A-Za-z]|\^)',posMut)
        strNumerodel=''.join(numeroDel)
        for m in motifDel.finditer(posMut):
            dictDel[m.group()]=parsingPosMut.index('^')-1
        posMut=posMut.replace(deletion,"")
        parsingPosMut=re.findall('(\d+|[A-Za-z])',posMut)
        parsingPosMut.insert(dictDel[deletion],numeroDel[0])
    else:
        parsingPosMut=re.findall('(\d+|[A-Za-z])',posMut)
    while i<len(parsingPosMut):
        matchInt=re.match(r'\d+',parsingPosMut[i])
        if matchInt:
            parsingPosMut[i]=int(parsingPosMut[i])
        i+=1
    i=0
    while i<len(parsingPosMut):
        if isinstance(parsingPosMut[i],int):
            pos+=parsingPosMut[i]
        else :
            pos+=1
        if pos==posCentraleRef:
            if isinstance(parsingPosMut[i],str):
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
##############################################################
##############################################################
def GetSequence(snpUp,snpLow):
    """Get reverse sequence"""
    seqLow=snpLow[9]
    seqUp=snpUp[9]
    if snpUp[1]=="16":
        i=0
        listSeqUp=list(seqUp)
        seqUp=''
        while i<len(listSeqUp):
            if seqUp!='':
                seqUp=str(ReverseComplement(listSeqUp[i]))+seqUp
            else :
                seqUp=str(ReverseComplement(listSeqUp[i]))
            i+=1
    if snpLow[1]=="16":
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

##############################################################
#dicoHeader snps : dicoHeader[key]=[posD,ntUp,ntLow] 
#dicoHeader indel : dicoHeader[key]=[posD,ind,amb]
##############################################################
def RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel):
    #INIT VALUES
    posSNPUp=None
    posSNPLow=None
    boolRefUp=None
    boolRefLow=None
    nucleoUp=None
    nucleoLow=None
    nucleoRefUp=None
    nucleoRefLow=None
    listPolymorphismePosUp=None
    listPolymorphismePosLow=None
    
    posModif=dicoHeaderUp["P_1"][0] #Position of the first variant gived by discosnp
    seqUp=list(snpUp[9])
    seqLow=list(snpLow[9])
    reverseUp=1
    reverseLow=1
    dicopolUp={}
    dicopolLow={}
    nucleoLow=None
    nucleoUp=None
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
#Shift by positions (insertion,deletion,sofclipping) and update of the position in alignment
    if int(snpUp[3])>0:
        #Check cigarCode : Presence of insertion, softclipping, deletion
        if snpUp[1]=="0": #Forward Strand
            posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],seqUp,listPos,posModif,indel) # Gets the positions of the variants with an eventual shift from the reference (insertion,deletion,soft clipping) ; in case of close snps return a list 
            listPolymorphismePosUp=listPos
            if len(listPos)==1 and indel==False: #simple snp
                nucleoUp=listnucleoUp[0] # Gets the nucleotide 
        elif snpUp[1]=="16":# Reverse Strand
            reverseUp=-1
            posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],seqUp,listPosR,posModif,indel)
            listPolymorphismePosUp=listPosR #List of all the reverse position
            if len(listPos)==1 and indel==False:
                nucleoUp=listnucleoUpR[0]
    else:
        listPolymorphismePosUp=listPos
        posCentraleUp=listPos
        shiftUp=None
        reverseUp="."
    if int(snpLow[3])>0:
        #Check cigarCode : Presence of insertion, softclipping, deletion
        if snpLow[1]=="0":
            posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],seqLow,listPos,posModif,indel)
            listPolymorphismePosLow=listPos
            if len(listPos)==1 and indel==False:
                nucleoLow=listnucleoLow[0]
        elif snpLow[1]=="16":
            reverseLow=-1
            posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],seqLow,listPosR,posModif,indel)
            listPolymorphismePosLow=listPosR
            if len(listPos)==1 and indel==False:
                nucleoLow=listnucleoLowR[0]
    else:
        listPolymorphismePosLow=listPos
        posCentraleLow=listPos
        shiftLow=None
        reverseLow="."
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
    ####Up Case : Check for all the variants if there are identical to the reference and proc 
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
        positionSnpUp = posCentraleUp
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
####Low 
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
        positionSnpLow = posCentraleLow

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
    if len(listPos)==1:
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
        ##Check the strand : if the alternative nucleotide is not on the same strand as the reference : reverse it
        if boolRefLow==True and reverseLow==-1 and reverseUp==1:
               nucleoUp=ReverseComplement(nucleoUp)
        elif boolRefLow==True and reverseLow==1 and reverseUp==-1:
               nucleoUp=ReverseComplement(nucleoUp)
        if boolRefUp==True and reverseUp==-1 and reverseLow==1:
               nucleoLow=ReverseComplement(nucleoLow)
        elif boolRefUp==True and reverseUp==1 and reverseLow==-1:
               nucleoLow=ReverseComplement(nucleoLow)
        return(nucleoLow,posSNPLow,nucleoUp,posSNPUp,boolRefLow,boolRefUp,reverseUp,reverseLow,nucleoRefUp,nucleoRefLow)
    else:
        return(dicopolUp,dicopolLow,listPolymorphismePosUp,listPolymorphismePosLow)
##############################################################
##############################################################
def MismatchChecker(snpUp,posUp,snpLow,posLow,nucleoRefUp,nucleoRefLow,nucleoUp,nucleoLow,boolRefUp,boolRefLow,indel):
    """In case of divergent main position (case snp whose two paths are mapped ) to define the reference = > check the number of mismatch
( If the number of mismatch is the same in both cases it is the lower lexicographical SNP which is selected for reference .
The Boolean allows to know the reference SNP ) """
    posRef=0
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
    snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp, listCoverageUp, listCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
    snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow, listCoverageLow,listClow,nb_polLow,lnlow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
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
            else:
              table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterfield,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)  
#---------------------------------------------------------------------------------------------------------------------------
##Case : Lower path mapped and upper path unmapped      
    elif int(snpUp[3])<=0 and int(snpLow[3])>0:
        table=FillVCF(table,numSNPUp,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterfield,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,boolRefLow)
        ##If the number of mismatch is the max : that means the second path will be mapped at number of mismatch +1
        if dmax:
                table[4]=nucleoUp
        else:
                table[4]='.'
        if snpLow[10]=="*":
                table[5]="."
        else:
                table[5]=snpLow[10]  
#---------------------------------------------------------------------------------------------------------------------------
##Case : Upper path mapped and lower path unmapped        
    elif int(snpUp[3])>0 and int(snpLow[3])<=0:
        table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterfield,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,boolRefLow)
        if dmax:
                table[4]=nucleoLow
        else:
                table[4]='.'
        if snpUp[10]=="*":
                table[5]="."
        else:
                table[5]=snpUp[10]
#---------------------------------------------------------------------------------------------------------------------------
##Case : Both paths are unmapped       
    elif int(snpUp[3])<=0 and int(snpLow[3])<=0:
        if nucleoLow<nucleoUp:
            table=FillVCF(table,numSNPLow,".",".",nucleoLow,nucleoUp,snpLow[10],filterfield,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,".",".",geno,nbGeno,phased,listCovGeno,boolRefLow)
        else:
            table=FillVCF(table,numSNPUp,".",".",nucleoUp,nucleoLow,snpUp[10],filterfield,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,".",".",geno,nbGeno,phased,listCovGeno,boolRefLow)
    return(table)
##############################################################
##############################################################   
def printOneline(table,VCF):
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
    ##Sorts the list of position to get the smallest one and its position in the list unsorted!
    listSortedPos=list(listPolymorphismePos)
    listSortedPos.sort()
    indexSmallestPos=listPolymorphismePos.index(listSortedPos[0])
    ##Keep the close snps to sort them : indeed all the lists and dictionnaries :listnucleoUp,listPolymorphismePosUp,listPolymorphismePosLow,listnucleoLow dicoUp,dicoLow are classified according to dicoHeader so if we start by sorting we lose the correspondence between data
    tablebis = []
    k=0
    #Variables
    snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp, listCoverageUp, listCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
    snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow, listCoverageLow, listCLow,nb_polLow,lnLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
#---------------------------------------------------------------------------------------------------------------------------
##Case : two mapped paths
    if int(snpUp[3])>0 and int(snpLow[3])>0:
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
        for l in range(len(tablebis)):
                printOneline(tablebis[l],VCF)
#---------------------------------------------------------------------------------------------------------------------------
##Case : Both paths are unmapped     
    elif int(snpUp[3])<=0 and int(snpLow[3])<=0:
        reverseUp="."
        i=0
        for i in range(len(listPolymorphismePos)):
            positionSnpUp="+"+str(listPolymorphismePos[i])
            positionSnpLow="+"+str(listPolymorphismePos[i])
            nucleoUp=listnucleoUp[i]
            nucleoRefUp="."
            nucleoLow=listnucleoLow[i]
            nucleoRefLow="."
            table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,0)
            tablebis.append(list(table))
        tablebis=sorted(tablebis, key=lambda colonnes: colonnes[1])
        l=0
        for l in range(len(tablebis)):
                printOneline(tablebis[l],VCF)
#---------------------------------------------------------------------------------------------------------------------------
##Case : Upper path mapped and lower path unmapped  
    elif int(snpUp[3])>0 and int(snpLow[3])<=0:
        comptPol=0
        reverseUp=dicoUp[listPolymorphismePosUp[0]][4]
        for comptPol in range(len(listPolymorphismePos)):
            nucleoLow="."
            nucleoRefLow="."
            positionSnpUp=dicoUp[listPolymorphismePosUp[comptPol]][5]
            nucleoUp=dicoUp[listPolymorphismePosUp[comptPol]][3]
            nucleoRefUp=dicoUp[listPolymorphismePosUp[comptPol]][1]
            table=FillVCF(table,numSNPUp,snpUp[2],positionSnpUp,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,multi,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,geno,nbGeno,phased,listCovGeno,0)
            if (int(reverseUp)==-1) and dmax:
                table[4]=ReverseComplement(seqLow[int(listPolymorphismePosLow[comptPol])-1])
            elif int(reverseUp)==1 and dmax:
                table[4]=seqLow[int(listPolymorphismePosLow[comptPol])-1]
            else:
                table[4]='.'
            tablebis.append(list(table))
        tablebis=sorted(tablebis, key=lambda colonnes: colonnes[1])
        l=0
        for l in range(len(tablebis)):
                printOneline(tablebis[l],VCF)
#---------------------------------------------------------------------------------------------------------------------------
##Case : Lower path mapped and upper path unmapped            
    elif int(snpUp[3])<=0 and int(snpLow[3])>0:
        reverseLow=dicoLow[listPolymorphismePosLow[0]][4]
        for comptPol in range(len(listPolymorphismePos)):
            nucleoUp="."
            nucleoRefUp="."
            positionSnpLow=dicoLow[listPolymorphismePosLow[comptPol]][5]
            nucleoLow=dicoLow[listPolymorphismePosLow[comptPol]][3]
            nucleoRefLow=dicoLow[listPolymorphismePosLow[comptPol]][1]
            table=FillVCF(table,numSNPLow,snpLow[2],positionSnpLow,nucleoLow,nucleoUp,snpLow[10],filterField,tp,valRankLow,multi,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,geno,nbGeno,phased,listCovGeno,0)
            if (int(reverseLow)==-1) and dmax:
                table[4]=ReverseComplement(seqUp[int(listPolymorphismePosUp[comptPol])-1])
            elif int(reverseLow)==1 and dmax:
                table[4]=seqUp[int(listPolymorphismePosUp[comptPol])-1]
            else:
                table[4]='.'
            tablebis.append(list(table))
        tablebis=sorted(tablebis, key=lambda colonnes: colonnes[1])
        l=0
        for l in range(len(tablebis)):
                printOneline(tablebis[l],VCF)   
            
##############################################################
# geno: [key][0] : genotype (0/0, 0/1 or 1/1)
# geno: [key][1] : likelihood for the 3 possible genotypes ("x,y,z")
##############################################################
def GetGenotype(geno,boolRefLow,table,nbGeno,phased,listCovGeno):
    j=0
    genotypes=""
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
def PrintVCFGhost(table,numSNPUp,tp,valRankUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,ntUp,ntLow,geno,nbGeno,phased,listCovGeno,VCF):
    '''Without samfile some fields will have default value : CHROM table[0],POS table[1], QUAL table[5], FILTER [6]''' 
    table[0]="."
    table[1]="."  
    table[5]="."
    table[6]="."
    
    table[2]=numSNPUp
    table[3]=ntUp
    table[4]=ntLow
    table[7]="Ty:"+str(tp)+";"+"Rk:"+str(valRankUp)+";"+"UL:"+str(unitigLeftUp)+";"+"UR:"+str(unitigRightUp)+";"+"CL:"+str(contigLeftUp)+";"+"CR:"+str(contigRightUp)+";"+str(covUp)
    table=GetGenotype(geno,0,table,nbGeno,phased,listCovGeno)
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
    table[7]="Ty:"+str(tp)+";"+"Rk:"+str(valRank)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeft)+";"+"UR:"+str(unitigRight)+";"+"CL:"+str(contigLeft)+";"+"CR:"+str(contigRight)+";"+str(cov)+";"+"Genome:"+str(nucleoRef)+";"+"Sd:"+str(reverse)
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













