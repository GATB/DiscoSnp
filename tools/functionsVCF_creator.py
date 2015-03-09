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
		snp=snp.rstrip('\n')
		snp=snp.split('\t')
	listeCouverture=[]
	listeC=[]
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
	INDEL=0
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
			INDEL=1
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
			couverture=i.split('_')
			j=0
			for j in range(len(couverture)):
				matchInt=re.match(r'^\d+$',couverture[j])
				if matchInt:
					listeCouverture.append(str(couverture[j]))
					listeC.append(str(couverture[j-1]))
		elif "nb_pol" in i :
			nb_pol=i.split('_')
			j=0
			for j in nb_pol:
				matchInt=re.match(r'^\d+$',j)
				if matchInt:
					nb_pol=j
		elif "G" in i:
			i=i.replace("_",":")
			listegeno=i.split(":")
			if len(listegeno)>2:
				listelikelihood=listegeno[2].split(",")
			else:
				listelikelihood=0
			dicoGeno[listegeno[0]]=[listegeno[1],listelikelihood]
			
	if boolNum==0:				
		return(snp,numSNP,unitigLeft,unitigRight,contigLeft,contigRight,valRank,listeCouverture,listeC,nb_pol,ln,posD,ntUp,ntLow,dicoGeno,dicoHeader)
	else:
		return(numSNP)
##############################################################
##############################################################
def isCouple(numSNP1,numSNP2):
	""" Verifies that the IDs of SNPs Up and Low correspond to the same path"""
	getcouple=1
	if numSNP1==numSNP2:
		getcouple=1
	else : getcouple=0
	return(getcouple)
##############################################################
##############################################################
def GetCouple(snp1,snp2):
	"""Retrieves for each path alignment information in a list ; retrieves a dictionary with all the positions of a path and the number of associated mismatch ; retrieves a Boolean indicating whether the SNP is multimapped or not"""
	posUp = {}
	posLow = {}
	i=0
	j=0
	#Boolean snps mapped
	boolMapUp=1
	boolMapLow=1
	#Boolean snps mulimapped
	boolXAUp=0
	boolXALow=0
	if "lower" in snp1[0]:
		snpLow=snp1
	else:
		snpUp=snp1
	if "lower" in snp2[0]:
		snpLow=snp2
	else:
		snpUp=snp2
	#Error list with mapping position very close 
	listeerreurUp=set([(int(snpUp[3])-1),(int(snpUp[3])+1),(int(snpUp[3])+2),(int(snpUp[3])+3),(int(snpUp[3])-3),(int(snpUp[3])-2),int(snpUp[3])])
	listeerreurLow=set([(int(snpLow[3])-1),(int(snpLow[3])+1),(int(snpLow[3])+2),(int(snpLow[3])+3),(int(snpLow[3])-3),(int(snpLow[3])-2),int(snpLow[3])])
	#Creation of a dict with mapping position associate with number of mismatch
	if 'XA' in ''.join(snpUp):
		i=0
		#parsing XA tag
		listeXA=snpUp[19].split(';')
		strXA = ','.join(listeXA)
		listeXA = strXA.split(',')
		listeXA.pop()
		position=listeXA[1:]
		while i<len(position):
			if abs(int(position[i])) not in listeerreurUp :
				boolXAUp=1
				posUp[abs(int(position[i]))]=int(position[i+2])
			i+=4
	if 'XA' in ''.join(snpLow):
		i=0
		listeXA=snpLow[19].split(';')
		strXA = ','.join(listeXA)
		listeXA = strXA.split(',')
		listeXA.pop()
		position=listeXA[1:]
		while i<len(position):
			if abs(int(position[i])) not in listeerreurLow :
				boolXALow=1
				posLow[abs(int(position[i]))]=int(position[i+2])
			i+=4
	#Add main position to the dict if the path is mapped
	if abs(int(snpUp[3]))>0:
		poubelle,poubelle,nbMismatchUp=snpUp[12].split(":")
		posUp[abs(int(snpUp[3]))]=int(nbMismatchUp)
		
	else : boolMapUp=0
	if abs(int(snpLow[3]))>0:
		poubelle,poubelle,nbMismatchLow=snpLow[12].split(":")
		posLow[abs(int(snpLow[3]))]=int(nbMismatchLow)	
	else : boolMapLow=0
	return(posUp,posLow,snpLow,snpUp,boolMapUp,boolMapLow,boolXAUp,boolXALow)


##############################################################
##############################################################
def GetCoverage(listeCUp,listeCLow,listeCouvertureUp,listeCouvertureLow):
	i=0
	couvUp=''
	listeCouvGeno=[]
	if len(listeCouvertureUp)>1 and len(listeCouvertureLow)>1:
		while i<len(listeCouvertureUp):
			listeCouvGeno.append(int(listeCouvertureUp[i])+int(listeCouvertureUp[i]))
			if couvUp!='':
				couvUp=str(couvUp)+';'+str(listeCUp[i])+':'+str(listeCouvertureUp[i])+"|"+str(listeCouvertureLow[i])
			else:
				couvUp=str(listeCUp[i])+':'+str(listeCouvertureUp[i])+"|"+str(listeCouvertureLow[i])
			i+=1
		i=0
		couvLow=''
		while i<len(listeCouvertureLow):
			if couvLow!='':
				couvLow=str(couvLow)+';'+str(listeCLow[i])+':'+str(listeCouvertureLow[i])+"|"+str(listeCouvertureUp[i])
			else:
				couvLow=str(listeCLow[i])+':'+str(listeCouvertureLow[i])+"|"+str(listeCouvertureUp[i])
			i+=1
	else:
		couvUp=str(listeCUp[0])+':'+str(listeCouvertureUp[0])
		couvLow=str(listeCLow[0])+':'+str(listeCouvertureLow[0])
	return(couvUp,couvLow,listeCouvGeno)
##############################################################
##############################################################
def addPosition(position,ensemble,delta):
	if len(ensemble)==0:
		ensemble.append(position)
		return(ensemble)
	if abs(position-ensemble[0])>delta:
		ensemble.append(position)
		return(ensemble)
	return(ensemble)
		
##############################################################
##############################################################
def ValidationSNP(snpLow,posLow,snpUp,posUp,NM):
	couple=None
	listePos= []
	ensemble=None
	delta=3
	if abs(int(snpUp[3]))<=0 and abs(int(snpLow[3]))<=0:
		couple="unmapped"
		return(couple)
	else:
		for position,nbMismatch in posUp.items():
			if nbMismatch==int(NM):
				listePos=addPosition(position,listePos,delta)
				ensemble=set(listePos)
				if abs(len(ensemble))>1:
					couple="multiple"
					return(couple)	
		for position,nbMismatch in posLow.items():
			if nbMismatch==int(NM):
				listePos=addPosition(position,listePos,delta)
				ensemble=set(listePos)
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
def CigarCodeChecker(cigarcode,seq,listepol,posModif,indel):
	"""Function which allows to recover the position of the SNPs according to the position of the deletions or insertions relative to the reference sequence (CigarCode Parsing : checks for insertion deletion or soft clipping"""
	parsingCigarCode=re.findall('(\d+|[A-Za-z])',cigarcode)
	listePosRef=[]
	listeShift=[]
	somme=0
	shift=0
	pos=0
	i=1
	j=0
	#Close snps
	if len(listepol)>1:
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
			while int(pos)>=int(listepol[j]):
				posRef=int(listepol[j])+shift
				listePosRef.append(posRef)
				listeShift.append(shift)
				if j<(len(listepol)-1):
					j+=1
				else:
					return(listePosRef,listeShift)
			i+=2
	#Lonely snp
	else:
		lenDemiSeq=int(listepol[0])
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
			if len(listepol)==1:
				if pos>=int(lenDemiSeq):
					posCentraleRef=int(lenDemiSeq)+shift
					return(posCentraleRef,shift)
			i+=2
	return(listePosRef,listeShift)
##############################################################
##############################################################
def ReferenceChecker(shift,posMut,posCentraleRef):
	"""MD tag parsing; checks if path nucleotide is identical to the reference nucleotide """
	pos=0
	i=0
	matchInt=None
	lettreRef=None
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
				boolEgalRef=0
				lettreRef=parsingPosMut[i]
			else:
				boolEgalRef=1
			break
		if pos>posCentraleRef:
			boolEgalRef=1
			break
		i+=1
	return(boolEgalRef,lettreRef)
##############################################################
##############################################################
def GetSequence(snpUp,snpLow):
	"""Get reverse sequence"""
	seqLow=snpLow[9]
	seqUp=snpUp[9]
	if snpUp[1]=="16":
		i=0
		listeSeqUp=list(seqUp)
		seqUp=''	
		while i<len(listeSeqUp):
			if seqUp!='':
				seqUp=str(ReverseComplement(listeSeqUp[i]))+seqUp
			else :
				seqUp=str(ReverseComplement(listeSeqUp[i]))
			i+=1
	if snpLow[1]=="16":
		i=0
		listeSeqLow=list(str(seqLow))
		seqLow=''	
		while i<len(listeSeqLow):
			if seqLow!='':
				seqLow=str(ReverseComplement(listeSeqLow[i]))+seqLow
			else :
				seqLow=str(ReverseComplement(listeSeqLow[i]))
			i+=1
	return(seqUp,seqLow)

##############################################################
##############################################################

def RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel):
	seqUp=list(snpUp[9])
	seqLow=list(snpLow[9])
	reverseUp=1
	reverseLow=1
	dicopolUp={}
	dicopolLow={}
	lettreLow=None
	lettreUp=None
	bug=0
	erreur=0
	i=0
	if int(snpUp[3])>0:
		posRef=int(snpUp[3])
	else:
		posRef=int(snpLow[3])	
#---------------------------------------------------------------------------------------------------------------------------
	#List of snps
	if indel==0:
		listePos,listeLettreUp,listeLettreLow,listePosR,listeLettreUpR,listeLettreLowR=GetPolymorphisme(dicoHeaderUp,snpUp[9],indel)
	else:
		seqUp=snpUp[9]
		seqLow=snpLow[9]
		if len(seqUp)<len(seqLow):
			seq=seqLow
		else:
			seq=seqUp
		listePos,listePosR,insert,ntStart,ambiguityPos=GetPolymorphisme(dicoHeaderUp,seq,indel)
	if indel==1:
		posModif=int(listePos[0])
	else:
		posModif=0
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#Shift by positions (insertion,deletion,sofclipping) and update of the position in alignment 
 	if int(snpUp[3])>0:
		#Check cigarCode : Presence of insertion, softclipping, deletion
		if snpUp[1]=="0":
			posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],seqUp,listePos,posModif,indel)
			listePolymorphismePosUp=listePos
			if len(listePos)==1 and indel==0:
				lettreUp=listeLettreUp[0]
		elif snpUp[1]=="16":
			reverseUp=-1
			posCentraleUp,shiftUp=CigarCodeChecker(snpUp[5],seqUp,listePosR,posModif,indel)
			listePolymorphismePosUp=listePosR
			if len(listePos)==1 and indel==0:
				lettreUp=listeLettreUpR[0]
	else:
		posCentraleUp=listePos
		shiftUp=None
		reverseUp="."
	if int(snpLow[3])>0:
		#Check cigarCode : Presence of insertion, softclipping, deletion
		if snpLow[1]=="0":
			posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],seqLow,listePos,posModif,indel)	
			listePolymorphismePosLow=listePos
			if len(listePos)==1 and indel==0:
				lettreLow=listeLettreLow[0]
		elif snpLow[1]=="16":
			reverseLow=-1
			posCentraleLow,shiftLow=CigarCodeChecker(snpLow[5],seqLow,listePosR,posModif,indel)
			listePolymorphismePosLow=listePosR
			if len(listePos)==1 and indel==0:
				lettreLow=listeLettreLowR[0]
	else:
		posCentraleLow=listePos
		shiftLow=None
		reverseLow="."
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
	####UP : Test position et SNP : identique ou non à la séquence de référence
	if int(snpUp[3])>0:
		#Extraction tag MD
		MD,Z,posMutUp = snpUp[18].split(":")
		if len(listePos)>1:
			#Dictionnaire contenant la position des polymorphismes(clés) et sous forme de liste boolRefUp,lettreRefUp,posCentraleUp[i],listeLettreUp[i],reverseUp
			if reverseUp==1:
				i=0
				for i in range(len(listePos)):
					boolRefUp,lettreRefUp=ReferenceChecker(shiftUp[i],posMutUp,posCentraleUp[i])
					if lettreRefUp==None:
						lettreRefUp=listeLettreUp[i]
					dicopolUp[listePos[i]]=[boolRefUp,lettreRefUp,posCentraleUp[i],listeLettreUp[i],reverseUp,(int(snpUp[3])+posCentraleUp[i])]
			elif reverseUp==-1:
				i=0
				for i in range(len(listePosR)):
					boolRefUp,lettreRefUp=ReferenceChecker(shiftUp[i],posMutUp,posCentraleUp[i])
					if lettreRefUp==None:
						lettreRefUp=listeLettreUpR[i]
					dicopolUp[listePosR[i]]=[boolRefUp,lettreRefUp,posCentraleUp[i],listeLettreUpR[i],reverseUp,(int(snpUp[3])+posCentraleUp[i])]
		else:
			boolRefUp,lettreRefUp=ReferenceChecker(shiftUp,posMutUp,posCentraleUp)
			if lettreRefUp==None:
				lettreRefUp=seqUp[posModif]

	elif int(snpUp[3])<=0 and len(listePos)==0 : 
		lettreUp = seqUp[int(listePos[0])-1]
		positionSnpUp = posCentraleUp
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
####Low : Test position et SNP : identique ou non à la séquence de référence
	if int(snpLow[3])>0:
		#Extraction tag MD
		MD,Z,posMutLow = snpLow[18].split(":")
		if len(listePos)>1:
			#Dictionnaire contenant la position des polymorphismes(clés) et sous forme de liste boolRefLow,lettreRefLow,posCentraleLow[i],listeLettreLow[i],reverseLow
			if reverseLow==1:
				i=0
				for i in range(len(listePos)):
					boolRefLow,lettreRefLow=ReferenceChecker(shiftLow[i],posMutLow,posCentraleLow[i])
					if lettreRefLow==None:
						lettreRefLow=listeLettreLow[i]
					dicopolLow[listePos[i]]=[boolRefLow,lettreRefLow,posCentraleLow[i],listeLettreLow[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
			elif reverseLow==-1:
				i=0
				for i in range(len(listePosR)):
					boolRefLow,lettreRefLow=ReferenceChecker(shiftLow[i],posMutLow,posCentraleLow[i])
					if lettreRefLow==None:
						lettreRefLow=listeLettreLowR[i]
					dicopolLow[listePosR[i]]=[boolRefLow,lettreRefLow,posCentraleLow[i],listeLettreLowR[i],reverseLow,(int(snpLow[3])+posCentraleLow[i])]
		else:
			boolRefLow,lettreRefLow=ReferenceChecker(shiftLow,posMutLow,posCentraleLow)
			if lettreRefLow==None:
				lettreRefLow=seqLow[int(posModif)-1]

	elif int(snpLow[3])<=0 and len(listePos)==0 : 
		lettreLow = seqLow[listePos[0]]
		positionSnpLow = posCentraleLow

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
####Two paths mapped at different main position (first given by bwa)
	if (int(snpUp[3])>0 and int(snpLow[3])>0) and (int(snpUp[3])!=int(snpLow[3])):
		if len(listePos)>1:
			i=0
			for i in range(len(listePos)):
				boolRefLow,boolRefUp,lettreRefUp,lettreRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,dicopolUp[listePolymorphismePosUp[i]][1],dicopolLow[listePolymorphismePosLow[i]][1],dicopolUp[listePolymorphismePosUp[0]][3],dicopolLow[listePolymorphismePosLow[0]][3],dicopolUp[listePolymorphismePosUp[i]][0],dicopolLow[listePolymorphismePosLow[i]][0],indel)
				dicopolUp[listePolymorphismePosUp[i]][0]=boolRefUp
				dicopolUp[listePolymorphismePosUp[i]][1]=lettreRefUp
				dicopolLow[listePolymorphismePosLow[i]][0]=boolRefLow
				dicopolLow[listePolymorphismePosLow[i]][1]=lettreRefLow
				dicopolLow[listePolymorphismePosLow[i]][5]=int(dicopolLow[listePolymorphismePosLow[i]][2])+int(posRef)
				dicopolUp[listePolymorphismePosUp[i]][5]=int(dicopolUp[listePolymorphismePosUp[i]][2])+int(posRef)
		else:
				boolRefLow,boolRefUp,lettreRefUp,lettreRefLow,posRef=MismatchChecker(snpUp,posUp,snpLow,posLow,lettreRefUp,lettreRefLow,lettreUp,lettreLow,boolRefUp,boolRefLow,indel)
#---------------------------------------------------------------------------------------------------------------------------	#---------------------------------------------------------------------------------------------------------------------------
	#SNPs positions
	if len(listePos)==1:
		if indel==1:
			lettreLow="."
			lettreUp="."
			lettreRefUp="."
			lettreRefLow="."
		if (int(snpUp[3])>0 and int(snpLow[3])>0):
			if indel==1:
				posSNPUp=posCentraleUp+posRef-int(ambiguityPos)
				posSNPLow=posCentraleLow+posRef-int(ambiguityPos)
				if snpUp[9]<snpLow[9]:
					boolRefUp=1
					boolRefLow=0
				else:
					boolRefUp=0
					boolRefLow=1
			else:
				posSNPUp=posCentraleUp+posRef
				posSNPLow=posCentraleLow+posRef
		elif int(snpUp[3])<=0 and int(snpLow[3])>0:
			if indel==1:
				boolRefUp=0
				boolRefLow=1
				posSNPUp=None
				posSNPLow=posCentraleLow+posRef-int(ambiguityPos)
			else:
				posSNPUp=None
				posSNPLow=posCentraleLow+posRef
		elif int(snpUp[3])>0 and int(snpLow[3])<=0:
			if indel==1:
				boolRefUp=1
				boolRefLow=0
				posSNPLow=None
				posSNPUp=posCentraleUp+posRef-int(ambiguityPos)
			else:
				posSNPLow=None
				posSNPUp=posCentraleUp+posRef
		return(lettreLow,posSNPLow,lettreUp,posSNPUp,boolRefLow,boolRefUp,bug,erreur,reverseUp,reverseLow,lettreRefUp,lettreRefLow)
	else:	
		return(dicopolUp,dicopolLow,listePolymorphismePosUp,listePolymorphismePosLow)
##############################################################
##############################################################
def MismatchChecker(snpUp,posUp,snpLow,posLow,lettreRefUp,lettreRefLow,lettreUp,lettreLow,boolRefUp,boolRefLow,indel):
	"""In case of divergent main position (case snp whose two paths are mapped ) to define the reference = > check the number of mismatch
( If the number of mismatch is the same in both cases it is the lower lexicographical SNP which is selected for reference .
The Boolean allows to know the reference SNP ) """
	posRef=0
	if int(snpUp[3])>0 and int(snpLow[3])>0:
		nmUp=posUp[int(snpUp[3])]
		nmLow=posLow[int(snpLow[3])]
		if nmUp<nmLow:
			boolRefLow=0
			boolRefUp=1
			lettreRefLow=lettreRefUp
			posRef=int(snpUp[3])
		elif nmUp>nmLow : 
			boolRefLow=1
			boolRefUp=0
			lettreRefUp=lettreRefLow
			posRef=int(snpLow[3])
		elif  nmUp==nmLow:
			if indel==0:
				if lettreUp<lettreLow:
					boolRefLow=0
					boolRefUp=1
					lettreRefLow=lettreRefUp
					posRef=int(snpUp[3])
				elif lettreUp>lettreLow:
					boolRefLow=1
					boolRefUp=0
					lettreRefUp=lettreRefLow
					posRef=int(snpLow[3])
				else :
					if int(snpUp[3])<int(snpLow[3]):
						boolRefLow=0
						boolRefUp=1
						lettreRefLow=lettreRefUp
						posRef=int(snpUp[3])
					else:
						boolRefLow=1
						boolRefUp=0
						lettreRefUp=lettreRefLow
						posRef=int(snpLow[3])
					
			else:
				if int(snpUp[3])<int(snpLow[3]):
					boolRefLow=0
					boolRefUp=1
					lettreRefLow=lettreRefUp
					posRef=int(snpUp[3])
				else:
					boolRefLow=1
					boolRefUp=0
					lettreRefUp=lettreRefLow
					posRef=int(snpLow[3])
	return(boolRefLow,boolRefUp,lettreRefUp,lettreRefLow,posRef)
##############################################################
##############################################################
def GetPolymorphisme(dicoHeader,seq,indel):
	#dicoHeader[key]=[posD,ntUp,ntLow]
	#Forward
	listePos=[]
	listeLettreUp=[]
	listeLettreLow=[]
	#Reverse
	listePosR=[]
	listeLettreUpR=[]
	listeLettreLowR=[]
	tailleSeq=len(seq)
	if indel==0:
		for key,(posD,ntUp,ntLow) in dicoHeader.items():
			listePos.append(posD)
			listePosR.append(tailleSeq-int(posD)+1)
			listeLettreUp.append(ntUp)
			listeLettreUpR.append(ReverseComplement(ntUp))
			listeLettreLow.append(ntLow)
			listeLettreLowR.append(ReverseComplement(ntLow))
		return(listePos,listeLettreUp,listeLettreLow,listePosR,listeLettreUpR,listeLettreLowR)
	else:
		for key,(posD,ind,amb) in dicoHeader.items():
			listePos.append(posD)
			listePosR.append(tailleSeq-int(posD)+1)
			insert=seq[(int(posD-1)-1):(int(posD-1)+int(ind))]
			ntStart=seq[(int(posD-1)-1)]
			ambiguityPos=amb
		return(listePos,listePosR,insert,ntStart,ambiguityPos)
##############################################################
##############################################################
def RemplissageLigneVCF(snpUp,snpLow,lettreLow,positionSnpLow,lettreUp,positionSnpUp,boolRefLow,boolRefUp,table,nbSnp,line,bug,erreur,dmax):
	""" Fills the different fields of vcf based on boolean ( on whether the SNP is identical or not the reference) , if neither is identical to the reference = > we take the SNP comes first in the lexicographical order"""	
	snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp,listeCouvertureUp,listeCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
	snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow,listeCouvertureLow,listeClow,nb_polLow,lnlow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
	seqUp=list(snpUp[9])
	seqLow=list(snpLow[9])
	if int(snpUp[3])>0 and int(snpLow[3])>0:
		table[line][0]=snpLow[2]
		if boolRefLow==1 and boolRefUp==0:	
			table[line][1]=positionSnpLow
			table[line][2]=numSNPLow
			table[line][3]=lettreLow
			table[line][4]=lettreUp
			if snpLow[10]=="*":
				table[line][5]="."
			else:
				table[line][5]=snpLow[10]
		elif boolRefUp==1 and boolRefLow==0:
			table[line][1]=positionSnpUp
			table[line][2]=numSNPUp
			table[line][3]=lettreUp
			table[line][4]=lettreLow
			if snpUp[10]=="*":
				table[line][5]="."
			else:
				table[line][5]=snpUp[10]
		elif boolRefUp==0 and boolRefLow==0:
			if lettreUp<lettreLow:
				table[line][1]=positionSnpUp
				table[line][2]=numSNPUp
				table[line][3]=lettreUp
				table[line][4]=lettreLow
				boolRefUp=1
				boolRefLow=0
				if snpUp[10]=="*":
					table[line][5]="."
				else:
					table[line][5]=snpLow[10]
			else:
				table[line][1]=positionSnpLow
				table[line][2]=numSNPLow
				table[line][3]=lettreLow
				table[line][4]=lettreUp
				boolRefUp=0
				boolRefLow=1
				if snpLow[10]=="*":
					table[line][5]="."
				else:
					table[line][5]=snpLow[10]
		else : 
			#SNP 
			table[line][1]=erreur
			table[line][2]='.'
			table[line][3]='.'
			table[line][4]='.'
			table[line][5]='.'

	elif int(snpUp[3])<=0 and int(snpLow[3])>0:
		table[line][0]=snpLow[2]
		table[line][1]=positionSnpLow
		table[line][2]=numSNPLow
		table[line][3]=lettreLow
		boolRefUp=0
		boolRefLow=1
		if dmax==1:
			table[line][4]=seqUp[(len(seqUp)/2)]
		else:
			table[line][4]='.'
		if snpLow[10]=="*":
			table[line][5]="."
		else:
			table[line][5]=snpLow[10]
	elif int(snpUp[3])>0 and int(snpLow[3])<=0:
		table[line][0]=snpUp[2]
		table[line][1]=positionSnpUp
		table[line][2]=numSNPUp
		table[line][3]=lettreUp
		boolRefUp=1
		boolRefLow=0
		if dmax==1:
			table[line][4]=seqLow[(len(seqLow)/2)]
		else:
			table[line][4]='.'
		if snpUp[10]=="*":
			table[line][5]="."
		else:
			table[line][5]=snpUp[10]
	elif int(snpUp[3])<=0 and int(snpLow[3])<=0:
		if lettreLow<lettreUp:
			table[line][0]="."
			table[line][1]=positionSnpLow
			table[line][2]=numSNPLow
			table[line][3]=lettreLow
			table[line][4]=lettreUp
			boolRefUp=0
			boolRefLow=1
			if snpLow[10]=="*":
				table[line][5]="."
			else:
				table[line][5]=snpLow[10]
		else:
			table[line][0]="."
			table[line][1]=positionSnpUp
			table[line][2]=numSNPUp
			table[line][3]=lettreUp
			table[line][4]=lettreLow
			boolRefUp=1
			boolRefLow=0
			if snpLow[10]=="*":
				table[line][5]="."
			else:
				table[line][5]=snpUp[10]			
	return(table,boolRefUp,boolRefLow)
##############################################################
##############################################################
def RemplissageVCFSNPproches(dicoUp,dicoLow,table,champFilter,dmax,snpUp,snpLow,listePolymorphismePosUp,listePolymorphismePosLow,listePolymorphismePos,line,multi,ok,couvUp,couvLow,listeLettreUp,listeLettreLow,geno,nbGeno,listeCouvGeno):
	info=''
	champAlt=0
	comptPol=0
	seqUp=snpUp[9]
	seqLow=snpLow[9]
	tp="SNP"
	phase=1
	#Variables
	snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp,listeCouvertureUp,listeCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
	snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow,listeCouvertureLow,listeCLow,nb_polLow,lnLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
	#SNPs proches chemins mappés
	if int(snpUp[3])>0 and int(snpLow[3])>0:
		#choix de la position de ref:
		boolRefUp=int(dicoUp[listePolymorphismePosUp[0]][0])
		boolRefLow=int(dicoLow[listePolymorphismePosLow[0]][0])
		lettreUp1=dicoUp[listePolymorphismePosUp[0]][3]
		lettreLow1=dicoLow[listePolymorphismePosLow[0]][3]
		positionSnpUp1=dicoUp[listePolymorphismePosUp[0]][5]
		positionSnpLow1=dicoLow[listePolymorphismePosLow[0]][5]	
		reverseLow=dicoLow[listePolymorphismePosLow[0]][4]
		reverseUp=dicoUp[listePolymorphismePosUp[0]][4]
		for comptPol in range(len(listePolymorphismePos)):
			positionSnpUp=dicoUp[listePolymorphismePosUp[comptPol]][5]
			positionSnpLow=dicoLow[listePolymorphismePosLow[comptPol]][5]
			lettreUp=dicoUp[listePolymorphismePosUp[comptPol]][3]
			lettreRefUp=dicoUp[listePolymorphismePosUp[comptPol]][1]
			lettreLow=dicoLow[listePolymorphismePosLow[comptPol]][3]
			lettreRefLow=dicoLow[listePolymorphismePosLow[comptPol]][1]
			#Remplissage VCF
			if boolRefLow==1 and boolRefUp==0:
				table=GetGenotype(geno,boolRefUp,boolRefLow,table,line,nbGeno,phase,listeCouvGeno)
				table[line][0]=snpLow[2]	
				table[line][1]=positionSnpLow
				table[line][2]=numSNPLow
				table[line][3]=lettreLow
				if snpLow[10]=="*":
					table[line][5]="."
				else:
					table[line][5]=snpLow[10]
				if (int(reverseUp)==-1 and int(reverseLow)==1) or (int(reverseUp)==1 and int(reverseLow)==-1):	
					table[line][4]=ReverseComplement(lettreUp)
				else:
					table[line][4]=lettreUp
			elif boolRefUp==1 and boolRefLow==0:
				table=GetGenotype(geno,boolRefUp,boolRefLow,table,line,nbGeno,phase,listeCouvGeno)
				table[line][0]=snpUp[2]
				table[line][1]=positionSnpUp
				table[line][2]=numSNPUp
				table[line][3]=lettreUp
				if snpUp[10]=="*":
					table[line][5]="."
				else:
					table[line][5]=snpUp[10]
				if (int(reverseUp)==1 and int(reverseLow)==-1) or (int(reverseUp)==-1 and int(reverseLow)==1):
					table[line][4]=ReverseComplement(lettreLow)
				else:
					table[line][4]=lettreLow

			elif boolRefUp==0 and boolRefLow==0:
				if lettreUp1<lettreLow1:
					table[line][0]=snpUp[2]
					table[line][1]=positionSnpUp
					table[line][2]=numSNPUp
					table[line][3]=lettreUp
					table=GetGenotype(geno,1,0,table,line,nbGeno,phase,listeCouvGeno)
					if snpUp[10]=="*":
						table[line][5]="."
					else:
						table[line][5]=snpUp[10]
					if (int(reverseUp)==1 and int(reverseLow)==-1) or (int(reverseUp)==-1 and int(reverseLow)==1):
						table[line][4]=ReverseComplement(lettreLow)
					else:
						table[line][4]=lettreLow
				elif  lettreUp1>lettreLow1:
					table[line][0]=snpLow[2]	
					table[line][1]=positionSnpLow
					table[line][2]=numSNPLow
					table[line][3]=lettreLow
					table=GetGenotype(geno,0,1,table,line,nbGeno,phase,listeCouvGeno)
					if snpLow[10]=="*":
						table[line][5]="."
					else:
						table[line][5]=snpLow[10]
					if (int(reverseUp)==-1 and int(reverseLow)==1) or (int(reverseUp)==1 and int(reverseLow)==-1):	
						table[line][4]=ReverseComplement(lettreUp)
					else:
						table[line][4]=lettreUp
				elif positionSnpUp1<positionSnpLow1:
					table[line][0]=snpUp[2]
					table[line][1]=positionSnpUp
					table[line][2]=numSNPUp
					table[line][3]=lettreUp
					table=GetGenotype(geno,1,0,table,line,nbGeno,phase,listeCouvGeno)
					if snpUp[10]=="*":
						table[line][5]="."
					else:
						table[line][5]=snpUp[10]
					if (int(reverseUp)==1 and int(reverseLow)==-1) or (int(reverseUp)==-1 and int(reverseLow)==1):
						table[line][4]=ReverseComplement(lettreLow)
					else:
						table[line][4]=lettreLow
				elif positionSnpUp1>positionSnpLow1:
					table[line][0]=snpLow[2]	
					table[line][1]=positionSnpLow
					table[line][2]=numSNPLow
					table[line][3]=lettreLow
					table=GetGenotype(geno,0,1,table,line,nbGeno,phase,listeCouvGeno)
					if snpLow[10]=="*":
						table[line][5]="."
					else:
						table[line][5]=snpLow[10]
					if (int(reverseUp)==-1 and int(reverseLow)==1) or (int(reverseUp)==1 and int(reverseLow)==-1):	
						table[line][4]=ReverseComplement(lettreUp)
					else:
						table[line][4]=lettreUp
			if boolRefUp==1:
				info="Ty:"+str(tp)+";"+"Rk:"+str(valRankUp)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftUp)+";"+"UR:"+str(unitigRightUp)+";"+"CL:"+str(contigLeftUp)+";"+"CR:"+str(contigRightUp)+";"+str(couvUp)+";"+"Genome:"+str(lettreRefUp)+";"+"Sd:"+str(reverseUp)
			else:
				info="Ty:"+str(tp)+";"+"Rk:"+str(valRankLow)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftLow)+";"+"UR:"+str(unitigRightLow)+";"+"CL:"+str(contigLeftLow)+";"+"CR:"+str(contigRightLow)+";"+str(couvLow)+";"+"Genome:"+str(lettreRefLow)+";"+"Sd:"+str(reverseLow)
			if dmax==1:
				table[line][7]=info+";"+"DMax:"+str(dmax)
			else:
				table[line][7]=info
			table[line][6]=champFilter
			line+=1
	#SNPs proches : chemins unmapped
	elif int(snpUp[3])<=0 and int(snpLow[3])<=0:
		i=0
		for i in range(len(listePolymorphismePos)):
			positionSnpUp=listePolymorphismePos[i]
			positionSnpLow=listePolymorphismePos[i]
			lettreUp=listeLettreUp[i]
			lettreRefUp=None
			lettreLow=listeLettreLow[i]
			lettreRefLow=None
			table[line][0]=snpUp[2]
			table[line][1]=positionSnpUp
			table[line][2]=numSNPUp
			table[line][3]=lettreUp
			table=GetGenotype(geno,1,0,table,line,nbGeno,phase,listeCouvGeno)
			if snpUp[10]=="*":
				table[line][5]="."
			else:
				table[line][5]=snpUp[10]
			info="Type:"+str(tp)+";"+"Rk:"+str(valRankUp)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftUp)+";"+"UR:"+str(unitigRightUp)+";"+"CL:"+str(contigLeftUp)+";"+"CR:"+str(contigRightUp)+";"+str(couvUp)+";"+"Genome:"+str(lettreRefUp)
			table[line][6]=champFilter
			table[line][7]=info
			line+=1
	
	elif int(snpUp[3])>0 and int(snpLow[3])<=0:
		reverseUp=dicoUp[listePolymorphismePosUp[0]][4]
		for comptPol in range(len(listePolymorphismePos)):
			positionSnpUp=dicoUp[listePolymorphismePosUp[comptPol]][5]
			lettreUp=dicoUp[listePolymorphismePosUp[comptPol]][3]
			lettreRefUp=dicoUp[listePolymorphismePosUp[comptPol]][1]
			table[line][0]=snpUp[2]
			table[line][1]=positionSnpUp
			table[line][2]=numSNPUp
			table[line][3]=lettreUp
			table=GetGenotype(geno,1,0,table,line,nbGeno,phase,listeCouvGeno)
			if snpUp[10]=="*":
				table[line][5]="."
			else:
				table[line][5]=snpUp[10]
			if (int(reverseUp)==-1) and dmax==1:
				table[line][4]=ReverseComplement(seqLow[int(listePolymorphismePosLow[comptPol])-1])
			elif int(reverseUp)==1 and dmax==1:
				table[line][4]=seqLow[int(listePolymorphismePosLow[comptPol])-1]
			else:
				table[line][4]='.'
			info="Ty:"+str(tp)+";"+"Rk:"+str(valRankUp)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftUp)+";"+"UR:"+str(unitigRightUp)+";"+"CL:"+str(contigLeftUp)+";"+"CR:"+str(contigRightUp)+";"+str(couvUp)+";"+"Genome:"+str(lettreRefUp)+";"+"Sd:"+str(reverseUp)
			table[line][7]=info
			table[line][6]=champFilter
			line+=1
	elif int(snpUp[3])<=0 and int(snpLow[3])>0:
		reverseLow=dicoLow[listePolymorphismePosLow[0]][4]
		for comptPol in range(len(listePolymorphismePos)):
			positionSnpLow=dicoLow[listePolymorphismePosLow[comptPol]][5]
			lettreLow=dicoLow[listePolymorphismePosLow[comptPol]][3]
			lettreRefLow=dicoLow[listePolymorphismePosLow[comptPol]][1]
			table[line][0]=snpLow[2]
			table[line][1]=positionSnpLow
			table[line][2]=numSNPLow
			table[line][3]=lettreLow
			table=GetGenotype(geno,0,1,table,line,nbGeno,phase,listeCouvGeno)
			if snpLow[10]=="*":
				table[line][5]="."
			else:
				table[line][5]=snpLow[10]
			if (int(reverseLow)==-1) and dmax==1:
				table[line][4]=ReverseComplement(seqUp[int(listePolymorphismePosUp[comptPol])-1])
			elif int(reverseLow)==1 and dmax==1:
				table[line][4]=seqUp[int(listePolymorphismePosUp[comptPol])-1]
			else:
				table[line][4]='.'
			info="Ty:"+str(tp)+";"+"Rk:"+str(valRankLow)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftLow)+";"+"UR:"+str(unitigRightLow)+";"+"CL:"+str(contigLeftLow)+";"+"CR:"+str(contigRightLow)+";"+str(couvLow)+";"+"Genome:"+str(lettreRefLow)+";"+"Sd:"+str(reverseLow)
			table[line][7]=info
			table[line][6]=champFilter
			line+=1
	return(table,line,boolRefUp,boolRefLow)

##############################################################
##############################################################
def GetGenotype(geno,boolRefUp,boolRefLow,table,line,nbGeno,phase,listeCouvGeno):
	j=0
	for i in range(0,nbGeno):
		gen="G"+str(i+1)
		if boolRefLow==1:
			if "1/1" in geno[gen][0]:
				geno[gen][0]=geno[gen][0].replace("1/1","0/0")
			elif "0/0" in geno[gen][0]:
				geno[gen][0]=geno[gen][0].replace("0/0","1/1")
		if j<=nbGeno-1:
			if phase==1:
				genotype=geno[gen][0].replace("/","|")
			else:
				genotype=geno[gen][0]
			if isinstance(geno[gen][1],list):
				likelihood=str(','.join(geno[gen][1]))
			else:
				likelihood=str(geno[gen][1])
			table[line][9+int(j)]=str(genotype)+":"+str(listeCouvGeno[i])+":"+str(likelihood)
			j+=1
	table[line][8]="GT:DP:PL"
	return(table)
##############################################################
##############################################################



















