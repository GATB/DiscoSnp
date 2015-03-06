#!/bin/python
# -*- coding: utf-8 -*-
###############################################
#
#
#
import os
import sys
import subprocess
import re
import time
import getopt
from functionsVCF_creator import *
#Help
def usage():
	usage= """
	################################
		Run VCF_creator
	################################
	
	-h --help : print this message
	-s --sam_file : <file>.sam of the alignment
	-n --mismatch : number of differences allowed in the mapper (BWA)
	-o --output : vcf file 
	
	"""
	print usage

try:
	opts, args = getopt.getopt(sys.argv[1:],"h:s:n:o",["help","disco_file","mismatch","output="])
	if not opts:
		usage()
		sys.exit(2)
except getopt.GetoptError, e:
	print e
	usage()
	sys.exit(2)
for opt, arg in opts : 
	if opt in ("-h", "--help"):
		usage()
		sys.exit(2)
	elif opt in ("-s","--sam_file"):
		fichier = arg
		if ".sam" in fichier:
			samfile=open(fichier,'r')
		else:
			print "...Unknown file extension for the input. Try with a <file>.sam..."
			sys.exit(2)
	elif opt in ("-n","--mismatch"):
		 nbMismatchBWA= arg
	elif opt in ("-o","--output"):
		if ".vcf" in arg:
			VCF = open(arg,'w')
		else :
			print "...Unknown file extension for the output. Try with a <file>.vcf..."
			sys.exit(2)
	else:
		print("Unkwnown option {} ".format(opt))
		usage()
		sys.exit(2)


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
###Header of the VCF file 
today=time.localtime()
date=str(today.tm_year)+str(today.tm_mon)+str(today.tm_mday)
VCF.write('##fileformat=VCFv4.1\n')
VCF.write('##filedate='+str(date)+'\n')
VCF.write('##REF=<ID=REF,Number=1,Type=String,Description="Allele of the path Disco aligned with the least mismatches ">\n')
VCF.write('##ALT=<ID=ALT,Number=1,Type=String,Description="Allele of the other path ">\n')
VCF.write('##FILTER=<ID=MULTIPLE,Number=1,Type=String,Description="Mapping type : PASS or MULTIPLE or \'.\' ">\n')
VCF.write('##INFO=<ID=Ty,Number=1,Type=Float,Description="SNP, INS, DEL or "."">\n')
VCF.write('##INFO=<ID=Rk,Number=1,Type=Float,Description="SNP rank">\n')
VCF.write('##INFO=<ID=DT,Number=1,Type=Integer,Description="Mapping distance with reference">\n')
VCF.write('##INFO=<ID=UL,Number=1,Type=Integer,Description="Lenght of the unitig left">\n')
VCF.write('##INFO=<ID=UR,Number=1,Type=Integer,Description="Lenght of the unitig right">\n')
VCF.write('##INFO=<ID=CL,Number=1,Type=Integer,Description="Lenght of the contig left">\n')
VCF.write('##INFO=<ID=CR,Number=1,Type=Integer,Description="Lenght of the contig right">\n')
VCF.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth accross samples (average)">\n')
VCF.write('##INFO=<ID=Genome,Number=1,Type=String,Description="Allele of the reference;for indel reference is <DEL> or <INS>">\n')
VCF.write('##INFO=<ID=Sd,Number=1,Type=Integer,Description="Reverse (1) or Forward (0) Alignement">\n')
nbGeno=0
nbSnp,nbGeno = Comptage(fichier)
nbCol=9+int(nbGeno)
table = [[0] * int(nbCol) for _ in range(int(nbSnp))]
if nbGeno==0:
	VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')
else:
	VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
	for i in range(0,int(nbGeno)):
		nomCol="G"+str(i+1)
		VCF.write(str(nomCol)+"\t")
		if i==int(nbGeno)-1:
			VCF.write(str(nomCol)+"\n")			
snpAvant=None
line=0
mul=0
nbOk=0
nbmulti=0
nbrupture1=0
nbrupture2=0
nbrupture3=0
nbrupturedmax=0
nbrupture1=0
nbrupture0=0
pb=0
nbunmapped=0
if ".sam" in fichier:
	for snp in samfile :
		if '@' != snp[0] :
			if snpAvant is not None:
				table[line][8]="."
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
				#Check if the snps are from the same path
				numSNP1=ParsingDiscoSNP(snp,1)
				numSNP2=ParsingDiscoSNP(snpAvant,1)			
				getcouple=isCouple(numSNP1,numSNP2)
				if getcouple==1:
					dicoHeaderUp={}
					dicoHeaderLow={}
					#Get DiscoSnp Variables
					if "higher" in snp:
						snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp,listeCouvertureUp,listeCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snp,0)
					else:
						snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow,listeCouvertureLow,listeCLow,nb_polLow,lnLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snp,0)
								
					if "higher" in snpAvant:					
						snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp,listeCouvertureUp,listeCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpAvant,0)
					else:
						snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow,listeCouvertureLow,listeCLow,nb_polLow,lnLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpAvant,0)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------				
					#Information on coverage by dataset
					couvUp,couvLow=GetCoverage(listeCUp,listeCLow,listeCouvertureUp,listeCouvertureLow)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------				
					#Variables
					comptPol=0
					info=None
					multi=None
					ok=None
					dmax=0
					indel=0
					phase=0
					champFilter='.'
					posUp,posLow,snpLow,snpUp,boolMapUp,boolMapLow,boolXAUp,boolXALow = GetCouple(snpUp,snpLow)
					seqUp=snpUp[9]
					seqLow=snpLow[9]
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------					
					#Validation of a couple of SNPs 
					NM=0
					for NM in range(0,int(nbMismatchBWA)+1):
						couple=ValidationSNP(snpLow,posLow,snpUp,posUp,NM)
						if couple== "ok" or couple == "multiple":
							rupture=NM
							break
					if rupture==nbMismatchBWA:
						if snpLow[1]=="4":
							dmax=1
						elif snpUp[1]=="4":
							dmax=1
						else:
							dmax==0
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
					#VCF champ INFO Multi
					if boolXAUp==1 and boolXALow==1:
						multi="multi"
						mul+=1
					elif boolXAUp==0 and boolXALow==0:
						multi="none"
					elif (boolXAUp==1 and boolXALow==0) or  (boolXAUp==0 and boolXALow==1):
						multi="one"
						mul+=1
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
					#VCF champs Filter
					if couple=="ok":
						champFilter="PASS"			
						ok=str(NM)
					
					elif couple == "unmapped":
						champFilter="."
						ok=-1
					elif couple == "multiple":
						champFilter="MULTIPLE"
						ok=-1
					else : 
						champFilter="probleme...."
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
					#Ecriture dans le vcf
					if "SNP" in snpUp[0] and "SNP" in snpLow[0]:
#---------------------------------------------------------------------------------------------------------------------------
						indel=0
						#Vérification du nombre de polymorphisme
						listePolymorphismePos=[]
						#Récupération des positions de mapping associées à leur nombre de mutation
						if (int(nb_polLow)>=2) or (int(nb_polUp)>=2):
							listePolymorphismePos,listeLettreUp,listeLettreLow,listePosR,listeLettreUpR,listeLettreLowR=GetPolymorphisme(dicoHeaderUp,seqUp,indel)
						#Remplissage du vcf selon 2 cas :
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
						##SNP simple
						if len(listePolymorphismePos)==0:
							tp="SNP"
							table[line][6]=champFilter
							posModif=len(snpUp[9])/2
							lettreLow,positionSnpLow,lettreUp,positionSnpUp,boolRefLow,boolRefUp,bug,erreur,reverseUp,reverseLow,lettreRefUp,lettreRefLow = RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel)
							#Création VCF
							table,boolRefUp,boolRefLow=RemplissageLigneVCF(snpUp,snpLow,lettreLow,positionSnpLow,lettreUp,positionSnpUp,boolRefLow,boolRefUp,table,nbSnp,line,bug,erreur,dmax)
							table=GetGenotype(genoUp,boolRefUp,boolRefLow,table,line,nbGeno,phase)
							#Remplissage du champ info en fonction de la référence
							if boolRefUp==1:
								info="Ty:"+str(tp)+";"+"Rk:"+str(valRankUp)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftUp)+";"+"UR:"+str(unitigRightUp)+";"+"CL:"+str(contigLeftUp)+";"+"CR:"+str(contigRightUp)+";"+str(couvUp)+";"+"Genome:"+str(lettreRefUp)+";"+"Sd:"+str(reverseUp)
							else:
								info="Ty:"+str(tp)+";"+"Rk:"+str(valRankLow)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftLow)+";"+"UR:"+str(unitigRightLow)+";"+"CL:"+str(contigLeftLow)+";"+"CR:"+str(contigRightLow)+";"+str(couvLow)+";"+"Genome:"+str(lettreRefLow)+";"+"Sd:"+str(reverseLow)
							if dmax==1:
								table[line][7]=info+";"+"DMax:"+str(dmax)
							else:
								table[line][7]=info	
				
						else: 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
							##SNPs proches
							indel=0
							dicoUp={}
							dicoLow={}
							posModif=None
							dicoUp,dicoLow,listePolymorphismePosUp,listePolymorphismePosLow=RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel)
							table,line,boolRefUp,boolRefLow=RemplissageVCFSNPproches(dicoUp,dicoLow,table,champFilter,dmax,snpUp,snpLow,listePolymorphismePosUp,listePolymorphismePosLow,listePolymorphismePos,line,multi,ok,couvUp,couvLow,listeLettreUp,listeLettreLow,genoUp,nbGeno)
							snpAvant=snp
							continue
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
					#Cas des Indels					
					if "INDEL" in snpUp[0] and "INDEL" in snpLow[0]:
						indel=1
						table[line][6]=champFilter
						seqInsert=0
						seqUp=snpUp[9]
						seqLow=snpLow[9]
						if len(seqUp)<len(seqLow):
							seq=seqLow
						else:
							seq=seqUp
						listePos,listePosR,insert,ntStart,ambiguity=GetPolymorphisme(dicoHeaderUp,seq,indel)
						snpUp,numSNPUp,unitigLeftUP,unitigRightUp,contigLeftUp,contigRightUp,valRankUp,listeCouvertureUp,listeCUp,nb_polUp,lnUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(snpUp,0)
						snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow,listeCouvertureLow,listeCLow,nb_polLow,lnLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(snpLow,0)
						lettreLow,positionSnpLow,lettreUp,positionSnpUp,boolRefLow,boolRefUp,bug,erreur,reverseUp,reverseLow,lettreRefUp,lettreRefLow= RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel)
						if len(seqUp)<len(seqLow) and reverseUp==0:
							lettreLow=insert
							lettreUp=ntStart
						else:
							lettreUp=insert
							lettreLow=ntStart
						if boolRefUp==1:
							if len(lettreUp)==len(insert):
								lettreRefUp="."
								tp="INS"
							else:
								lettreRefUp="."
								tp="DEL"
							table[line][0]=snpUp[2]					
							table[line][1]=int(positionSnpUp)-1
							table[line][2]=numSNPUp
							table[line][3]=lettreUp
							table[line][4]=lettreLow
							table=GetGenotype(genoUp,1,0,table,line,nbGeno,phase)
							if snpUp[10]=="*":
								table[line][5]="."
							else:
								table[line][5]=snpUp[10]
							info="Ty:"+str(tp)+";"+"Rk:"+str(valRankUp)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftUp)+";"+"UR:"+str(unitigRightUp)+";"+"CL:"+str(contigLeftUp)+";"+"CR:"+str(contigRightUp)+";"+str(couvUp)+";"+"Genome:"+str(lettreRefUp)+";"+"Sd:"+str(reverseUp)
						elif boolRefLow==1:
							if len(lettreLow)==len(insert):
								lettreRefLow="."
								tp="INS"
							else:
								lettreRefLow="."
								tp="DEL"
							table[line][0]=snpLow[2]	
							table[line][1]=int(positionSnpLow)-1
							table[line][2]=numSNPLow
							table[line][3]=lettreLow
							table[line][4]=lettreUp
							table=GetGenotype(genoUp,0,1,table,line,nbGeno,phase)
							if snpLow[10]=="*":
								table[line][5]="."
							else:
								table[line][5]=snpLow[10]
							info="Ty:"+str(tp)+";"+"Rk:"+str(valRankLow)+";"+"MULTI:"+str(multi)+";"+"DT:"+str(ok)+";"+"UL:"+str(unitigLeftLow)+";"+"UR:"+str(unitigRightLow)+";"+"CL:"+str(contigLeftLow)+";"+"CR:"+str(contigRightLow)+";"+str(couvLow)+";"+"Genome:"+str(lettreRefLow)+";"+"Sd:"+str(reverseLow)
						else:
							if lettreUp==refLettre or lettreUp==refLettreR:
								lettreRefUp="<DEL>"
							else:
								lettreRefUp="<INS>"
							table[line][0]=snpUp[2]					
							table[line][1]=pos
							table[line][2]=numSNPUp
							table[line][3]=lettreUp
							table[line][4]=lettreLow
							table=GetGenotype(genoUp,1,0,table,line,nbGeno,phase)
							table[line][5]=snpUp[10]
						table[line][7]=info
					line+=1
			snpAvant=snp	
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#Ecriture fichier
for line in table :
	for element in line:			
		VCF.write(str(element))
		VCF.write("\t")
	VCF.write('\n')
VCF.close()
samfile.close()

