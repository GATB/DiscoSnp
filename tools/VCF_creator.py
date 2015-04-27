#!/bin/python
# -*- coding: utf-8 -*-
#*****************************************************************************
#   VCF_Creator: mapping and VCF creation feature in DiscoSnp++
#   Copyright (C) 2015  INRIA
#   Author: C.Riou
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************
import os
import sys
import subprocess
import re
import time
import getopt
from functionsVCF_creator import *
#Default value
nbMismatchBWA=3
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
    opts, args = getopt.getopt(sys.argv[1:],"h:s:n:o:",["help","disco_file","mismatch","output="])
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
        boolmyname=False
        if os.path.isfile(arg):
               fileName = arg
               listNameFile=fileName.split(".")
               if "BWA_OPT" in listNameFile[0]: # when the samfile is created by run_VCF_creator.sh ; it adds BWA_OPT to separate the name of the file and the BWA options
                       boolmyname=True # Boolean to know if the file name contains the option of BWA
                       listName=listNameFile[0].split("BWA_OPT") # parsing of the file name listName[0] corresponds to the name of the discofile listName[1] corresponds to the bwa options
                       listName[1]=listName[1].replace("_", " ")
               samfile=open(fileName,'r')
        else : 
              print "!! No file :" +str(arg)+"!!"
              sys.exit(2)  
    elif opt in ("-n","--mismatch"):
         nbMismatchBWA= arg
    elif opt in ("-o","--output"):
        if arg!=None:
        #if ".vcf" in arg: (we don't care the out file name, user is free to call it foo.hey)
                VCF = open(arg,'w')
        else :
                print "!! No output !!"
                sys.exit(2)      
        #else :
        #    print "...Unknown file extension for the output. Try with a <file>.vcf..."
        #    sys.exit(2)
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
VCF.write('##source=VCF_creator\n')
nbGeno=0
nbSnp=0
nbSnp,nbGeno = Counting(fileName)
if boolmyname:
        VCF.write('##BWA_Options='+str(listName[1])+'\n')
        VCF.write('##SAMPLE=file://'+str(listName[0])+".fa"+'\n')
else:
        VCF.write('##SAMPLE=file://'+str(fileName)+'\n')
VCF.write('##REF=<ID=REF,Number=1,Type=String,Description="Allele of the path Disco aligned with the least mismatches">\n')
VCF.write('##FILTER=<ID=MULTIPLE,Description="Mapping type : PASS or MULTIPLE or .">\n')
VCF.write('##INFO=<ID=Ty,Number=1,Type=String,Description="SNP, INS, DEL or .">\n')
VCF.write('##INFO=<ID=Rk,Number=1,Type=Float,Description="SNP rank">\n')
VCF.write('##INFO=<ID=MULTI,Number=1,Type=String,Description="State of the mapping in BWA : path multiply mapped : True ; else : False">\n')
VCF.write('##INFO=<ID=DT,Number=1,Type=Integer,Description="Mapping distance with reference">\n')
VCF.write('##INFO=<ID=UL,Number=1,Type=Integer,Description="length of the unitig left">\n')
VCF.write('##INFO=<ID=UR,Number=1,Type=Integer,Description="length of the unitig right">\n')
VCF.write('##INFO=<ID=CL,Number=1,Type=Integer,Description="length of the contig left">\n')
VCF.write('##INFO=<ID=CR,Number=1,Type=Integer,Description="length of the contig right">\n')
VCF.write('##INFO=<ID=Genome,Number=1,Type=String,Description="Allele of the reference;for indel reference is <DEL> or <INS>">\n')
VCF.write('##INFO=<ID=Sd,Number=1,Type=Integer,Description="Reverse (-1) or Forward (1) Alignement">\n')
VCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
VCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">\n')
VCF.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled Genotype Likelihoods">\n')
VCF.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Depth of each allele by sample">\n')
table = [0] * 10 # create a 10 cols array

##Create the columns of the VCF File with all the fields + one field by genotypes/samples/individuals
if nbGeno==0: # Without genotypes
    VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')
else:
    i=0
    VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
    for i in range(0,int(nbGeno)):
        nomCol="G"+str(i+1)
        VCF.write(str(nomCol))
        if i<nbGeno-1 :
                VCF.write("\t" )# Add a \t except if this is the last genotype
        if i==int(nbGeno)-1:
            VCF.write("\n")
i=0            
############################################################### 
###VCF_creator NORMAL MODE : take a samfile and create a vcf###
###############################################################
#---------------------------------------------------------------------------------------------------------------------------
#VARIABLES
#--------------------------------------------------------------------------------------------------------------------------- 
#--------------------------------------------------------------------------------------------------------------------------- 
#Upper path
snpUp=None # samline stocks into a list
numSNPUp=None #ID of discoSnp++
unitigLeftUp=None #lenght of the unitig left
unitigRightUp=None #lenght of the unitig right
contigLeftUp=None #lenght of the contig left
contigRightUp=None #lenght of the contig right
valRankUp=None #Rank value
listCoverageUp=None #listCoverageUp=[23,45]
listCUp=None #listCup=["C1","C2"]
nb_polUp=None # number of polymorphisme in the disco path 
posDUp=None #position of the polymorphisme on the disco path
ntUp=None #allele of the upper path from the discosnp++ header
genoUp=None #dictionnary with the information of the genotype dicoGeno[listgeno[0]]=[listgeno[1],listlikelihood]
dicoHeaderUp=None #dictionnary of all the information from the header of discosnp++ : depending on the variant
posUp=None #dictionnary with all the mapping positions associated with their number of mismatches with the reference
boolXAUp=None #boolean to know if a path is multiple mapped in BWA
covUp=None #string with all the coverage by samples
listnucleoUpR=None #list of all the allele for the path in reverse direction
listnucleoUp=None #list of all the allele for the path in forward direction
dicoUp=None # dictionnary with all the informations for close snps : boolean to know if the allele is identical to the reference,position on the variant, reverse nucleotide for the allele, if the path is reverse or not, mapping position
#dicoUp=dicopolUp[listPosR[i]]=[boolRefUp,nucleoRefUp,posCentraleUp[i],listnucleoUpR[i],reverseUp,(int(snpUp[3])+posCentraleUp[i])]
listPolymorphismePosUp=None # list of the position of the allele on the upper path
nucleoUp=None #allele of the upper path in the samfile
positionSnpUp=None #mapping postion of the snp on the upper path
reverseUp=None #"boolean" to know if the path is mapped in forward (1) or reverse (-1) direction 
nucleoRefUp=None #allele at the variant position in the reference genome
boolRefUp=None #boolean to know if the allele on the upper path is identical to the reference 
#---------------------------------------------------------------------------------------------------------------------------
#Lower path
snpLow=None
numSNPLow=None
unitigLeftLow=None
unitigRightLow=None
contigLeftLow=None
contigRightLow=None
valRankLow=None
listCoverageLow=None
listClow=None
nb_polLow=None
posDLow=None
ntLow=None
genoLow=None
dicoHeaderLow=None
posLow=None
boolXALow=None
covLow=None
listnucleoLowR=None
listnucleoLow=None
dicoLow=None
listPolymorphismePosLow=None
nucleoLow=None
positionSnpLow=None
reverseLow=None
nucleoRefLow=None
boolRefLow=None 
#---------------------------------------------------------------------------------------------------------------------------
# Other variables
couple=None #string to know if the variant (couple of path) is mapped at a uniq genomic position for a given distance with the reference
NM=0 #distance with the reference 
rupture=0 #distance with the reference for which the validation algorithm stopped
filterfield=None #string with the result for the filter field of the VCF
listCovGeno=None #list with all the sum of the coverage by sample 
listPolymorphismePos=None #list with all the position for a path
listPosR=None #list with all the reverse position for a path
insert=None #sequence of the insertion + the nucleotide just before
ntStart=None #nucleotide just before the insertion
ambiguity=None #possible ambiguity for the position of the insertion/deletion on the path
key=None # key in dictionnary
numberCloseSNp=0
ntSet="ATCG"
if ".sam" in fileName: #checks if it's a samfile
    while True:
        line1=samfile.readline()
        if not line1: break # end of file
        
        if line1.startswith('@'): continue # we do not read headers
        
        line2=samfile.readline() # read couple of lines
        ##snpUp and snpLow are lists of the line in the samfile file
        
        discoNameUp,snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp,listCoverageUp,listCUp,nb_polUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(line1,0)
        discoNameLow,snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow,listCoverageLow,listCLow,nb_polLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(line2,0)
        #Verifies that the samfile is formatted
        if 0 in [c in "ACGT" for c in snpUp[9]] or len(snpUp)<11:#Checks if it is really a sequence at the 9th place in the line
                print "WARNING wrong format for the variant : "+str(discoNameUp)
                print "...We skipped it..."
                continue
        if 0 in [c in "ACGT" for c in snpLow[9]] or len(snpLow)<11:
                print "WARNING wrong format for the variant : "+str(discoNameLow)
                print "...We skipped it..."
                continue
        if numSNPLow != numSNPUp:
            print "WARNING two consecutive lines do not store the same variant id: "
            print line1
            print line2
            sys.exit(1)
        table = [0] * 10 # create a 10 cols array;init it every line   
        #Information on coverage by dataset
        covUp,covLow,listCovGeno=GetCoverage(listCUp,listCLow,listCoverageUp,listCoverageLow)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------                
        #Variables
        comptPol=0     # number of SNPs in case of close SNPs - useless for indels
        info=None      # info vcf field
        ok=None        # distance for which the SNP is mapped, -1 if not mapped or if multiple mapped
        indel=False    # boolean to know if it is an indel
        phased=False    # am I phased?
        filterField='.' # init the vcf field filter
        posUp,posLow,boolXAUp,boolXALow = GetCouple(snpUp,snpLow) # get all the positions of mapping for one variant with the associated number of mapping errors
        seqUp=snpUp[9]   # sequences
        seqLow=snpLow[9] # sequences
        
    
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------                    
        #Validation of a couple of SNPs SEE DOC.
        NM=0
        rupture=0
        for NM in range(0,int(nbMismatchBWA)+1): 
            couple=ValidationSNP(snpLow,posLow,snpUp,posUp,NM) 
            if couple== "ok" or couple == "multiple" or couple=="unmapped":
                rupture=NM
                break
        #VCF champs Filter
        if couple=="ok":
            filterField="PASS"            
            ok=str(NM)
        
        elif couple == "unmapped":
            filterField="."
            ok="."
        elif couple == "multiple":
            filterField="MULTIPLE"
            ok=-1
        else : 
            filterField="bug"

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

        if "SNP" in snpUp[0] :
#---------------------------------------------------------------------------------------------------------------------------
            indel=False
            listPolymorphismePos=[]
            listPos,listnucleoUp,listnucleoLow,listPosR,listnucleoUpR,listnucleoLowR=GetPolymorphisme(dicoHeaderUp,snpUp[9],indel,False)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
            ##one SNP
            if int(nb_polUp)==1:
                tp="SNP"
                nucleoLow,positionSnpLow,nucleoUp,positionSnpUp,boolRefLow,boolRefUp,reverseUp,reverseLow,nucleoRefUp,nucleoRefLow = RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel)
                #Creation VCF
                table=fillVCFSimpleSnp(snpUp,snpLow,nucleoLow,positionSnpLow,nucleoUp,positionSnpUp,boolRefLow,boolRefUp,table,nbSnp,filterField,ok,tp,phased,listCovGeno,nucleoRefUp,nucleoRefLow,reverseUp,reverseLow,genoUp,nbGeno,covUp,covLow)
                printOneline(table,VCF)
                continue
    
            else:
                numberCloseSNp+=int(nb_polUp)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
                ##Close SNPs
                indel=False
                dicoUp={}
                dicoLow={}
                dicoUp,dicoLow,listPolymorphismePosUp,listPolymorphismePosLow=RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel)
                # this function comptutes the VCF and prints it!!
                printVCFSNPclose(dicoUp,dicoLow,table,filterField,snpUp,snpLow,listPolymorphismePosUp,listPolymorphismePosLow,listPolymorphismePos,ok,covUp,covLow,listnucleoUp,listnucleoLow,genoUp,nbGeno,listCovGeno,VCF) 
                continue # 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
        ##Case of INDEL                   
        if "INDEL" in snpUp[0] :
            indel=True
            seqInsert=0
            seqUp=snpUp[9]
            seqLow=snpLow[9]
            # To know in which sequence the insertion is and to get it
            #Reverse the sequence to get the good mapped insert if it needs to
            seqUp,seqLow=GetSequence(snpUp,snpLow)
            if len(seqUp)<len(seqLow):
                seq=seqLow
            else:
                seq=seqUp
#---------------------------------------------------------------------------------------------------------------------------
            #Get from the dicoHeader[key]=[posD,ind,amb]: Position of the insertion ; the insertion with the nucleotide just before ; the nucleotide just before ; the possible ambiguity for the position of the insertion 
            listPos,listPosR,insert,ntStart,ambiguity=GetPolymorphisme(dicoHeaderUp,seq,indel,False)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
            #Get the positons of the variant by taking into account the shift of mapping
            nucleoLow,positionSnpLow,nucleoUp,positionSnpUp,boolRefLow,boolRefUp,reverseUp,reverseLow,nucleoRefUp,nucleoRefLow= RecupPosSNP(snpUp,snpLow,posUp,posLow,nb_polUp,nb_polLow,dicoHeaderUp,indel)
            #Checks the strand (forward or reverse) to have the right sequence of insert
            if boolRefUp==True and reverseUp==-1:
                insert=ReverseSeq(insert)
                ntStart=ReverseComplement(ntStart)
            elif boolRefLow==True and reverseLow==-1:
                insert=ReverseSeq(insert)
                ntStart=ReverseComplement(ntStart)
            #Checks if the insert correpond to the upper path or to the lower path
            if len(seqUp)<len(seqLow):
                nucleoLow=insert
                nucleoUp=ntStart
            else:
                nucleoUp=insert
                nucleoLow=ntStart
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
            ## Fill the VCF if the upper path is considered as the reference
            if boolRefUp==True:
                if len(nucleoUp)==len(insert):
                    nucleoRefUp="."
                    tp="INS"
                else:
                    nucleoRefUp="."
                    tp="DEL"
                table=FillVCF(table,numSNPUp,snpUp[2],int(positionSnpUp)-1,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,genoUp,nbGeno,phased,listCovGeno,boolRefLow)
            ## Fill the VCF if the lower path is considered as the reference
            elif boolRefLow==True:
                if len(nucleoLow)==len(insert):
                    nucleoRefLow="."
                    tp="INS"
                else:
                    nucleoRefLow="."
                    tp="DEL"
                table=FillVCF(table,numSNPLow,snpLow[2],int(positionSnpLow)-1,nucleoLow,nucleoUp,snpLow[10],filterField,tp,valRankLow,ok,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,covLow,nucleoRefLow,reverseLow,genoLow,nbGeno,phased,listCovGeno,boolRefLow)
            else:
                  if len(nucleoLow)==len(insert):
                    nucleoRefLow="."
                    tp="DEL"
                    nucleoRefUp="."
                    
                  elif len(nucleoUp)==len(insert):
                    nucleoRefUp="."
                    tp="INS"
                    nucleoRefLow="."  
                  table=FillVCF(table,numSNPUp,discoNameUp,int(positionSnpUp)-1,nucleoUp,nucleoLow,snpUp[10],filterField,tp,valRankUp,ok,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,nucleoRefUp,reverseUp,genoUp,nbGeno,phased,listCovGeno,boolRefLow)      
            printOneline(table,VCF)
            continue
            
#---------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
##################################################################### 
###VCF_creator GHOST MODE : take a fasta in input and create a vcf###
#####################################################################   
else:
    while True:
        line1=samfile.readline()
        if not line1: break # end of file
        seq1=samfile.readline() # read the seq associate to the SNP
        line2=samfile.readline() # read a couple of line
        seq2=samfile.readline()
        line1=line1.rstrip('\n')
        line1=line1.strip('>')
        line2=line2.rstrip('\n')
        line2=line2.strip('>')
        table = [0] * 10
        posUnmappedUp=None
        pos=None
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------        
        #Variables
        comptPol=0     # number of SNPs in case of close SNPs - useless for indels
        
        discoNameUp,snpUp,numSNPUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,valRankUp,listCoverageUp,listCUp,nb_polUp,posDUp,ntUp,ntLow,genoUp,dicoHeaderUp=ParsingDiscoSNP(line1,0)
        discoNameLow,snpLow,numSNPLow,unitigLeftLow,unitigRightLow,contigLeftLow,contigRightLow,valRankLow,listCoverageLow,listCLow,nb_polLow,posDLow,ntUp,ntLow,genoLow,dicoHeaderLow=ParsingDiscoSNP(line2,0)
        posUnmappedUp=CheckContigUnitig(unitigLeftUp,contigLeftUp)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------        
        if numSNPLow != numSNPUp:
            print "WARNING two consecutive lines do not store the same variant id: "
            print line1
            print line2
            sys.exit(1)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------            
        #Information on coverage by dataset
        covUp,covLow,listCovGeno=GetCoverage(listCUp,listCLow,listCoverageUp,listCoverageLow)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
        ##one SNP
        if "SNP" in line1:
            tp="SNP"
            indel=False
            listPolymorphismePos=[]
            #Gets the position and the nucleotide of the variants by header parsing
            if (int(nb_polLow)>=2) or (int(nb_polUp)>=2):
                listPolymorphismePos,listnucleoUp,listnucleoLow,listPosR,listnucleoUpR,listnucleoLowR=GetPolymorphisme(dicoHeaderUp,seq1,indel,False)
	#dicoHeader[key]=[posD,ntUp,ntLow]
	    if len(listPolymorphismePos)==0:
                ntLow=dicoHeaderUp["P_1"][2]
                ntUp=dicoHeaderUp["P_1"][1]
                phased=False
                pos=(int(dicoHeaderUp["P_1"][0])+int(posUnmappedUp))
                PrintVCFGhost(table,numSNPUp,discoNameUp,pos,tp,valRankUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,ntUp,ntLow,genoUp,nbGeno,phased,listCovGeno,VCF)
                continue
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
            ##Close SNPs
	    else:
                phased=True
                ID=1
                for comptPol in range(0,len(listPolymorphismePos)):
                        key="P_"+str(comptPol+1)
                        ntLow=dicoHeaderUp[key][2]
                        ntUp=dicoHeaderUp[key][1]
                        pos=(int(dicoHeaderUp[key][0])+int(posUnmappedUp))
                        PrintVCFGhost(table,str(numSNPUp)+"_"+str(ID),discoNameUp,pos,tp,valRankUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,ntUp,ntLow,genoUp,nbGeno,phased,listCovGeno,VCF)
                        ID+=1
                continue
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
        ##Case of INDEL
        elif "INDEL" in line1:
            indel=True
            phased=False
            if len(seq1)<len(seq2):
                seq=seq1
            else:
                seq=seq2
            listPos,listPosR,insert,ntStart,ambiguity=GetPolymorphisme(dicoHeaderUp,seq,indel,True)                
            if seq==seq1:
                ntLow=ntStart
                ntUp=insert
            else:
                ntUp=ntStart
                ntLow=insert
            pos=(int(dicoHeaderUp["P_1"][0])+int(posUnmappedUp))
            PrintVCFGhost(table,numSNPUp,discoNameUp,pos,tp,valRankUp,unitigLeftUp,unitigRightUp,contigLeftUp,contigRightUp,covUp,ntUp,ntLow,genoUp,nbGeno,phased,listCovGeno,VCF)
            continue     
        	
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

VCF.close()
samfile.close()




