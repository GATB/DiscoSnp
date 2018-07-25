#!/bin/python
# -*- coding: utf-8 -*-
###############################################
import os
import sys
import subprocess
import re
import time
from ClassVCF_creator import *
#############################################################################################
#Function_________________________________________________________________________________________________________
#INDEX_____________________________________________________________________________________________________________
#      InitVariant(line1,line2):"""Initialization of the variant by taking into accoutn its type"""
#      MappingTreatement(variant_object,vcf_field_object,nbGeno):
#      UnmappedTreatement(variant_object,vcf_field_object,nbGeno,seq1,seq2):"""Fills VCFfile in ghost mode (from a fasta file)"""
#      CounterGenotype(fileToRead)
#      CheckAtDistanceXBestHits(upper_path,lower_path):"""Prediction validation : check if the couple is validated with only one mapping position """
#      PrintVCFHeader(VCF,listName,fileName,boolmyname):    
#############################################################################################
def InitVariant(line1,line2,fileName,dicoIndex):
        """Initialization of the variant by taking into account its type"""
        #Object Creation    
        if "SNP" in line1 and "|nb_pol_1|" in line1:
                variant_object=SNP(line1,line2)
        elif "SNP" in line1 and "|nb_pol_1|" not in line1:
                variant_object=SNPSCLOSE(line1,line2)
        elif "INDEL" in line1:
                #Supression of indel when the ambiguity is greater than 20
                #if int(line1.split("\t")[0].split("|")[1].split("_")[3])>=20:
                #        return 1,1
                variant_object=INDEL(line1,line2)
        else :
                print("!!!!Undefined Variant!!!!")
                return (1,1)                
        variant_object.setDicoIndex(dicoIndex)                
        #VCF object Creation and filling variant's attribut   
        vcf_field_object=VCFFIELD()
        variant_object.FillInformationFromHeader(vcf_field_object)
        return (variant_object, vcf_field_object)   
        
#############################################################################################
#############################################################################################
def MappingTreatement(variant_object,vcf_field_object,nbGeno):
        """Treatement of the mapping informations and filling of the vcf object """
        table = [0] * 10 #Creates a 10 cols array
        #Fills information of the variant object with the informations of the discosnp++ header
        variant_object.RetrievePolymorphismFromHeader() 
        #Gets the coverage and quality for each path
        variant_object.upper_path.RetrieveCoverage(variant_object.dicoIndex)
        variant_object.lower_path.RetrieveCoverage(variant_object.dicoIndex)
        variant_object.upper_path.RetrieveQualityFQ(variant_object.dicoIndex)
        variant_object.lower_path.RetrieveQualityFQ(variant_object.dicoIndex)
        #Gives the position corrected for each path by taking into account the shift, deletion, insertion ....
        dicoCloseUp=variant_object.upper_path.CheckPosVariantFromRef(vcf_field_object)
        dicoCloseLow=variant_object.lower_path.CheckPosVariantFromRef(vcf_field_object)
        #Checks if the variant has close SNPs
        if int(variant_object.nb_pol)>1:#Close SNPs
                variant_object.RetrieveDicoClose(dicoCloseUp,dicoCloseLow)#fills the dictionnary with all informations for each snps of the path
                table=variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File
        else:#Indel simple SNP
                variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File  
        #Defines the mapping position for the couple
        variant_object.RetrieveMappingPositionCouple()
        #Checks if the couple is validated with only one mapping position       
        vcf_field_object.RetrieveFilterField(CheckAtDistanceXBestHits(variant_object.upper_path,variant_object.lower_path))
        #Defines the genotype for the couple
        variant_object.RetrieveGenotypes(nbGeno,vcf_field_object)
        #Defines variant with multiple mapping : return the XA tag in the vcf file : case of multiply mapped variant
        variant_object.upper_path.RetrieveXA(vcf_field_object)
        variant_object.lower_path.RetrieveXA(vcf_field_object)
        return(table)
        
#############################################################################################
#############################################################################################
def UnmappedTreatement(variant_object,vcf_field_object,nbGeno,seq1,seq2):
        """Treatement of the couple of path without alignment"""
        
        seq1=seq1.rstrip('\n')
        seq2=seq2.rstrip('\n')
        #Gets the sequence of each path
        variant_object.upper_path.RetrieveSeq(seq1)
        variant_object.lower_path.RetrieveSeq(seq2)
        #Fills information of the variant object with the informations of the discosnp++ header        
        variant_object.RetrievePolymorphismFromHeader()
        #Gets the coverage for each path        
        variant_object.upper_path.RetrieveCoverage(variant_object.dicoIndex)
        variant_object.lower_path.RetrieveCoverage(variant_object.dicoIndex)
        variant_object.upper_path.RetrieveQualityFQ(variant_object.dicoIndex)
        variant_object.lower_path.RetrieveQualityFQ(variant_object.dicoIndex)
        dicoCloseUp=variant_object.upper_path.CheckPosVariantFromRef(vcf_field_object)
        dicoCloseLow=variant_object.lower_path.CheckPosVariantFromRef(vcf_field_object)
        #Checks if the variant has close SNPs        
        if int(variant_object.nb_pol)>1:#Close SNPs
                variant_object.RetrieveDicoClose(dicoCloseUp,dicoCloseLow)#fills the dictionnary with all informations for each snps of the path
                table=variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File
        else:#Indel simple SNP
                variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File
                table = [0]*10
        #Defines the mapping position for the couple        
        variant_object.RetrieveMappingPositionCouple()
        #Defines the genotype for the couple
        variant_object.RetrieveGenotypes(nbGeno,vcf_field_object)
        return(table)         
#############################################################################################
#############################################################################################
def CounterGenotype(fileName):
        samfile=open(fileName,'r')
        nbGeno=0
        while True:
                line=samfile.readline()
                if not line: break #End of file
                if line.startswith('@'): continue #We do not read headers
                #>SNP_higher_path_3|P_1:30_C/G|high|nb_pol_1|left_unitig_length_86|right_unitig_length_261|left_contig_length_166|right_contig_length_761|C1_124|C2_0|Q1_0|Q2_0|G1_0/0:10,378,2484|G2_1/1:2684,408,10|rank_1
                if nbGeno==0:
                        line=line.rstrip('\r').rstrip('\n').split('\t')
                        for i in line:
                                if 'SNP' or 'INDEL' in i:
                                        nomDisco=i.split('|')
                                        for k in nomDisco:
                                                if k[0]=='G':
                                                        nbGeno+=1
                                        samfile.close()
                                        return(nbGeno)
#############################################################################################
#############################################################################################
def GetIndex(fileName):
       stream_file=open(fileName,'r')
       while True:                
                line=stream_file.readline()
                if not line: break #End of file
                if line.startswith('@'): continue #We do not read headers
                if ".fa" in fileName:
                        line=line.strip('>')
                        listLine=line.split("|")
                elif ".sam" in fileName:
                        listLine=line.split("\t")[0].split("|")
                #Init dictionnary
                dicoIndex={}
                if "C1_" in line:
                        dicoIndex["C"]=[]
                if "G1_" in line:
                        dicoIndex["G"]=[]
                if "Q1_" in line:
                        dicoIndex["Q"]=[]
                if "unitig" in line:
                       dicoIndex["unitig"]=[]
                if "contig" in line :
                       dicoIndex["contig"]=[]     
                for i in range(len(listLine)):
                        if 'P_1' in listLine[i]:#P_1:30_A/G => {'P_1': ['30', 'A', 'G']} or P_1:30_A/G,P_2:31_G/A
                                dicoIndex["P_"]=int(i)                         
                        elif "unitig" in listLine[i]:                                
                                if "left" in listLine[i]:
                                        dicoIndex["unitig"].append(int(i))
                                if "right" in listLine[i]:
                                        dicoIndex["unitig"].append(int(i))                                     
                        elif "contig" in listLine[i]:
                                if "left" in listLine[i]:
                                       dicoIndex["contig"].append(int(i)) 
                                if "right" in listLine[i]:
                                       dicoIndex["contig"].append(int(i)) 
                        elif "rank" in listLine[i]:
                                dicoIndex["rank"]=int(i)
                        elif "nb_pol" in listLine[i]:
                                dicoIndex["nb_pol"]=int(i)                                                             
                        elif "G" in listLine[i]: #Gets the genotype and likelihood by samples
                               matchG=re.match(r'^G',listLine[i])# finds the genotype in the item of the dicoSnp++ header
                               if matchG:
                                       dicoIndex["G"].append(int(i))
                        elif "C" in listLine[i]:                                
                                matchC=re.match(r'^C',listLine[i])
                                if matchC:
                                       dicoIndex["C"].append(int(i))
                                       
                        elif "Q" in listLine[i]:
                              matchQ=re.match(r'^Q',listLine[i])
                              if matchQ:
                                dicoIndex["Q"].append(int(i))  
                break
       stream_file.close()                              
       return(dicoIndex)                     

#############################################################################################
#############################################################################################
def CheckAtDistanceXBestHits(upper_path,lower_path):
        """Prediction validation : checks if the couple is validated with only one mapping position """
        
        posUp=upper_path.dicoMappingPos
        posLow=lower_path.dicoMappingPos
        # get the best mapping distance for upper path 
        best_up=1024
        if int(upper_path.mappingPosition)==0 and int(lower_path.mappingPosition)==0:#Checks if paths are unmappped
                return(".")
        for position,(nbMismatch,cigarcode) in posUp.items(): 
                if nbMismatch<best_up:
                        best_up=nbMismatch

        # get the best mapping distance for lower path 
        best_low=1024
        for position,(nbMismatch,cigarcode) in posLow.items(): 
                if nbMismatch<best_low:
                        best_low=nbMismatch
    
        # get the union of the mapping position at the best mapping positions
        position_set = set()
        for position,(nbMismatch,cigarcode) in posUp.items():
                if nbMismatch == best_up:
                        position_set.add(position)
                if len(position_set) > 1: 
                        return("MULTIPLE")

        for position,(nbMismatch,cigarcode) in posLow.items():
                if nbMismatch == best_low:
                        position_set.add(position)
                if len(position_set) > 1: 
                        return("MULTIPLE")
    
        if len(position_set) == 1: 
                return("PASS")
    
        return(".")

#############################################################################################
#############################################################################################
def PrintVCFHeader(VCF,listName,fileName,boolmyname):
        ###Header of the VCF file 
        samfile=open(fileName,'r')
        boolQuality=False
        while True:                
                line=samfile.readline()
                if not line: break #End of file
                if line.startswith('@'): continue #We do not read headers
                if "Q1_" in line : 
                        boolQuality=True
                        break
                else: break 
        today=time.localtime()
        date=str(today.tm_year)+str(today.tm_mon)+str(today.tm_mday)
        VCF.write('##fileformat=VCFv4.1\n')
        VCF.write('##filedate='+str(date)+'\n')
        VCF.write('##source=VCF_creator\n')
        nbGeno=0
        nbSnp=0
        nbGeno = CounterGenotype(fileName)
        VCF.write('##SAMPLE=file://'+str(fileName)+'\n')
        VCF.write('##REF=<ID=REF,Number=1,Type=String,Description="Allele of the path Disco aligned with the least mismatches">\n')
        VCF.write('##FILTER=<ID=MULTIPLE,Description="Mapping type : PASS or MULTIPLE or .">\n')
        VCF.write('##INFO=<ID=Ty,Number=1,Type=String,Description="SNP, INS, DEL or .">\n')
        VCF.write('##INFO=<ID=Rk,Number=1,Type=Float,Description="SNP rank">\n')
        VCF.write('##INFO=<ID=UL,Number=1,Type=Integer,Description="length of the unitig left">\n')
        VCF.write('##INFO=<ID=UR,Number=1,Type=Integer,Description="length of the unitig right">\n')
        VCF.write('##INFO=<ID=CL,Number=1,Type=Integer,Description="length of the contig left">\n')
        VCF.write('##INFO=<ID=CR,Number=1,Type=Integer,Description="length of the contig right">\n')
        VCF.write('##INFO=<ID=Genome,Number=1,Type=String,Description="Allele of the reference;for indel reference is . ">\n')
        VCF.write('##INFO=<ID=Sd,Number=1,Type=Integer,Description="Reverse (-1) or Forward (1) Alignement">\n')        
        VCF.write('##INFO=<ID=XA,Number=0/1,Type=String,Description="Other mapping positions (chromosome_position). Position is negative in case of Reverse alignment. The position designs the starting position of the alignment, not the position of the variant itself.">\n')        
        
        

        ##Creates the columns of the VCF File with all the fields + one field by genotypes/samples/individuals
        if nbGeno==0: # Without genotypes
            VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        else:
            i=0
            VCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            VCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">\n')
            VCF.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled Genotype Likelihoods">\n')
            VCF.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Depth of each allele by sample">\n')
            if boolQuality:
                VCF.write('##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">\n')
            VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
            for i in range(0,int(nbGeno)):
                nomCol="G"+str(i+1)
                VCF.write(str(nomCol))
                if i<nbGeno-1 :
                        VCF.write("\t" )# Adds a \t except if this is the last genotype
                if i==int(nbGeno)-1:
                    VCF.write("\n")
        samfile.close()

#############################################################################################
#############################################################################################


        
        



