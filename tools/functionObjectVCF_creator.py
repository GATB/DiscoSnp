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
#      UnmappingTreatement(variant_object,vcf_field_object,nbGeno,seq1,seq2):"""Fills VCFfile in ghost mode (from a fasta file)"""
#      Counting(fileToRead)
#      AddPosition(position,setPositions,delta)
#      CheckAtDistanceXBestHits(upper_path,lower_path):"""Prediction validation : check if the couple is validated with only one mapping position """
#      PrintVCFHeader(VCF,listName,fileName,boolmyname):    
#############################################################################################
def InitVariant(line1,line2):
        """Initialization of the variant by taking into accoutn its type"""
        #Object Creation
        if "SNP" in line1 and "nb_pol_1" in line1:
                variant_object=SNP(line1,line2)
        elif "SNP" in line1 and "nb_pol_1" not in line1:
                variant_object=SNPSCLOSE(line1,line2)
        elif "INDEL" in line1:
                variant_object=INDEL(line1,line2)
        else :
                print "!!!!Undefined Variant!!!!"
                return 1                
                        
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
        variant_object.GetPolymorphismFromHeader() 
        #Gets the coverage for each path
        variant_object.upper_path.GetCoverage()
        variant_object.lower_path.GetCoverage()
        #Gives the position corrected for each path by taking into account the shift, deletion, insertion ....
        dicoCloseUp=variant_object.upper_path.CheckPosVariantFromRef(vcf_field_object)
        dicoCloseLow=variant_object.lower_path.CheckPosVariantFromRef(vcf_field_object)
        #Checks if the variant has close SNPs
        if int(variant_object.nb_pol)>1:#Close SNPs
                variant_object.GetDicoClose(dicoCloseUp,dicoCloseLow)#fills the dictionnary with all informations for each snps of the path
                table=variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File
        else:#Indel simple SNP
                variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File  
        #Defines the mapping position for the couple
        variant_object.GetMappingPositionCouple()
        #Checks if the couple is validated with only one mapping position       
        vcf_field_object.GetFilterField(CheckAtDistanceXBestHits(variant_object.upper_path,variant_object.lower_path))
        #Defines the genotype for the couple
        variant_object.GetGenotypes(nbGeno,vcf_field_object)
        return table
#############################################################################################
#############################################################################################
def UnmappingTreatement(variant_object,vcf_field_object,nbGeno,seq1,seq2):
        """Treatement of the couple of path without alignment"""
        table = [0] * 10 #Creates a 10 cols array
        seq1=seq1.rstrip('\n')
        seq2=seq2.rstrip('\n')
        #Gets the sequence of each path
        variant_object.upper_path.GetSeq(seq1)
        variant_object.lower_path.GetSeq(seq2)
        #Fills information of the variant object with the informations of the discosnp++ header        
        variant_object.GetPolymorphismFromHeader()
        #Gets the coverage for each path        
        variant_object.upper_path.GetCoverage()
        variant_object.lower_path.GetCoverage()
        dicoCloseUp=variant_object.upper_path.CheckPosVariantFromRef(vcf_field_object)
        dicoCloseLow=variant_object.lower_path.CheckPosVariantFromRef(vcf_field_object)
        #Checks if the variant has close SNPs        
        if int(variant_object.nb_pol)>1:#Close SNPs
                variant_object.GetDicoClose(dicoCloseUp,dicoCloseLow)#fills the dictionnary with all informations for each snps of the path
                table=variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File
        else:#Indel simple SNP
                variant_object.WhichPathIsTheRef(vcf_field_object)#Checks which path will be considered as reference in the VCF File
        #Defines the mapping position for the couple        
        variant_object.GetMappingPositionCouple()
        #Defines the genotype for the couple
        variant_object.GetGenotypes(nbGeno,vcf_field_object)
        return table         
#############################################################################################
#############################################################################################
def Counting(fileName):
        samfile=open(fileName,'r')
        nbGeno=0
        while True:
                line=samfile.readline()
                if not line: break #End of file
                if line.startswith('@'): continue #We do not read headers
                if nbGeno==0:
                        line=line.rstrip('\r').rstrip('\n').split('\t')
                        for i in line:
                                if 'SNP' or 'INDEL' in i:
                                        nomDisco=i.split('|')
                                        for k in nomDisco:
                                                match=re.match(r'^G',k)
                                                if match:
                                                        nbGeno+=1
                                        samfile.close()
                                        return(nbGeno)                        
#############################################################################################
#position : current position to add at the ensemble
#delta : minimum number of difference allowed between two positions
#############################################################################################
def AddPosition(position,setPositions,delta):
        """Add a position to a set of positions only if the difference between the two is greater than delta"""
        if len(setPositions)==0: #Case : it's the first position add to the set
                setPositions.append(position)
                return(setPositions)
        if abs(position-setPositions[0])>delta: #Checks in the postions to test have a difference greatest than delta with the first position of the set
                setPositions.append(position)
                return(setPositions)
        return(setPositions)
#############################################################################################
#############################################################################################
def CheckAtDistanceXBestHits(upper_path,lower_path):
        """Prediction validation : checks if the couple is validated with only one mapping position """
        
        posUp=upper_path.dicoMappingPos
        posLow=lower_path.dicoMappingPos
        # get the best mapping distance for upper path 
        best_up=1024
        if int(upper_path.mappingPosition)==0 and int(lower_path.mappingPosition)==0:#Checks if paths are unmappped
                return "."
        for position,nbMismatch in posUp.items(): 
                if nbMismatch<best_up:
                        best_up=nbMismatch

        # get the best mapping distance for lower path 
        best_low=1024
        for position,nbMismatch in posLow.items(): 
                if nbMismatch<best_low:
                        best_low=nbMismatch
    
        # get the union of the mapping position at the best mapping positions
        position_set = set()
        for position,nbMismatch in posUp.items():
                if nbMismatch == best_up:
                        position_set.add(position)
                if len(position_set) > 1: 
                        return "MULTIPLE"

        for position,nbMismatch in posLow.items():
                if nbMismatch == best_low:
                        position_set.add(position)
                if len(position_set) > 1: 
                        return "MULTIPLE"
    
        if len(position_set) == 1: 
                return "PASS"
    
        return "."

#############################################################################################
#############################################################################################
def PrintVCFHeader(VCF,listName,fileName,boolmyname):
        ###Header of the VCF file 
        today=time.localtime()
        date=str(today.tm_year)+str(today.tm_mon)+str(today.tm_mday)
        VCF.write('##fileformat=VCFv4.1\n')
        VCF.write('##filedate='+str(date)+'\n')
        VCF.write('##source=VCF_creator\n')
        nbGeno=0
        nbSnp=0
        nbGeno = Counting(fileName)
        if boolmyname:
                VCF.write('##BWA_Options='+str(listName[1])+'\n')
                VCF.write('##SAMPLE=file://'+str(listName[0])+".fa"+'\n')
        else:
                VCF.write('##SAMPLE=file://'+str(fileName)+'\n')
        VCF.write('##REF=<ID=REF,Number=1,Type=String,Description="Allele of the path Disco aligned with the least mismatches">\n')
        VCF.write('##FILTER=<ID=MULTIPLE,Description="Mapping type : PASS or MULTIPLE or .">\n')
        VCF.write('##INFO=<ID=Ty,Number=1,Type=String,Description="SNP, INS, DEL or .">\n')
        VCF.write('##INFO=<ID=Rk,Number=1,Type=Float,Description="SNP rank">\n')
        VCF.write('##INFO=<ID=UL,Number=1,Type=Integer,Description="length of the unitig left">\n')
        VCF.write('##INFO=<ID=UR,Number=1,Type=Integer,Description="length of the unitig right">\n')
        VCF.write('##INFO=<ID=CL,Number=1,Type=Integer,Description="length of the contig left">\n')
        VCF.write('##INFO=<ID=CR,Number=1,Type=Integer,Description="length of the contig right">\n')
        VCF.write('##INFO=<ID=Genome,Number=1,Type=String,Description="Allele of the reference;for indel reference is "." ">\n')
        VCF.write('##INFO=<ID=Sd,Number=1,Type=Integer,Description="Reverse (-1) or Forward (1) Alignement">\n')
        VCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        VCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">\n')
        VCF.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled Genotype Likelihoods">\n')
        VCF.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Depth of each allele by sample">\n')


        ##Creates the columns of the VCF File with all the fields + one field by genotypes/samples/individuals
        if nbGeno==0: # Without genotypes
            VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')
        else:
            i=0
            VCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
            for i in range(0,int(nbGeno)):
                nomCol="G"+str(i+1)
                VCF.write(str(nomCol))
                if i<nbGeno-1 :
                        VCF.write("\t" )# Adds a \t except if this is the last genotype
                if i==int(nbGeno)-1:
                    VCF.write("\n")

#############################################################################################
#############################################################################################
