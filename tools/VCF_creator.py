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
from functionObjectVCF_creator import *
from ClassVCF_creator import *



#Default value
nbMismatchBWA=3 #number of mismatch allowed for the alignment in BWA
listName=[] #List from the file name used to create the VCF
nbSnp=0 #number of SNP in the input file
nbGeno=0 #number of genotype for every path
filtered_sam=False
filtered_sam_file=""
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
    -f --output_filtered_SAM : if provided, a SAM file in which uncorrectly mapped prediction (corresponding to filter '.' in the provided VCF) are removed is output in this file.
    
    """
    print usage
###OPTIONS 
try:
    opts, args = getopt.getopt(sys.argv[1:],"h:s:n:o:f:",["help","disco_file","mismatch","output=","output_filtered_SAM"])
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
        if os.path.isfile(arg):#checks if the file exists
               fileName = arg
               listNameFile=fileName.split(".")
               if "BWA_OPT" in listNameFile[0]: #When the samfile is created by run_VCF_creator.sh ; it adds BWA_OPT to separate the name of the file and the BWA options
                       boolmyname=True #Boolean to know if the file name contains the option of BWA
                       listName=listNameFile[0].split("BWA_OPT") #Parsing of the file name listName[0] corresponds to the name of the discofile ; listName[1] corresponds to the bwa options
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
                VCFFile = open(arg,'w')
        else :
                print "!! No output !!"
                sys.exit(2)      
    elif opt in ("-f","--output_filtered_SAM"):
        if arg!=None:
            filtered_sam = True
            filtered_sam_file = open(arg,'w')
        else:
            print "!! No filtered sam output !!"
            sys.exit(2)      
    else:
        print("Unkwnown option {} ".format(opt))
        usage()
        sys.exit(2)

#Creates the VCF file and prints the header into it 
PrintVCFHeader(VCFFile,listName,fileName,boolmyname)
nbGeno=Counting(fileName)
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

#Start to read the file two lines by two lines
if ".sam" in fileName: #Checks if it's a samfile
        while True:
                line1=samfile.readline()
                if not line1: break #End of file
                if line1.startswith('@'): 
                    if filtered_sam:
                        filtered_sam_file.write(line1)
                    continue #We do not read headers
                line2=samfile.readline() #Read couple of lines
                #Initializes variant object with the samline
                variant_object, vcf_field_object=InitVariant(line1,line2) #Fills the object with the line of the samfile        
                if variant_object.CheckCoupleVariantID()==1: #Checks whether the two lines are from the same path
                        sys.exit(1)
                #Checks the mapping on reference and determines the shift with the reference, which path is the reference ...
                table=MappingTreatement(variant_object,vcf_field_object,nbGeno)
                #Added by Pierre Peterlongo, Sept 2015. Outputs a new SAM file corresponding to mapped sequences (PASS or MULTIPLE)
                if(filtered_sam and vcf_field_object.filterField!="."):
                    filtered_sam_file.write(line1)
                    filtered_sam_file.write(line2)
                
                #Fills the VCF file
                variant_object.FillVCF(VCFFile,nbGeno,table,vcf_field_object)
                  

            

elif ".fa" in fileName: #Treatement of the fasta file (no mapping information)
        while True:
                line1=samfile.readline()
                if not line1: break #End of file
                seq1=samfile.readline() #Reads the seq associate to the SNP
                line2=samfile.readline() #Reads a couple of line
                seq2=samfile.readline()
                #Initializes variant object with the samline
                variant_object, vcf_field_object=InitVariant(line1,line2)
                if variant_object.CheckCoupleVariantID()==1:
                        sys.exit(1) 
                table=UnmappingTreatement(variant_object,vcf_field_object,nbGeno,seq1,seq2)
                variant_object.FillVCF(VCFFile,nbGeno,table,vcf_field_object)


VCFFile.close()
samfile.close()
if filtered_sam:
    filtered_sam_file.close()





