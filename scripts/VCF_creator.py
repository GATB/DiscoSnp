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


                

#Help
def usage():
        usage= """
        ################################
                Run VCF_creator
        ################################
        
        -h --help : print this message
        -s --sam_file : <file>.sam of the alignment

        -o --output : vcf file 
        -f --output_filtered_SAM : if provided, a SAM file in which uncorrectly mapped prediction (corresponding to filter '.' in the provided VCF) are removed is output in this file.
        
        """
        print(usage)
        
def main():
        #Default value
        listName=[] #List from the file name used to create the VCF
        nbGeno=0 #number of genotype for every path
        filtered_sam=False
        filtered_sam_file=""
        VCFFile = None
        ###OPTIONS 
        try:
                opts, args = getopt.getopt(sys.argv[1:],"s:o:f:h",["help","sam_file=","output=","output_filtered_SAM="])
                if not opts:
                        usage()
                        sys.exit(2)
        except getopt.GetoptError as e:
                sys.stderr.write(f"--- {e}\n")
                usage()
                sys.exit(2)
        for opt, arg in opts : 
                if opt in ("-h", "--help"):
                        usage()
                        sys.exit(2)
                elif opt in ("-s","--sam_file"):
                        boolmyname=False
                        fileName=arg
                        if os.path.isfile(fileName):#checks if the file exists
                                listNameFile=fileName.split(".")
                                if "BWA_OPT" in listNameFile[0]: #When the stream_file is created by run_VCF_creator.sh ; it adds BWA_OPT to separate the name of the file and the BWA options
                                        boolmyname=True #Boolean to know if the file name contains the option of BWA
                                        listName=listNameFile[0].split("BWA_OPT") #Parsing of the file name listName[0] corresponds to the name of the discofile ; listName[1] corresponds to the bwa options
                                        listName[1]=listName[1].replace("_", " ")
                                stream_file=open(fileName,'r')
                        else : 
                                sys.stderr.write(f"File \"{arg}\" does not exist")
                                sys.exit(2)  
                elif opt in ("-o","--output"):
                        if not arg:
                                sys.stderr.write("-o needs a value\n")
                                sys.exit(2)
                        
                        VCFFile = open(arg,'w')
                elif opt in ("-f","--output_filtered_SAM"):
                        if arg!=None:
                                filtered_sam = True
                                filtered_sam_file = open(arg,'w')
                        else:
                                sys.stderr.write("!! No filtered sam output !!\n")
                                sys.exit(2)          
                else:
                        sys.stderr.write(f"Unknown option {format(opt)}\n")
                        usage()
                        sys.exit(2)

        if not VCFFile:
                VCFFile = sys.stdout
                
        #Creates the VCF file and prints the header into it 
        PrintVCFHeader(VCFFile,listName,fileName,boolmyname)
        nbGeno=CounterGenotype(fileName)
        #---------------------------------------------------------------------------------------------------------------------------
        #---------------------------------------------------------------------------------------------------------------------------
        #Generates the field index (first occurrence of "contig", and so on)
        fieldIndex=GetIndex(fileName)


        #---------------------------------------------------------------------------------------------------------------------------
        #---------------------------------------------------------------------------------------------------------------------------
        #Start to read the file two lines by two lines
        if ".sam" in fileName: #Checks if it's a stream_file
                while True:
                        line1=stream_file.readline() 
                        if not line1: break #End of file  
                        line1 = line1.strip()
                        if line1.startswith('@'):
                                if filtered_sam:
                                        filtered_sam_file.write(line1) 
                                continue #We do not read headers                                
                        
                        while True:
                                listline1=line1.split("\t") 
                                if len(listline1) > 1 and int(listline1[1]) & 2048 :#checks if it's not a secondary alignment => means splitted aligned sequence 
                                        line1=stream_file.readline()
                                else: break
                        line2=stream_file.readline().strip() #Read couple of lines
                        while True:
                                listline2=line2.split("\t") 
                                if len(listline2) > 1 and int(listline2[1]) & 2048 :#checks if it's not a secondary alignment => means splitted aligned sequence 
                                        line2=stream_file.readline()
                                        if not line2: break # happens when second part of the previous pair (last one) add an unread secondary alignment
                                else: break           
                        #Initializes variant object with the samline
                        #if InitVariant(line1,line2)[0]==1:
                        #                continue
                        if not line2: break
                        variant_object, vcf_field_object=InitVariant(line1,line2,fileName, fieldIndex) #Fills the object with the line of the stream_file                
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
                        line1=stream_file.readline()
                        if not line1: break #End of file
                        seq1=stream_file.readline() #Reads the seq associate to the variant
                        line2=stream_file.readline() #Reads a couple of line
                        seq2=stream_file.readline()
                        #Initializes variant object with the samline
                        variant_object, vcf_field_object=InitVariant(line1,line2,fileName, fieldIndex)
                        table=UnmappedTreatement(variant_object,vcf_field_object,nbGeno,seq1,seq2)
                        variant_object.FillVCF(VCFFile,nbGeno,table,vcf_field_object)


        VCFFile.close()
        stream_file.close()
        if filtered_sam:
                filtered_sam_file.close()





if __name__ == '__main__':
        main()
        sys.stderr.write(f"VCF file had been generated. Beware, this is a zero-based VCF file.\n")
        sys.stderr.write(f"Please use 'zero2one.py' script if you wish to obtain a one-based VCF.\n")
        
