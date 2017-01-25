#!/bin/python
# -*- coding: utf-8 -*-
###################################
#Transforms a 1-based vcf into 0-based vcf
# usage : one2zeroBased_vcf.py vcf_file
import sys
from functionsVCF_creator import *
file=open(sys.argv[1],"r")
VCFfile=open("VCFone2zeroBAsed.vcf","w")

for line in file:
    if "#" == line[0]: 
        print (line,)
        continue
    line=line.rstrip('\n')
    #SNP_higher_path_3       117     3       C       G       .       .       Ty=SNP;Rk=1.00000;UL=86;UR=261;CL=.;CR=.;C1=124,0;C2=0,134      GT:DP:PL        0/0:124:10,378,2484     1/1:134:2684,408,10
    splitted_line = line.split("\t")
    
    # decrease the position value if it is a number:
    if str.isdigit(splitted_line[1].strip()):
        splitted_line[1]=int(splitted_line[1])#-1

    # prints the line
    printOneline(splitted_line,VCFfile)
