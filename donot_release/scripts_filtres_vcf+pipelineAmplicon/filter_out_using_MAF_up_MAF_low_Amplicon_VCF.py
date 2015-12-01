#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import gzip
import re
#Version UP and Low MAF
if len(sys.argv)<3:
    print "This tool filters out discoSnp prediction having a minor allele frequency lower than a provided threshold for ALL datasets."
    print "python filter_out_using_MAF.py \".fa from discoSnp\" \"MAF threshold\""
    sys.exit()


if "gz" in sys.argv[1]:
    coherent_file=gzip.open(sys.argv[1],"r")
else: 
    coherent_file=open(sys.argv[1],"r")

try:
        maf_threshold_low = float(sys.argv[2])
except IndexError :
        maf_threshold_low = 0       
try:
        maf_threshold_up = float(sys.argv[3])
except IndexError :
       maf_threshold_up = 0  
first_coverage_field_id=0
while True:
    coverage_ref=[]
    coverage_alt= [] 

    if ".vcf" in str(sys.argv[1]):
        line = coherent_file.readline() 
        if not line: break
        if line.startswith("#") : 
                print line.rstrip()
                continue
        listLine=line.rstrip('\r').rstrip('\n').split('\t')
        for item in listLine:
                matchGeno=re.match(r'^[01][/\|]',item)#finds genotypes fields 1/1 or 0/0 etc ...
                if matchGeno:
                        listItem=item.split(":")
                        coverage_ref.append(float(listItem[3].split(",")[0])) # adds AD to the list : GT:DP:PL:AD	0/1:205:800,71,2217:138,67	0/1:208:812,72,2249:140,68
                        coverage_alt.append(float(listItem[3].split(",")[1]))
    #print comment1,
    #print coverage_high
    #print coverage_low
    to_output=False
    if sum(coverage_alt)==0: continue
    if len(coverage_alt)==1:
        if coverage_alt[0]/(coverage_alt[0]+coverage_ref[0])>=maf_threshold_low and (coverage_alt[0]+coverage_ref[0])>=150:
                to_output=True
    else:
            try:
                    if float(coverage_alt[0])/(float(coverage_alt[0])+float(coverage_ref[0]))>=float(maf_threshold_low) and (float(coverage_alt[0])+float(coverage_ref[0]))>=150:
                        #if coverage_alt[1]/(coverage_alt[1]+coverage_ref[1])>=maf_threshold_low and (coverage_alt[1]+coverage_ref[1])>=150:#AND
                        to_output=True
                    elif float(coverage_alt[1])/(float(coverage_alt[1])+float(coverage_ref[1]))>=float(maf_threshold_low) and (float(coverage_alt[1])+float(coverage_ref[1]))>=150:
                        to_output=True
            except ZeroDivisionError:
                to_output=False
            
            
    
    if to_output:
        if ".fa" in str(sys.argv[1]): 
                print comment1,path1,comment2,path2,
        elif ".vcf" in str(sys.argv[1]):
                print line.rstrip()
                
coherent_file.close() 
