#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#To apply on sorted vcf by position
#Usage : python filterOnBestDP_multiple_variant_at_same_pos.py <vcf_for_igv> > output.vcf
import sys
import gzip
import re
#Filtre on coverage alt
vcf_for_igv=open(sys.argv[1],"r")
pos_prev=0
compt=0
line_suiv=""
line=""
out=False
while True:
        to_output=False
        line_to_print=None
        line_suiv=None
        DP_max=0
        coverage_ref=[]
        coverage_alt=[]
        list_line_same_pos=[]
        to_output=False        
        if compt==0:#Start to read the first line
                line = vcf_for_igv.readline()
                if not line: break
        if line.startswith("#"): #We just print the header
                print (line.rstrip())
                continue
        list_line_same_pos.append(line) #list with all the lines 
        pos=int(line.split("\t")[1])  
        while True:#Get all the vcf line with the same mapping position     
                line_suiv = vcf_for_igv.readline()
                if not line_suiv : 
                        out=True
                        break
                pos_suiv=int(line_suiv.split("\t")[1])
                if pos!=pos_suiv:
                        line=line_suiv 
                        break
                if pos==pos_suiv:
                     list_line_same_pos.append(line_suiv)
        if len(list_line_same_pos)==1:# If there is only one position in the list : we will print the variant
                line_to_print=list_line_same_pos[0]
                to_output=True
        
        else: #We get the variant with the best DP and we will print it
                i=0
                for i in range(len(list_line_same_pos)):
                        listLine=list_line_same_pos[i].rstrip('\r').rstrip('\n').split('\t') 
                        for item in listLine:
                                matchGeno=re.match(r'^[01][/\|]',item)#finds genotypes fields 1/1 or 0/0 etc ...
                                if matchGeno:
                                        listItem=item.split(":")
                                        listItem=item.split(":")
                                        coverage_ref.append(float(listItem[3].split(",")[0])) # adds AD to the list : GT:DP:PL:AD	0/1:205:800,71,2217:138,67	0/1:208:812,72,2249:140,68
                                        coverage_alt.append(float(listItem[3].split(",")[1]))
                                        DP=(int(listItem[1])) # gets DP: GT:DP:PL:AD	0/1:205:800,71,2217:138,67	0/1:208:812,72,2249:140,68
                                        if int(DP)>int(DP_max):
                                                line_to_print=list_line_same_pos[i]
                                                DP_max=DP
        compt+=1
        if to_output:
                print (line_to_print.rstrip())
        elif line_to_print :
                print (line_to_print.rstrip())
        if out:break         
vcf_for_igv.close() 

