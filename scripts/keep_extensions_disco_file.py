#!/bin/python
# -*- coding: utf-8 -*-
###################################
# change extensions in uppercase and replace relative positions of SNPs in the header
# usage : keep_extensions_disco_file.py <disco_snps_file>.fa

import sys


def contigOrUnitig(fa_file):
    #return "unitig" or "contig"
    mode = "unitig"

    for line in fa_file:
        headerComp = line.split("\n")[0].split("|")
    
        for elt in headerComp:
            if elt.startswith("left_contig"): 
                mode = "contig"
                break

        #back to the first line
        fa_file.seek(0, 0)
        return mode

def findShift(line, mode):
    headerComp = line.split("\n")[0].split("|")
    left_shift = 0    
    prefix = 'left_'+mode

    for elt in headerComp:
        if not elt.startswith(prefix): continue
        left_shift = int(elt.split("_")[-1])
        break

    return left_shift

def replacePos(position, left_shift):
    # SNP: P_1:30_C/T
    # INDEL: P_1:30_2_11
    new_pos = int(position.split(":")[1].split("_")[0]) + left_shift
    res = position.split(":")[0] + ":" + str(new_pos)
    for i in range (1,len(position.split(":")[1].split("_"))):
        res += "_" + position.split(":")[1].split("_")[i]
    return res

def replaceHeader(elt, left_shift):
    if elt.startswith("P_"):
        elt = ",".join([replacePos(position, left_shift) for position in elt.split(",")])
    return elt




fa_file_name = sys.argv[1]
output_file_name = sys.argv[2]

fa_file = open(fa_file_name,"r")
fa_file_output = open(output_file_name, "w")
left_shift = 0

mode = contigOrUnitig(fa_file)


for line in fa_file:
    if not line.startswith(">"):
        fa_file_output.write(line.upper())
        continue
    
    left_shift = findShift(line, mode)

    if left_shift == 0:
        fa_file_output.write(line)
        continue
    new_header = "|".join([replaceHeader(elt, left_shift) for elt in line.split("|")])

    fa_file_output.write(new_header)

fa_file.close()
fa_file_output.close()




