#!/bin/python
# -*- coding: utf-8 -*-
###################################
# change extensions in uppercase and replace relative positions of SNPs in the header
# usage : keep_extensions_disco_file.py <disco_snps_file>.fa

import os
import sys
import getopt
import time


fichier = sys.argv[1]
output = sys.argv[2]

snp_file = open(fichier,"r")
snp_discobis = open(output, "w")
left_shift = 0

def findShift(line):
    headerComp = line.split("\n")[0].split("|")
    for elt in headerComp:
        if not elt.startswith("left"): continue
        left_shift = int(elt.split("_")[-1])
        break
    return left_shift

def replacePos(snp, left_shift):
    new_pos = int(snp.split(":")[1].split("_")[0]) + left_shift
    return snp.split(":")[0] + ":" + str(new_pos) + "_" + snp.split(":")[1].split("_")[1]

def replaceHeader(elt, left_shift):
    if elt.startswith("P_"):
        elt = ",".join([replacePos(snp, left_shift) for snp in elt.split(",")])
    return elt


for line in snp_file:
    if not line.startswith(">"):
        snp_discobis.write(line.upper())
        continue
        
    left_shift = findShift(line)

    if left_shift == 0:
        snp_discobis.write(">" + header)
        continue

    new_header = "|".join([replaceHeader(elt, left_shift) for elt in line.split("|")])

    snp_discobis.write(">" + new_header)

snp_file.close()
snp_discobis.close()
