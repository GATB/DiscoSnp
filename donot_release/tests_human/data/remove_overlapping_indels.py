#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

import sys



filein= open("ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes_indivs_HG00096_HG00100.vcf","r")




def check_word(word):
    for l in word:
        if l != 'A' and l != 'C' and l!= 'G' and l!='T' and l!='N':
            return False
    return True

next_valid_pos=0
for line in filein: 
    if not line: break
    if line.startswith('#'):
        print line,
        continue
    
    pos=int(line.split(",")[1].split("'")[1])-1
    ref=line.split(",")[3].split("'")[1]
    alt=line.split(",")[4].split("'")[1]
    
    if not check_word(ref) or not check_word(alt):
        continue
    
    if pos>next_valid_pos: # remove overlapping indels
        if max(len(ref)-len(alt), len(alt)-len(ref))<=10: # remove large indectectable indels
            print line,
            next_valid_pos=pos+max(len(alt),len(ref))
    
    