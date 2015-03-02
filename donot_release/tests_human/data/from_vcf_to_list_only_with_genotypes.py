#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
            

print "read genomes"
filein= open("ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes_indivs_HG00096_HG00100_no_overlapping_indels.vcf","r")






f = file('ref_human_genotype', 'w')

print "mutate genomes"

prev_pos=0
for line in filein: 
    if not line: break
    if line.startswith('#'):
        continue 
    
    # 
    pos=int(line.split(",")[1].split("'")[1])-1
    ref=line.split(",")[3].split("'")[1]
    alt=line.split(",")[4].split("'")[1]
    mod_96=line.split("]")[1].split(" ")[1].split(":")[0]
    mod_100=line.split("]")[1].split(" ")[2].split(":")[0]

    

    var="SNP"
    if len(ref) != len(alt):
        var="INDEL"
    
    # print the position in the reference list
    print >> f, var, str(pos),"1", ref, alt, mod_96, mod_100
    
   
f.close()


