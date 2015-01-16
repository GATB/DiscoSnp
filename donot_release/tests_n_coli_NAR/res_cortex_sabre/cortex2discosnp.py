#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys

def cortex2disco(cortexFile):
    filin = open(cortexFile, 'r')
    snp_id=0
    indel_id=0
    while 1:
        line=filin.readline() #>var_18_5p_flank length:1031 INFO:KMER=31
        if not line: break
        k=int(line.split("=")[-1])
        id_pol=int(line.split("_")[1])

        
        start=filin.readline()[-k:-1]# get the last k characters from the 5p sequence
        filin.readline() # >var_17_branch_1 length:32 kmer:31 we don't care
        var1=filin.readline()
        filin.readline() # >var_17_branch_2 length:32 kmer:31 we don't care
        var2=filin.readline()
        filin.readline() # >var_17_3p_flank length:8 kmer:31 we don't care
        filin.readline() # 3p sequence, we don't care (the 2 var1 and var2 already contain the 3p kmer)
        
        if len(var1) == len(var2): # SNP
            print ">SNP_higher_path_"+str(snp_id)
            print start+var1,
            print ">SNP_lower_path_"+str(snp_id)
            print start+var2,
            snp_id+=1
        else:
            print ">INDEL_higher_path_"+str(indel_id)
            print start+var1,
            print ">INDEL_lower_path_"+str(indel_id)
            print start+var2,
            indel_id+=1
            
            
            
cortex2disco(sys.argv[1])
        
       