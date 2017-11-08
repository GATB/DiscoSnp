#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


#
#input: disco.vcf
#output : nouveau disco.vcf avec un SNP par cluster, celui ayant le moins de génotypes manquants
#
#utilisation : python filter_paralogs.py disco.vcf
#



def store_info(vcf_file):

    dict = {}  ## {num_cluster : [num_SNP ayant le moins de géno manquants, nb_missgeno]}
    filin = open(vcf_file, 'r')
    nb_SNP_tot = 0

    while True:
        """cluster_4651_size_14_SNP_higher_path_826058	31	826058	A	G	.	.	Ty=SNP;Rk=1;UL=1;UR=2;CL=.;CR=.;Genome=.;Sd=.	GT:DP:PL:AD:HQ	1/1:7:144,25,5:0,7:0,72	1/1:14:284,46,5:0,14
:0,71	1/1:8:164,28,5:0,8:0,73	1/1:59:1184,182,7:0,59:0,70	1/1:144:2884,438,11:0,144:0,72	1/1:30:604,95,6:0,30:0,70	1/1:37:744,116,6:0,37:0,72	./.:1:.,.,.:0,1:0,71	1/1:31:624,98,6:0,31
:0,69	1/1:10:204,34,5:0,10:0,72	1/1:45:904,140,6:0,45:0,72	1/1:31:624,98,6:0,31:0,73	1/1:11:224,37,5:0,11:0,66	1/1:133:2664,405,10:0,133:0,72	1/1:26:524,83,5:0,26:0,68	1/1:
43:864,134,6:0,43:0,72	1/1:27:544,86,5:0,27:0,72	1/1:6:124,22,5:0,6:0,70	1/1:32:644,101,6:0,32:0,73	1/1:55:1104,170,7:0,55:0,71	1/1:20:404,64,5:0,20:0,73	./.:0:.,.,.:0,0:0,0	./.:
0:.,.,.:0,0:0,0	1/1:8:164,28,5:0,8:0,69	1/1:58:1164,179,7:0,58:0,72	1/1:49:984,152,6:0,49:0,72	1/1:3:64,13,4:0,3:0,69	1/1:78:1564,239,8:0,78:0,72	1/1:37:744,116,6:0,37:0,70	1/1:32:644,1
01,6:0,32:"""

        line = filin.readline()
        if not line: break
        if line.startswith("#"): continue
        if not line.startswith("cluster"): continue

        nb_SNP_tot += 1

        num_cluster = int(line.split("\t")[0].split("_")[1])
        
        if num_cluster == -1: continue
        id_SNP = line.split("\t")[2]

        if num_cluster not in dict:
            dict[num_cluster] = ["x", 290]

        genotypes = [i.split(":")[0] for i in line.split("\t")[9:]]
        nb_missing = 0

        for ind in genotypes:
            #print(ind)
            if ind[0] != "." : continue 
            nb_missing += 1
       
        if nb_missing >= dict[num_cluster][1]: continue
        dict[num_cluster] = [id_SNP, nb_missing]

    filin.close()
    return dict, nb_SNP_tot, len(dict)



def output_newvcf(vcf_file, dict) :
    
    filin = open(vcf_file, 'r')
    new_vcf = open("1SNPlocus_" + str(vcf_file), 'w')

    while True:

        line = filin.readline() 
        if not line: break
        if line.startswith("#"): 
             new_vcf.write(line)        
             continue

        if not line.startswith("cluster"): continue
        cluster = int(line.split("\t")[0].split("cluster_")[1].split("_size")[0]) 
        if cluster == -1 :  continue
 
        id_SNP = line.split("\t")[2]       
        if id_SNP != dict[cluster][0]: continue
       
        new_vcf.write(line)

    filin.close()


vcf_file=sys.argv[1]

dict, nb_SNP_before, nb_SNP_after = store_info(vcf_file)
output_newvcf(vcf_file, dict)

print("nb SNP before : " + str(nb_SNP_before))
print("nb SNP after  : " + str(nb_SNP_after))
