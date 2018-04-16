#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


#
#input: disco.vcf
#output : new disco.vcf transforming "cluster_4651_size_14_SNP_higher_path_826058	31	826058_2	A	G	" into "cluster_4651_size_14	.	SNP_higher_path_826058_2	A	G"	
#Warning: the mapping position becomes erroneous as we are not able to map on a cluster. This is why ti is silented
#
#usage : python format_VCF_with_cluster_ids.py disco.vcf
#



def output_newvcf(vcf_file):

    filin = open(vcf_file, 'r')

    while True:
        """cluster_4651_size_14_SNP_higher_path_826058	31	826058	A	G	.	.	Ty=SNP;Rk=1;UL=1;UR=2;CL=.;CR=.;Genome=.;Sd=.	GT:DP:PL:AD:HQ	1/1:7:144,25,5:0,7:0,72	1/1:14:284,46,5:0,14
:0,71	1/1:8:164,28,5:0,8:0,73	1/1:59:1184,182,7:0,59:0,70	1/1:144:2884,438,11:0,144:0,72	1/1:30:604,95,6:0,30:0,70	1/1:37:744,116,6:0,37:0,72	./.:1:.,.,.:0,1:0,71	1/1:31:624,98,6:0,31
:0,69	1/1:10:204,34,5:0,10:0,72	1/1:45:904,140,6:0,45:0,72	1/1:31:624,98,6:0,31:0,73	1/1:11:224,37,5:0,11:0,66	1/1:133:2664,405,10:0,133:0,72	1/1:26:524,83,5:0,26:0,68	1/1:
43:864,134,6:0,43:0,72	1/1:27:544,86,5:0,27:0,72	1/1:6:124,22,5:0,6:0,70	1/1:32:644,101,6:0,32:0,73	1/1:55:1104,170,7:0,55:0,71	1/1:20:404,64,5:0,20:0,73	./.:0:.,.,.:0,0:0,0	./.:
0:.,.,.:0,0:0,0	1/1:8:164,28,5:0,8:0,69	1/1:58:1164,179,7:0,58:0,72	1/1:49:984,152,6:0,49:0,72	1/1:3:64,13,4:0,3:0,69	1/1:78:1564,239,8:0,78:0,72	1/1:37:744,116,6:0,37:0,70	1/1:32:644,1
01,6:0,32:"""

        line = filin.readline()
        if not line: break
        if line.startswith("#"): 
            print (line.rstrip())
            continue
        if not line.startswith("cluster"): continue
        splitline = line.rstrip().split()
        
        cluster_and_size_and_path_ids=splitline[0].split("_") # cluster 4651 size 14 SNP higher path 826058
        # generate the cluster id + size only: 
        cluster_id_and_size=cluster_and_size_and_path_ids[0]
        for i in range(1,4):
            cluster_id_and_size+="_"+cluster_and_size_and_path_ids[i]
        
        # generate the path id
        path_id = ""
        for i in range(4,len(cluster_and_size_and_path_ids)-1):
            path_id+=cluster_and_size_and_path_ids[i]+"_"
        path_id+=splitline[2]
            
        
        print (cluster_id_and_size+"\t.\t"+path_id, end='')
        for i in range (3,len(splitline)):
            print ("\t"+splitline[i],end='')
        print()
    filin.close()
    

vcf_file=sys.argv[1]
output_newvcf(vcf_file)

