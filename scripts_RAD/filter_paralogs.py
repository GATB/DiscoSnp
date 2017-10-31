#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


#
#input: disco.vcf
#output : new disco.vcf without SNP located in clusters with more than y% SNP heterozygous in more than x% individuals
#
#Usage : python filter_paralogs.py disco.vcf x y
#



def store_info(vcf_file, x):

    dict = {}  ## {num_cluster : [nb SNP avec 0/1 > x%, nb SNP total]}

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
        if line.startswith("#"): continue
        if not line.startswith("cluster"): continue

        num_cluster = int(line.split("\t")[0].split("_")[1])
        
       	if num_cluster not in dict : 
            dict[num_cluster] = [0,0]

        dict[num_cluster][1] += 1        
        genotypes = [i.split(":")[0] for i in line.split("\t")[9:]]
        nb_het = 0
        nb_geno = 0

        for ind in genotypes:
            #print(ind)
            if ind[0] == "." : continue  ##si géno absent
            if ind[0] == "1" or ind[2] == "0" : ##si homo muté ou ref
                nb_geno += 1
                continue

            nb_het += 1    ##si hétérozygote
            nb_geno += 1
        
        if float(nb_het)/float(nb_geno) < x: continue
        dict[num_cluster][0] += 1

    filin.close()
    return dict, len(dict)


def discard_clusters(dict, y):

    clusters_to_keep = []
    for i in dict:
        if float(dict[i][0])/float(dict[i][1]) >= y : 
            continue

        clusters_to_keep.append(i)

    return clusters_to_keep, len(clusters_to_keep)


def output_newvcf(vcf_file, clusters_to_keep) :
    
    filin = open(vcf_file, 'r')
    new_vcf = open("para_" + str(x) + "_" + str(y) + "_" + str(vcf_file), 'w')

    while True:

        line = filin.readline() 
        if not line: break
        if line.startswith("#"): 
             new_vcf.write(line)        
             continue

        if not line.startswith("cluster"): continue
        cluster = int(line.split("\t")[0].split("cluster_")[1].split("_size")[0]) 
        if cluster not in clusters_to_keep: continue
        new_vcf.write(line)

    filin.close()


vcf_file=sys.argv[1]
x = float(sys.argv[2])
y = float(sys.argv[3])

dict, nb_cluster_tot = store_info(vcf_file, x)
clusters_to_keep, nb_clusters_kept = discard_clusters(dict, y)
output_newvcf(vcf_file, clusters_to_keep)

print(str(nb_clusters_kept) + " on " + str(nb_cluster_tot) + " clusters had less than " + str(y*100) + "% of SNP with less than " + str(x*100) + "% heterygous genotypes") 
