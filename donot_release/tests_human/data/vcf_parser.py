#! /usr/bin/python
# -*- coding: ISO-8859-1 -*-
import sys

""" Cette fonction permet d'extraire un fichier vcf lu sur l'entree standdard les SNPs différenciellements présents chez deux jeux de données (dont les identifiants sont ecrits en durs ci dessous. """


ind1="HG00096"
ind2="HG00100"

 


## Trouver les indices des individus d'intéret:
for line in sys.stdin:
    values = line.split('\t')
    if len(values) > 11 and values[0]=="#CHROM": # header line
        id_ind1=values.index(ind1)
        id_ind2=values.index(ind2)
        print values[0:9], values[id_ind1], values[id_ind2] # print restricted header 
        break; # Continue with content


def homozygotes_diff(snp1, snp2):
    if snp1=="0|0" and snp2=="1|1":
        return True
    if snp1=="1|1" and snp2=="0|0":
        return True
    return False

def heterozygotes_diff(snp1, snp2):
    if snp1=="0|0" and snp2=="0|0":
        return False
    if snp1=="1|1" and snp2=="1|1":
        return False
    return True


for line in sys.stdin:
    #lecture ligne par ligne de stdin
    values = line.split('\t')
    if len(values) > 11:
        ind1_snp=values[id_ind1].split(':')[0]
        ind2_snp=values[id_ind2].split(':')[0]
#        if homozygotes_diff(ind1_snp, ind2_snp): # VERSION JANICE/VINCENT
        if heterozygotes_diff(ind1_snp, ind2_snp): # VERSION PIERRE
            print values[0:9], values[id_ind1], values[id_ind2]
        

        

