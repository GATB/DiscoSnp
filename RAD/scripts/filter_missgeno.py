#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys


#
#  elimine les predictions ayant plus de x (fraction de 1) données manquantes dans un fichier disco.fa
#  python filter_missgeno.py disco.fa x
#  output : Xmissing_disco.fa
#

def filter_MAF(fasta_file, max_missing):

    simple_name = str(fasta_file)

    if "/" in str(fasta_file):
        simple_name = str(fasta_file).split("/")[-1]
     
    output = open(str(max_missing) + "missing_" + simple_name, 'w')

    filin = open(fasta_file, 'r')    
    keep_seq = 0

    while True :

        line = filin.readline()
        if not line: break
        if not line.startswith(">"):
            if keep_seq != 1: continue
            output.write(line)
            keep_seq = 0
            continue

        header = [x for x in line.split("|")]
        cpt_geno = 0
        cpt_missing = 0    
 
        for i in header:
            if not i.startswith("G"): continue
            geno = i.split(":")[0].split("_")[1]
            if geno.startswith("."):
                 cpt_geno += 1
                 cpt_missing += 1
                 continue

            cpt_geno += 1

        missing_percent = float(cpt_missing)/float(cpt_geno)
        if missing_percent > max_missing: continue

        output.write(line)
        keep_seq = 1
         
    output.close()
    filin.close()

def main(fasta_file, max_missing):
    filter_MAF(fasta_file, max_missing)

if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]))
