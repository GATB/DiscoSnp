#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

import sys



filein= open("ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes_indivs_HG00096_HG00100_no_overlapping_indels.vcf","r")
refgenome=str(list(SeqIO.parse("humch1.fasta", "fasta"))[0].seq)
print "read genomes"
genomes=[""]*4
genomes[0]=""
genomes[1]=""
genomes[2]=""
genomes[3]=""

f = file('ref_human', 'w')

print "mutate genomes"

prev_pos=0
for line in filein: 
    if not line: break
    if line.startswith('#'):
        continue 


    pos=int(line.split(",")[1].split("'")[1])-1


    ref=line.split(",")[3].split("'")[1]
    alt=line.split(",")[4].split("'")[1]
    mod_96=line.split("]")[1].split(" ")[1].split(":")[0]
    do_genomes=[False]*4
    if mod_96.split("|")[0]=="1":
        do_genomes[0]=True # 96_1
    if mod_96.split("|")[1]=="1":
        do_genomes[1]=True # 96_2
    
    mod_100=line.split("]")[1].split(" ")[2].split(":")[0]

    if mod_100.split("|")[0]=="1":
        do_genomes[2]=True # 100_1
    if mod_100.split("|")[1]=="1":
        do_genomes[3]=True # 100_2

    
    var="SNP"
    if len(ref) != len(alt):
        var="INDEL"
        # print pos,do_genomes,ref,alt
    if refgenome[pos:pos+len(ref)] != ref:
        print "WARNING, pos",int(pos)," ref declared is",ref,"while in the line it is",ref[pos:pos+len(ref)]
    
    print >> f, var, str(pos),"1"
    for genome_id in range(4):
        genomes[genome_id]+=refgenome[prev_pos:pos]
        if do_genomes[genome_id]:
            genomes[genome_id]+=alt
        else:
            genomes[genome_id]+=ref
    prev_pos=pos+len(ref)


f.close()
f = file('humch1_00096_1.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00096_1 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
print genomes[0]
f.close()

f = file('humch1_00096_2.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00096_2 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
print genomes[1]
f.close()

f = file('humch1_00100_1.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00100_1 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
print genomes[2]
f.close()


f = file('humch1_00100_2.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00100_2 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
print genomes[3]
f.close()


