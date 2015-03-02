#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

import sys

def check_word(word):
    for l in word:
        if l != 'A' and l != 'C' and l!= 'G' and l!='T' and l!='N':
            return False
    return True

class genome:
    def __init__(self, size):
        self.size=int(size)
        self.table=["\0"]*int(size)
        self.pos=0
        
    
    def add_word(self,word):
        check_word(word)
        for l in word:
            self.table[self.pos]=l
            self.pos+=1
            if self.pos==self.size:
                return False
        return True
    
    def print_genome(self):
        for i in xrange(self.pos):
            sys.stdout.write(self.table[i])
            


print "read genomes"
filein= open("ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes_indivs_HG00096_HG00100_no_overlapping_indels.vcf","r")
refgenome=str(list(SeqIO.parse("humch1.fasta", "fasta"))[0].seq)

genomes=[""]*4
genomes[0]=genome(len(refgenome)*1.1)
genomes[1]=genome(len(refgenome)*1.1)
genomes[2]=genome(len(refgenome)*1.1)
genomes[3]=genome(len(refgenome)*1.1)



f = file('ref_human', 'w')

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

    if not check_word(ref) or not check_word(alt):
        print line, ref, alt
    var="SNP"
    if len(ref) != len(alt):
        var="INDEL"
        # print pos,do_genomes,ref,alt
    if refgenome[pos:pos+len(ref)] != ref:
        print "WARNING, pos",int(pos)," ref declared is",ref,"while in the line it is",ref[pos:pos+len(ref)]
    
    # print the position in the reference list
    print >> f, var, str(pos),"1", ref, alt
    
    # fill each of the 4 genomes up to the end of the polymorphism
    for genome_id in xrange(4):
        
        # copy what is missing before the variant position:
        genomes[genome_id].add_word(refgenome[prev_pos:pos])
        
        # copy what is in the alt if the genome is polymorphic at this position
        if do_genomes[genome_id]:
            genomes[genome_id].add_word(alt)
        # copy what is in the ref if the genome is not polymorphic at this position
        else:
            genomes[genome_id].add_word(ref)

    prev_pos=pos+len(ref) # next position on the genome
    
    ############### DEBUG ####################
    # print
   #  print str(pos), ref+" --> "+alt
   #  print do_genomes
   #  print "reference:"
   #  print refgenome[pos-10:pos+10]
    ############## FIN DEBUG #################

#ENDING THE GENOMES:
for genome_id in range(4):
    # copy what is missing before the end of the genome
    genomes[genome_id].add_word(refgenome[prev_pos:len(refgenome)])
f.close()

print "output muted genomes"

f = file('humch1_00096_1.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00096_1 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
genomes[0].print_genome()
f.close()

f = file('humch1_00096_2.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00096_2 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
genomes[1].print_genome()
f.close()

f = file('humch1_00100_1.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00100_1 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
genomes[2].print_genome()
f.close()


f = file('humch1_00100_2.fasta', 'w')
sys.stdout = f
print ">gi|224384768|gb|CM000663.1|HG00100_2 Homo sapiens chromosome 1, GRCh37 primary reference assembly"
genomes[3].print_genome()
f.close()


