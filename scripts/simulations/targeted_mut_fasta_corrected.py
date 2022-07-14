#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from time import time


  
def index_reference_mut(mut_file):

    index = {}

    filin = open(mut_file, 'r')    
 
    while True:
        line = filin.readline() 
        if not line: break
        if line.startswith("#") : continue
        
        chr = line.split("\t")[0].strip()
        pos = int(line.split("\t")[1])
        ref = line.split("\t")[2].strip()
        alt = line.split("\t")[3].strip()
        if len(ref) == len(alt):
            if chr not in index:
                index[chr] = {}
            index[chr][pos] = {}
            index[chr][pos]["ref"] = ref
            index[chr][pos]["alt"] = alt 

    #print(len(index[chr]))   
    filin.close()
    return index


def targeted_mut(fasta_file, index):
    
    genome = SeqIO.parse(fasta_file, "fasta")
    fasta_file_mut = ("{}_mut".format(fasta_file))
    nb_mut = 0
    sequences_mut = []
  
    for record in genome:
        chr = record.id
        mutable_seq = record.seq.tomutable()
        #print(mutable_seq[0])
        if chr in index: 
            for pos in index[chr]:
                if mutable_seq[pos - 1] == index[chr][pos]["ref"]:
                    mutable_seq[pos - 1] = index[chr][pos]["alt"]
                    nb_mut += 1
       
        record_mut=SeqRecord(mutable_seq, record.id,'','')
        sequences_mut.append(record_mut)
    SeqIO.write(sequences_mut, fasta_file_mut, "fasta")

    print("{} bases have been mutated".format(nb_mut))

def main(fasta_file, mutation_file):
    index = index_reference_mut(mutation_file)
    targeted_mut(fasta_file, index)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
