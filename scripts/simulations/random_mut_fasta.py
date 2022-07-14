#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import random, choice
from Bio.Alphabet import IUPAC
import sys


#
#  Mute aléatoirement les bases d'une fichier fasta
#  selon la fréquence choisie  
#  genère une copie mutée du fichier fasta réference
#  rend les positions des bases mutées
#
#

def check_word(word):
    for l in word:
        if l != 'A' and l != 'C' and l!= 'G' and l!='T' and l != 'a' and l != 'c' and l!= 'g' and l!='t' : 
            return False
    return True


def random_mut(fasta_file, mutation_freq):
        
        sequences_ref = []
        sequences_mut = []

        for record in SeqIO.parse(fasta_file, "fasta"):
   
            #### recopie record dans un nouveau genome ref pour que le format soit IUPAC.ambiguousDNA comme la copie mutée... (pas bien) ####

            sequence_ref=''

            for base in record.seq:
                sequence_ref = sequence_ref + base

            record_ref=SeqRecord(Seq(sequence_ref, IUPAC.ambiguous_dna), record.id,'','')
            sequences_ref.append(record_ref)
             
            sequence_mut=''
            
            base_id = 0
            ref = ''
            alt = ''
            
            #### genere record muté et print les bases mutées ########

            for base in record.seq:
                base_id += 1
                val = random()
                if val < mutation_freq:
                    if check_word(base):
                        #print(base)
                        ref = base 
                        if base == base.upper():
                            base = choice([x for x in "ACTG" if x != base])
                        else:
                            base = choice([x for x in "actg" if x != base])
                        alt = base
                        print("{}\t{}\t{}\t{}".format(record.id, base_id, ref, alt))

                sequence_mut = sequence_mut + base
                
            record_mut=SeqRecord(Seq(sequence_mut, IUPAC.ambiguous_dna), record.id,'','')
            sequences_mut.append(record_mut)
      
        SeqIO.write(sequences_ref, "genome_ref.fasta", "fasta")
        SeqIO.write(sequences_mut, "genome_mut.fasta", "fasta")



def main(fasta_file, mutation_freq):
    random_mut(fasta_file, mutation_freq)

if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]))
