#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
From raw facts, remove those for which the overlap is not perfect. 
Badly overlapping variants may occur when discovering bubbles with b1 or b2. 
      |-----A------|
     |              |
-----|               |
      |        |------T------|
       |----C-|               |----------
              |               |
               |------G------|

Bubble1 ---A/C-----T---
Bubble2 ---C---T/G-----

If a read contains ----A-----T----, it will be mapped both on 
* Bubble1A:  ---A-------T--- and on
* Bubble2T:  ---C-------T--- with one mismatch.

While we must not phase those two paths. 

This script removes from raw phased variants those whose overlapping is not perfect. 

Possible improvement todo: accept mutations on positions were no variant exist (sequencing errors).
'''

__author__ = "Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr"


import sys
import K3000_common as kc

    
def filter_and_print_phased_facts(sequences, raw_fact_file_name):
    '''
    Given a set of indexed sequences, the k value and some compacted facts, prints the corresponding fasta seequences
    Enables also to validate the compacted facts. 
    '''
    mfile = open(raw_fact_file_name)
    conserved=0
    removed=0
    nb_paired_facts=0 
    facts_to_print = {} # key = fact (eg "99l_0;-630l_-6;1129l_-22; -1176h_0;316l_-27;-225h_-27;"), value = abundance
    for line in mfile: 
        # 99l_0;-630l_-6;1129l_-22; -1176h_0;316l_-27;-225h_-27; => 1
        # or
        # 99l_0;-630l_-6;1129l_-22; => 1
        if line[0]=="#": continue
        nb_paired_facts+=1
        abundance = int(line.split()[-1])

        line = line.split("=")[0].strip() # 
        # from  "-1000l_0;-136l_-24;-254h_-18;493l_-16;  -577h_0;-977h_-26;1354h_-25;  =>  1"
        # to    "-1000l_0;-136l_-24;-254h_-18;493l_-16;  -577h_0;-977h_-26;1354h_-25;"

        # if the input fact contains a space it means that this is a paired fact, each part is treated

        printed_fact = ""
        for fact_as_ids in line.split():    
            toprint, _, _, _ = kc.line2seq(fact_as_ids, sequences, False)
            if toprint:
                conserved+=1
                if len(printed_fact)>0: printed_fact+=" "
                printed_fact+=fact_as_ids
            else: removed+=1

        # check if the remaining line is composed of a unique fact composed of a unique variant:
        if len(printed_fact.split())==1:
            if len(printed_fact.strip(";").split(";")) < 2:
                conserved-=1
                removed+=1
                printed_fact=""

        if len(printed_fact)>0: 
            if printed_fact not in facts_to_print: facts_to_print[printed_fact]=0
            facts_to_print[printed_fact]+=abundance

    for fact, abundance in facts_to_print.items():
        print(f"{fact} => {abundance}")       

    mfile.close()
    sys.stderr.write(f"Among {conserved+removed} facts (for {nb_paired_facts} paired facts), {conserved} were conserved and {removed} were removed\n")
    sys.stderr.write(f"Removing redundancies due to filtering, this represents a total of {len(facts_to_print)} (possible paired) facts\n")


def main():
    if len(sys.argv) !=4:
        sys.stderr.write(f"Usage: python {sys.argv[0].split('/')[-1]} disco_coherent.fa  disco_uncoherent.fa phased_alleles_read_set_id_i.txt> filtered_phased_alleles_read_set_id_i.txt\n")
        sys.exit(0)
    sequences=kc.index_sequences(sys.argv[1]) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    sequences=kc.index_sequences(sys.argv[2], sequences) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    filter_and_print_phased_facts(sequences, sys.argv[3])

if __name__ == "__main__":
     main()
