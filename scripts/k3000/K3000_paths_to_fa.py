#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Creation of a FA file from a compacted fact int file. 
'''


__author__ = "Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr"

import sys
import K3000_common as kc





    
def generate_sequence_paths(sequences, compacted_fact_file_name, int_facts_format):
    '''
    Given a set of indexed sequences, the k value and some compacted facts, prints the corresponding fasta seequences
    Enables also to validate the compacted facts. 
    '''
    mfile = open(compacted_fact_file_name)
    nb_non_writen=0
    for line in mfile: 
        if line[0] == "#": continue
        # * int_facts_format:
        #   38772_0;-21479_1;27388_3;-494_28;-45551_36;-11894_10;-50927_7;-66981_10;29405_22;34837_1;20095_5;
        # * not int_facts_format:
        #  -1000l_0;-136l_-24;-254h_-18;493l_-16;  -577h_0;-977h_-26;1354h_-25;  =>  1
        
        line = line.split("=")[0].strip() # in case a line contains an abundance, it is indicated by "  =>  1"
        # from  "-1000l_0;-136l_-24;-254h_-18;493l_-16;  -577h_0;-977h_-26;1354h_-25;  =>  1"
        # to    "-1000l_0;-136l_-24;-254h_-18;493l_-16;  -577h_0;-977h_-26;1354h_-25;"

        # if the input fact contains a space it means that this is a paired fact, each part is treated
        for fact_as_ids in line.split():
            # print(fact_as_ids)
            toprint, header, bubble_facts_position_start_stops, full_seq = kc.line2seq(fact_as_ids, sequences, int_facts_format)
            if toprint:
                print(f">{header}\t{bubble_facts_position_start_stops}\n{full_seq}")
            else: nb_non_writen+=1
            
    if nb_non_writen>0:
        sys.stderr.write("Warning, "+str(nb_non_writen)+" facts were removed as their sequence concatenation were not coherent or because they contained non coherent predictions\n")
        
    mfile.close()

def is_int_fact(file_name):
    '''
    Determines if a given file is under the format 
    2468_0;-2708_6;1954_-25;1154_-26; (called 'int_facts')
    or
    -577h_0;-977h_-26;1354h_-25;  =>  1
    '''
    with open(file_name) as my_file:
        while True:
            line = my_file.readline()
            if not line: 
                raise RuntimeError (f'input file {file_name} does not contain correctly formated facts')
            if line[0] == "#" : continue # comment
            if  ';' not in line:
                raise RuntimeError (f'input file {file_name} does not contain correctly formated facts. The line {line} should contain at least a \';\'')
            # here we are check the kind of file:
            line = line.split()[0]
            if 'h' in line or 'l' in line:
                return False
            return True



def main():
    '''
    Creation of a FA file from a compacted fact int file. 

    '''
    sequences=kc.index_sequences(sys.argv[1]) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    sequences=kc.index_sequences(sys.argv[2], sequences) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    int_facts_format = is_int_fact(sys.argv[3])
    generate_sequence_paths(sequences, sys.argv[3], int_facts_format)
    



if __name__ == "__main__":
     main()
