#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Creation of a FA file from a compacted fact int file. 
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
import K3000_common as kc

def index_sequences(compacted_facts_fa_file_name):
    '''
    Stores for each sequence header its position in the fa file. This is used latter to retrieve the corresponding sequences
    '''
    header_to_file_position = {}
    compacted_facts_fa_file=open(compacted_facts_fa_file_name)
    while(True):
        pos=compacted_facts_fa_file.tell()
        header_fa=compacted_facts_fa_file.readline()
        if not header_fa: break
        header_fa=header_fa.strip().split()[0]  # remove the starting and ending positions from the headers. TODO use them for provinding the overlap length between nodes. 
        sequence_fa=compacted_facts_fa_file.readline().strip()
        header_to_file_position[kc.generate_header(header_fa[1:])]=pos
        # print(header_fa[1:]," pos ",pos)
    compacted_facts_fa_file.close()
    return header_to_file_position


def modify_gfa_file(gfa_file_name, compacted_facts_fa_file_name, header_to_file_position):
    print ("H\t#################")
    print ("H\t# GFA of variants")
    print ("H\t#################")
    print ("H\t# Nodes are (compacted) facts with their read mapping coverage. Eg. \"S	2	ACGGACGGACCGT;	RC:i:144\".")
    print ("H\t# Three types of edges:")
    print ("H\t#   1. Overlap between facts. These links have a fake overlap length of 1. Eg, \"L	1	-	29384	+	1M\", with:")
    print ("H\t#       \"S	1	ACGGACGGACCGT	RC:i:24\", and")
    print ("H\t#       \"S	29384	CGGACCGTACGGACGATCCG;	RC:i:43\".")
    print ("H\t#   2. Facts linked by paired end reads.  Eg \"L	10735	+	29384	+	0M	FC:i:5\".")
    print ("H\t#       These links are non directed and do no validate the facts orientation. The coverage indicates the number of pairend read linking the two facts")
    print ("H\t#       These links have a fake overlap of length 0.")
    print ("H\t#   3. Facts linked by unitigs. The unitig finishing a fact overlaps the unitig starting another fact. Eg \"L	19946	+	11433	+	-1M\".")
    print ("H\t#       These links are directed and validate the facts orientation. ")
    print ("H\t#       These links have a fake overlap of length -1.")
    
    gfa_file=open(gfa_file_name)
    compacted_facts_fa_file=open(compacted_facts_fa_file_name)
    while(True):
        gfa_line=gfa_file.readline()
        if not gfa_line: break
        gfa_line.strip()
        if gfa_line[0]=='H': continue
        if gfa_line[0]!='S':
            print(gfa_line.strip())
            continue
        #S	0	102h;100l;168h;	RC:i:11
        gfa_line=gfa_line.split()
        assert gfa_line[2] in header_to_file_position, gfa_line[2]+" is not in header_to_file_position"
        compacted_facts_fa_file.seek(header_to_file_position[gfa_line[2]])
        header_fa=compacted_facts_fa_file.readline().strip()
        sequence_fa=compacted_facts_fa_file.readline().strip()
        # assert gfa_line[2] == allele_header,gfa_line[2]+" is not "+allele_header+" "+header_fa[1:]
        print(gfa_line[0]+"\t"+gfa_line[1]+"\t"+sequence_fa+"\t"+gfa_line[3])
        
    gfa_file.close()
    compacted_facts_fa_file.close()


def main():
    '''
    Produces a gfa file replacing the node content from int ids of alleles to their sequence
    '''
    if len(sys.argv) !=3:
        sys.stderr.write("Usage: python K3000_node_ids_to_node_sequences.py graph_plus.gfa compacted_facts.fa > graph_final.gfa\n")
        sys.exit(0)
    header_to_file_position = index_sequences(sys.argv[2])
    modify_gfa_file(sys.argv[1],sys.argv[2], header_to_file_position)
    



if __name__ == "__main__":
     main()
