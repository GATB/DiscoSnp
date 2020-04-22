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
    file_size=kc.file_size(compacted_facts_fa_file)
    step=0
    while(True):
        if step%10000==0:        kc.update_progress(compacted_facts_fa_file.tell()/file_size)
        step+=1
        pos=compacted_facts_fa_file.tell()
        header_fa=compacted_facts_fa_file.readline()
        if not header_fa: break
        header_fa=header_fa.strip().split()[0]  # remove the starting and ending positions from the headers. TODO use them for provinding the overlap length between nodes. 
        compacted_facts_fa_file.readline().strip() # dont care
        header_to_file_position[kc.generate_header(header_fa[1:])]=pos
        # print(header_fa[1:]," pos ",pos)
    kc.update_progress(1)
    compacted_facts_fa_file.close()
    return header_to_file_position

# def kmers_equals_wildchars(kmerA,kmerB):
#     for i in range(len(kmerA)):
#         if kmerA[i]==kmerB[i]:  continue
#         if kmerA[i]=='N':       continue
#         if kmerB[i]=='N':       continue
#         return False
#     return True

def yield_occurring_positions_reverse(k,kmer, seq):
    for i in range(len(seq)-k,-1,-1):
        if kc.hamming_near_perfect(seq[i:i+k],kmer):
            yield i

def overlap_length(seqA, seqB):
    '''
    For two UPPER CASE sequences that overlap at least by at least k characters, return the length of the largest overlap with hamming distance < 5
    1/ find on seqB all occurrence positions of the last kmer of seqA
    2/ check for each position (from the biggest) that the overlap is perfect
    3/ return the length of the biggest overlap. 
    '''
    k=13 
    last_seqA_kmer=seqA[-k:]
    # print("seqA", seqA)
    # print("seqB", seqB)
    # print("last_seqA_kmer", last_seqA_kmer)
    erro_code=-1
    for i in yield_occurring_positions_reverse(k,last_seqA_kmer, seqB):
        if len(seqA[-i-k:]) != len(seqB[:i+k]): # a sequence is included into another one, we do not print those edges
            erro_code=-2
        if i+k > len(seqA):
            erro_code=-2
        if kc.hamming_near_perfect(seqA[-i-k:], seqB[:i+k]):
            return i
    return erro_code
    
# print(overlap_length("TCAACTACTTATTTGTCGTACAAAACTGTCCCGTACATAGGATGATCTTATTCCCGTACCGGATTTCGTACACAATAACAGGAACAATGTCGATATAAAATTTTCTTCAAATGGCTTCAACCCTTACATTATTATGGCAGACGATGTAAACTCTCTAGTCTTCTCAACTCTATTAATAATACATAGTAGTAGCTATTCAGCCATTTTAAAAACGCAATACAACGTTTGTCCCGTAATATT","CTACTTATTTGTCGTACAAAACTGTCCCGTACATAGGATGATCTTATTCCTGTACCGGATTTCGTACACAATAACAGGAACAATGTCGATATAAAATTTTCTTCAAATGGCTTCAACCCTTACATTATTATGGCAGACGACGTAAACTCTCTAGTCTTCTCAACTCTATTAATAATACATAGTAGTAGCTATTCAGCCATTTTAAAAACGCAATACAACGTTTGTCCCGTAATAT"))

def modify_gfa_file(gfa_file_name, compacted_facts_fa_file_name, header_to_file_position):
    print ("H\t#################")
    print ("H\t# GFA of variants")
    print ("H\t#################")
    print ("H\t# Nodes are (compacted) facts with their read mapping coverage. Eg. \"S       1       gtGCAATAAGAATTGTCTTTCTTATAATAATTGTCCAACTTAGgGTCAATTTCTGTACaaacaaCACCATCCAAt     AS:-577h;-977l;1354l;   SP:0_44;11_64;32_75;    BP:0_41;-26_41;-25_41;  FC:i:52 min:17  max:410 mean:180.0 AC:113;17;410\".")
    print ("H\t#  * field AS stands for \"Allele Set\". It reminds the ids of the variants that compose this fact and that were used for reconstructing the genomic sequence")
    print ("H\t#  * field SP stands for \"Sequence Position\". It indicates for each allele of the fact its starting and ending position in the ACGT sequence")
    print ("H\t#  * field BP stands for \"Bubble Position\". For each allele of the fact it indicates:")
    print ("H\t#     - first: the relative position of first nucleotide of the bubble with respect to the position of the last nucleotide of the bubble of the previous allele. This value is equal to zero for the first allele")
    print ("H\t#     - second: the length of the bubble of the allele") 
    print ("H\t#  * field FC is the coverage of the fact, as provided by the total number of reads that phased at least two alleles of the fact")
    print ("H\t#     - first: the relative position of first nucleotide of the bubble with respect to the position of the last nucleotide of the bubble of the previous allele. This value is equal to zero for the first allele")
    print ("H\t#     - second: the length of the bubble of the allele") 
    print ("H\t#  * fields min, max, and mean stand resp. for the min, max and mean of the read coverage of all alleles")
    print ("H\t#  * field AC stands for \"Allele Coverage\". The number of reads that map each allele is given in the same order as the variant ids (eg. \"17;410;113;\" are respectively, the coverages of variants \"-577h;-977l;1354l\")")
   
    print ("H\t# Four types of edges:")
    print ("H\t#   1. Overlap between facts, Overlap length is >0. Eg, \"L	1	-	29384	+	8M\"")
    print ("H\t#       \"S	1	ACGGACGGACCGT	RC:i:24\", and")
    print ("H\t#       \"S	29384	CGGACCGTACGGACGATCCG;	RC:i:43\".")
    print ("H\t#   2. Facts linked by paired end reads.  Eg \"L	10735	+	29384	+	0M	FC:i:5\".")
    print ("H\t#       The coverage indicates the number of pairend read linking the two facts")
    print ("H\t#       These links have a fake overlap of length 0.")
    print ("H\t#   3. Facts linked by unitigs. The unitig finishing a fact overlaps the unitig starting another fact. Eg \"L	19946	+	11433	+	-1M\".")
    print ("H\t#       These links are directed and validate the facts orientation. ")
    print ("H\t#       These links have a fake overlap of length -1.")
    print ("H\t#   4. Facts sharing at least one SNP identifier.")
    print ("H\t#       These links have an overlap of length -2.")
    
    gfa_file=open(gfa_file_name)
    compacted_facts_fa_file=open(compacted_facts_fa_file_name)
    node_id_to_sequence={}

    file_size=kc.file_size(gfa_file)
    step=0
    while(True):
        if step%10000==0:        kc.update_progress(gfa_file.tell()/file_size)
        step+=1
        gfa_line=gfa_file.readline()
        if not gfa_line: break
        gfa_line.strip()
        if gfa_line[0]=='H': continue       #Header was changed
        if gfa_line[0]=='S':                #Deal with sequences
            #FROM: 
            #S       1       -577h;-977l;1354l;      SP:0_44;11_64;32_75;    BP:0_41;-26_41;-25_41;  FC:i:52 min:17  max:410 mean:180.0      AC:410;113;17;
            #TO
            #\"S       1       gtGCAATAAGAATTGTCTTTCTTATAATAATTGTCCAACTTAGgGTCAATTTCTGTACaaacaaCACCATCCAAt    SP:0_44;11_64;32_75;    BP:0_41;-26_41;-25_41;  FC:i:52   AS:-577h;-977l;1354l;      min:17    max:410 mean:180.0      AC:17;410;113;

            gfa_line=gfa_line.split()
            assert gfa_line[2] in header_to_file_position, gfa_line[2]+" is not in header_to_file_position"
            compacted_facts_fa_file.seek(header_to_file_position[gfa_line[2]])
            compacted_facts_fa_file.readline().strip() # dont care
            sequence_fa=compacted_facts_fa_file.readline().strip()
            # assert gfa_line[2] == allele_header,gfa_line[2]+" is not "+allele_header+" "+header_fa[1:]
            node_id_to_sequence[gfa_line[1]]=sequence_fa                                #TODO: optimize this to avoid memory usage. One may store the position of the node in the file and retreive the sequence latter
            
            print(gfa_line[0]+"\t"+gfa_line[1]+"\t"+sequence_fa+"\tAS:"+gfa_line[2]+"\t"+gfa_line[3]+"\t"+gfa_line[4]+"\t"+gfa_line[5]+"\t"+gfa_line[6]+"\t"+gfa_line[7]+"\t"+gfa_line[8]+"\t"+gfa_line[9]) 
            continue
        
        if gfa_line[0]=='L':
            split_gfa_line=gfa_line.split()
            if split_gfa_line[1] == split_gfa_line[3]: #no not print self loops
                continue
            # print (split_gfa_line)
            if split_gfa_line[5]=="0M" or split_gfa_line[5]=="-1M" or split_gfa_line[5]=="-2M": # non overlapping edges, we simply write them
                print(gfa_line.strip())
                continue
            # if we are here, this is a true overlapping edge: L	3	+	255	-	2M
            # we need to retreive the sequences of the two nodes
            # print(split_gfa_line)
            seqA = node_id_to_sequence[split_gfa_line[1]].upper()
            seqB = node_id_to_sequence[split_gfa_line[3]].upper()
            if split_gfa_line[2]=='-': seqA=kc.get_reverse_complement(seqA)
            if split_gfa_line[4]=='-': seqB=kc.get_reverse_complement(seqB)
            OL = overlap_length(seqA,seqB)
            # assert OL!=-1,seqA+" "+seqB
            if OL>-1:
                print (split_gfa_line[0]+"\t"+split_gfa_line[1]+"\t"+split_gfa_line[2]+"\t"+split_gfa_line[3]+"\t"+split_gfa_line[4]+"\t"+str(OL)+"M")
            continue
            
            
        print(gfa_line.strip()) # shold not happend

    kc.update_progress(1)
    gfa_file.close()
    compacted_facts_fa_file.close()


def main():
    '''
    Produces a gfa file replacing the node content from int ids of alleles to their sequence
    '''
    if len(sys.argv) !=3:
        sys.stderr.write("Usage: python K3000_node_ids_to_node_sequences.py graph_plus.gfa compacted_facts.fa > graph_final.gfa\n")
        sys.exit(0)
    sys.stderr.write("Indexing sequence positions\n")
    header_to_file_position = index_sequences(sys.argv[2])
    sys.stderr.write("Writing the updated gfa file\n")
    modify_gfa_file(sys.argv[1],sys.argv[2], header_to_file_position)
    



if __name__ == "__main__":
     main()
