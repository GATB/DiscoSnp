#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
from Bio import SeqIO
'''
Projet readMap - M1 info Rennes - 2014-2015 - C. Lemaitre C. Belleannee P. Peterlongo -
http://people.rennes.inria.fr/Pierre.Peterlongo/enseignements/projet-m1-bif-2013/
'''

import pythonDict

import sys
import getopt

import random
from time import time, strftime, gmtime

#indexLib=pythonDict

import resource
def using(point=""):
    usage=resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s
        '''%(point,usage[0],usage[1],
             (usage[2]*resource.getpagesize())/1000000.0 )


def number_of_substitutions(read,genome,position):
    ''' 
    returns the number of substitutions while aligning a read on a genome at a given position.
    Note that the read may starts before the genome (position<0) and end after
    '''
    i_read=0
    substitutions=0
    start=position
    stop=position+len(read)
    if position<0:
        substitutions+=-position # adding read positions that are before the beginning  of the genome as substitutions
        start=0
        i_read=-position
        
    if position+len(read)>len(genome):
        substitutions+=position+len(read)-len(genome) # adding read positions that are after the end of the genome as substitutions
        stop=len(genome)
    
    for i_genome in range(start,stop):
        if read[i_read]!=genome[i_genome]: substitutions+=1
        i_read+=1
    
    return substitutions

        
        
def revComp(sequence):
    baseComplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    S=''.join([baseComplement[base] for base in reversed(list(sequence))])
    return S
    
    
def align_a_read(read,genome,index,k,dmax):
    
    
    alignment_nb=0
    tested_positions = set() # avoid to test several times the same positions
    for i in range(0,len(read)-k+1):
        seed=read[i:i+k]
        seed_positions=pythonDict.queryIndex(index,seed)
        if seed_positions==-1:
            continue
        for pos in seed_positions:
            starting_on_genome=pos-i
            if starting_on_genome in tested_positions: continue
            tested_positions.add(starting_on_genome)
            nbsbst = number_of_substitutions(read,genome,starting_on_genome)
            if nbsbst <=dmax:
                alignment_nb+=1
    return alignment_nb

def align_a_read_or_rev_comp(read,genome,index,k,dmax):
    return align_a_read(read,genome,index,k,dmax)+align_a_read(revComp(read),genome,index,k,dmax)
    
    
def align_ref_bubbles(index, ref_sequence, bubble_file, k, dmax, nb_predicted_bubbles):
    pythonDict.queryIndex(index,ref_sequence[0:k])
    nb_ref_bubbles=0
    nb_mapped=0
    sequence_id=0
    all_sequences=[seq_record.seq for seq_record in SeqIO.parse(bubble_file, "fasta")]
    while sequence_id<len(all_sequences):
        sequence_up = str(all_sequences[sequence_id])
        sequence_down = str(all_sequences[sequence_id+1])
        sequence_id+=2
        nb_ref_bubbles+=1
        bubble_mapped=False;
        if align_a_read_or_rev_comp( sequence_up,  ref_sequence, index, k, dmax)>=1 and align_a_read_or_rev_comp(sequence_down,  ref_sequence, index, k, dmax)>=1: bubble_mapped=True
        if bubble_mapped: nb_mapped+=1
    print "nb_ref_bubbles="+str(nb_ref_bubbles)
    print "nb_predicted_bubbles="+str(nb_predicted_bubbles)
    print "nb_mapped="+str(nb_mapped)
    
    print "precision=%.2f "%(100*nb_mapped/float(nb_predicted_bubbles))
    print "recall=%.2f "%(100*nb_mapped/float(nb_ref_bubbles))
    
        
    

def usage():
    '''Usage'''
    print "-----------------------------------------------------------------------------"
    print sys.argv[0]," : mapper of reads on a reference genome"
    print "-----------------------------------------------------------------------------"
    print "usage: ",sys.argv[0]," -r read_file -g ref_genome_file -o output [-k kmer_size] [-d max_dist] [-s]"
    print "  -b: fasta file containing the reference bubbles"
    print "  -p: fasta file containing the predicted bubbles"
    print "  -k: kmer size for the seeds (default=15)"
    print "  -d: max distance (number of substitutions) to report an alignment (default=5)"
    print "  -h: help"
    print "-----------------------------------------------------------------------------"
    sys.exit(2)


def main():
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:p:k:d:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    k = 15
    dmax=5
    bubble_file=0
    predicted_file=0

    both_strands=1
    index="dict"
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-b"):
            bubble_file = arg
        elif opt in ("-k", "--kmer"):
            k = int(arg)
        elif opt in ("-p"):
            predicted_file = arg
        elif opt in ("-d", "--dmax"):
            dmax = int(arg)
        else:
            assert False, "unhandled option"


    if bubble_file==0 or predicted_file==0:
        print "Missing arguments"
        usage()
        sys.exit(2)
    

    ref_sequence=""
    nb_predictions=0
    for seq in SeqIO.parse(predicted_file, "fasta"):
        nb_predictions+=1
        ref_sequence+=str(seq.seq)+"$"
    
    nb_predicted_bubbles=nb_predictions/2
    
    t0=time()
    # GENERATE The indexes
    ref_index=pythonDict.indexGenome(str(ref_sequence),k)
    
    
    
    # pythonDict.queryIndex(ref_sequence[0:k],ref_index)

    #
    # print(using("After indexation"))

    t1=time()    #
    # print "time for indexing = "+ str(t1-t0)

    # basicHashTable.statsIndex(index,5)
    # ALIGN ALL READS
    align_ref_bubbles(ref_index, ref_sequence, bubble_file, k, dmax, nb_predicted_bubbles)
    
    t2=time()    #
    # print "time for mapping = "+ str(t2-t1)
    # print "total time = "+ str(t2-t0)


if __name__ == "__main__":
    main()




        
        
        
