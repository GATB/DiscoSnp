#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 


''' ***********************************************
    
    Script to filter and format discoSnp raw output file (.fa) into a vcf format file (.vcf) (ghost vcf)
    Author - Claire Lemaitre
    
    Usage:
    python3 create_filtered_vcf.py -i disco_bubbles_coherent.fa [-o disco_bubbles_coherent.vcf -m 0.95 -r 0.4]
    
    *********************************************** '''


import sys
import getopt
import random
import re #regular expressions
import time
from vcf_formatting_functions import *

def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]," : discoSnp output filtering and formatting in vcf (ghost vcf)")
    print("-----------------------------------------------------------------------------")
    print("usage: ",sys.argv[0]," -i disco_bubbles.fa")
    print("  -r: min rank value filter (default = 0)")
    print("  -m: max missing value filter (default = 1)")
    print("  -o: output vcf file path (default = stdout)")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:r:m:o:", ["help", "in=", "rank=", "miss=", "out="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    # Default parameters
    fasta_file = 0
    min_rank = 0
    max_miss = 1
    k = 31
    out_file = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--in"):
            fasta_file = arg
        elif opt in ("-r", "--rank"):
            min_rank = float(arg)
        elif opt in ("-m", "--miss"):
            max_miss = float(arg)
        elif opt in ("-o", "--out"):
            out_file = arg
        else:
            assert False, "unhandled option"

    if fasta_file == 0:
        print("option -i (--in) is mandatory")
        usage()
        sys.exit(2)
    else:
        today = time.localtime()
        date = str(today.tm_year) + str(today.tm_mon) + str(today.tm_mday)
        source = sys.argv[0]
        

        # First identifying what kind of fasta we have with the first line
        tig_type = 0 # 0 : no extension, 1 unitig, 2 contig (ie. unitig and contig length are output)
        nb_samples = 0
        nb_fixed_fields = 0 #nb fields before C1_X|C2_Y|... depends if unitig and/or contig lengths have been output
        with_cluster = False
        with open(fasta_file, 'r') as filin:
            for line in filin:
                splitted_1 = line.split("|")
                #cluster :
                if re.match(">cluster_",splitted_1[0]):
                    with_cluster = True
                #nb_samples:
                tig_type = len(re.findall("left_\w+_length",line))
                nb_fixed_fields = 4 + 2*tig_type
                nb_samples = (len(splitted_1) - (nb_fixed_fields + 1))/3
                if nb_samples % 1 != 0:
                    print(f"Warning: could not detect the correct nb of samples : {nb_samples}")
                    sys.exit(2)
                nb_samples = int(nb_samples)
                break

        sys.stdout.close = lambda: None  #make stdout unclosable, to use with and handle both `with open(…)` and `sys.stdout` nicely. cf. https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
        # Now going through all lines
        with open(fasta_file, 'r') as filin, (open(out_file,'w') if out_file else sys.stdout) as filout:

            # Writing Comments and Header
            filout.write(vcf_header(source,date,fasta_file,nb_samples))

            nb_kept_variants = 0
            nb_analyzed_variants = 0
            line_count = 0
            keep_variant = False
            for line in filin:
                line_count += 1
                
                if line_count == 1:
                    keep_variant = False
                    nb_analyzed_variants += 1

                    # Header higher path
                    line = line.strip()
                    splitted_1 = line.split("|")
                    #fasta_4lines = splitted_1[0] + "\n"  #simplified headers for fasta_only and src

                    ## FILTERING
                    #filter rank
                    rank = float(splitted_1[-1].split("rank_")[1])
                    if rank < min_rank:
                        continue
                    # filter missing genotype ratio
                    nb_missing = len(re.findall(r"G\d+_\./\.",line))
                    missing_ratio = nb_missing / nb_samples
                    if missing_ratio >= max_miss:
                        continue
                    
                    keep_variant = True  # for fasta_only mode
                    nb_kept_variants += 1
            
                if keep_variant and line_count == 3:
                    #Header lower path
                    line = line.strip()
                    splitted_2 = line.split("|")

                if line_count == 4:
                    line_count = 0
                    if keep_variant:
                        #now format in vcf format
                        line = line.strip()
                        filout.write(format_vcf(splitted_1, splitted_2, nb_samples, rank, line, ".", ".", tig_type))


        #print(f"{nb_lost_variants} variant bubbles filtered out")
        #print(f"{nb_kept_variants} variant bubbles output out of {nb_tot_variants} ({nb_analyzed_variants} analyzed)")
        sys.stderr.write(f"{nb_kept_variants} variant bubbles output out of {nb_analyzed_variants}\n")


if __name__ == "__main__":
    main()
                      


