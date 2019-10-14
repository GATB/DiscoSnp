#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 


''' ***********************************************
    
    Script to filter and format discoSnp raw output file (.fa) into a vcf format file (.vcf)
    Author - Claire Lemaitre
    
    Usage:
    python3 fasta_and_cluster_to_filtered_vcf.py  -i disco_bubbles_coherent.fa [-o disco_bubbles_coherent.vcf -m 0.95 -r 0.4]
    
    *********************************************** '''


import sys
import getopt
import random
import re #regular expressions
import time

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../scripts/"))
from vcf_formatting_functions import *


''' Usages in discoSnp pipeline scripts (RAD/clustering_scripts/discoRAD_clustering.sh):

    If no clustering:
    python3 fasta_and_cluster_to_filtered_vcf.py -i disco_bubbles_coherent.fa -o disco_bubbles_coherent.vcf -m 0.95 -r 0.4
    end of the pipeline
    
    If clustering:
    python3 fasta_and_cluster_to_filtered_vcf.py -i disco_bubbles_coherent.fa -f -o disco_bubbles_coherent_filtered.fa_removemeplease -m 0.95 -r 0.4
    # then clustering with file disco_bubbles_coherent_filtered.fa_removemeplease
    python3 fasta_and_cluster_to_filtered_vcf.py -i disco_bubbles_coherent_filtered.fa_removemeplease -c disco_bubbles_coherent_filtered_simpler.cluster -o disco_bubbles_coherent_clustered.vcf -s 150
    rm -f disco_bubbles_coherent_filtered.fa_removemeplease
    end of pipeline
    '''

def store_clusters(cluster_file):
    if cluster_file==None: return None, None
    clusters=open(cluster_file,"r")
    read_id_to_cluster_id={}
    cluster_id_to_cluster_size={}
    cluster_id=-1
    for cluster in clusters:
        # a line is "70166 70345 70409 70222 70406 70167 70223 69786 70407 69787 70408 70611 70610 70344 "
        cluster_id+=1
        cluster_id_to_cluster_size[cluster_id]=int(len(cluster.rstrip().split())/2)
        for read_id in cluster.rstrip().split():
            read_id_to_cluster_id[int(read_id.split('-')[0])]=cluster_id # A line can be formated as 70166 70345-info_about_similarity
    clusters.close()
    return read_id_to_cluster_id, cluster_id_to_cluster_size
    
def get_cluster_id_and_size(sequence_id, read_id_to_cluster_id, cluster_id_to_cluster_size):
    if not read_id_to_cluster_id: return ".", "."
    if sequence_id not in read_id_to_cluster_id:
        print("Warning, sequence id "+str(sequence_id)+" not in clusters",file=sys.stderr)
        return ".", "."
    return read_id_to_cluster_id[sequence_id], cluster_id_to_cluster_size[read_id_to_cluster_id[sequence_id]]
    

def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]," : discoSnp output filtering and formatting in vcf")
    print("-----------------------------------------------------------------------------")
    print("usage: ",sys.argv[0]," -i disco_bubbles.fa")
    print("  -r: min rank value filter (default = 0)")
    print("  -m: max missing value filter (default = 1)")
    print("  -s: max cluster size filter (default = 0, ie no size limit)")
    print("  -o: output vcf file path (default = stdout)")
    print("  -f: output a filtered fasta file instead of a vcf file")
    print("  -c: considers a cluster input file. In this situation, can filter on cluster size and prints the cluster_id and cluster_size in the INFO field of each variant")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:r:m:s:o:fc:", ["help", "in=", "rank=", "miss=", "size=", "out=", "fastaout", "cluster"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    # Default parameters
    fasta_file = 0
    fasta_only = 0
    min_rank = 0
    max_miss = 1
    max_cluster_size = 0
    k = 31
    out_file =      None
    cluster_file =  None
    with_cluster = False
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
        elif opt in ("-s", "--size"):
            max_cluster_size = int(arg)
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-f", "--fastaout"):
            fasta_only = 1
        elif opt in ("-c"):
            cluster_file = arg
            with_cluster = True
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
        
        ## LOAD clusters
        read_id_to_cluster_id, cluster_id_to_cluster_size = store_clusters(cluster_file)
        
        with open(fasta_file, 'r') as filin:
            for line in filin:
                splitted_1 = line.split("|")

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
            if not fasta_only:
                # Write vcf comment lines
                filout.write(vcf_header(source,date,fasta_file,nb_samples))

            nb_kept_variants = 0
            nb_analyzed_variants = 0
            line_count = 0
            sequence_id=-1
            fasta_4lines = ""  ## Remembering 4 consecutive lines for kept variants to write in a fasta file if fasta_only mode
            keep_variant = False
            cluster_id=-1
            cluster_size=0
            for line in filin:
                if line_count%2 == 0:   sequence_id+=1  # first sequence is 1, second (lower path of first variant) is 2, ...
                if line_count ==0:      
                    fasta_4lines=""                     #back to empty for incoming variant 
                    cluster_id,cluster_size = get_cluster_id_and_size(sequence_id, read_id_to_cluster_id, cluster_id_to_cluster_size)
                line_count += 1
                fasta_4lines += line
                if line_count == 1:
                    keep_variant = False
                    nb_analyzed_variants += 1

                    # Header higher path
                    line = line.strip()
                    splitted_1 = line.split("|")
                    #fasta_4lines = splitted_1[0] + "\n"  #simplified headers for fasta_only and src

                    ## FILTERING
                    #filter cluster size
                    if max_cluster_size >0 and cluster_id != ".":
                        if cluster_size > max_cluster_size:
                            continue
                    #filter rank
                    rank = float(splitted_1[-1].split("rank_")[1])
                    if rank < min_rank:
                        continue
                    # filter missing genotype ratio
                    if max_miss !=1:
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
                    #fasta_4lines += splitted_2[0] + "\n"  #simplified headers for fasta_only and src

                if line_count == 4:
                    line_count = 0
                    if keep_variant:
                        if fasta_only:
                            #fasta_4lines += line   #simplified headers for fasta_only and src
                            filout.write(fasta_4lines) # TODO: do we writte fasta variants if not in a cluster and a cluster file is provided?
                        else:
                            #now format in vcf format
                            line = line.strip()
                            filout.write(format_vcf(splitted_1, splitted_2, nb_samples, rank, line, cluster_id, cluster_size, tig_type))
                    


        #print(f"{nb_lost_variants} variant bubbles filtered out")
        #print(f"{nb_kept_variants} variant bubbles output out of {nb_tot_variants} ({nb_analyzed_variants} analyzed)")
        sys.stderr.write(f"{nb_kept_variants} variant bubbles output out of {nb_analyzed_variants}\n")


if __name__ == "__main__":
    main()
                      


