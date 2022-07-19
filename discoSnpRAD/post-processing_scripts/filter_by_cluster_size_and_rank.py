#!/usr/bin/env python3
# -*- coding: utf-8 -*-


''' ***********************************************

Script to filter out variants in a discoSnp vcf output file (.vcf), according to cluster size and/or rank value
Author - Claire Lemaitre, Pierre Peterlongo, Inria

Usage:
python3 filter_by_cluster_size_and_rank.py -i vcf_file [-o output_file -m 0 -M 150 -r 0.4]

Details:
filter a vcf file by keeping only variants such that :
  - the cluster (locus) they belong to contains x variants with x in [m,M]
  - with a rank value >= r
outputs a vcf
*********************************************** '''


import sys
import getopt


def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]+" : filter vcf variants based on their cluster size and/or rank value")
    print("-----------------------------------------------------------------------------")
    print("usage: "+sys.argv[0]+" -i vcf_file [-o output_file -m 0 -M 150 -r 0.4]")
    print("  -i: input vcf file [mandatory]")
    print("  -m: min cluster size (included)")
    print("  -M: max cluster size (included)")
    print("  -r; min rank (included)")
    print("  -o: output vcf file (default = stdout)")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)


def output_newvcf(in_file, out_file, min_cluster_size, max_cluster_size, rank_min):

    filin = open(in_file, 'r')
    if out_file:
        filout=open(out_file,'w')
    else:
        filout = sys.stdout

    for line in filin.readlines():
        line=line.strip()
        if line[0]=='#': 
            filout.write(line+"\n")
            continue
        
        #SNP_higher_path_3       199     3       C       G       .       .       Ty=SNP;Rk=1.0;UL=86;UR=261;CL=169;CR=764;Genome=.;Sd=.;Cluster=0;ClSize=3  ...
        if min_cluster_size>0 or max_cluster_size<sys.maxsize:          # make the split only if necessary
            try:
                cluster_size = int(line.split("ClSize=")[1].split()[0])
            except ValueError:
                print ("No cluster size information stored in the vcf, exit")
                sys.exit(1)
            if cluster_size < min_cluster_size or cluster_size > max_cluster_size: 
                continue # does not respect the cluster size filtering criteria
        
        if rank_min > 0:                                                # make the split only if necessary
            rank = float(line.split()[7].split(';')[1].split('=')[-1])  
            if rank < rank_min: 
                continue # does not respect the rank min filtering criteria
        
        filout.write(line+"\n") # all filters passed

    filin.close()
    filout.close()



def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:m:M:o:r:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))
        usage()
        sys.exit(2)
    
    # Default parameters
    min_cluster_size = 0
    max_cluster_size = sys.maxsize
    min_rank=0
    in_file = None
    out_file = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i"):
            in_file = arg
        elif opt in ("-m"):
            min_cluster_size =  float(arg)
        elif opt in ("-M"):
            max_cluster_size =  float(arg)
        elif opt in ("-r"):
            min_rank =          float(arg)
        elif opt in ("-o", "--out"):
            out_file = arg
        else:
            assert False, "unhandled option"

    if in_file==None:
        print("Error: option -i is mandatory")
        usage()
        sys.exit(2)
    
    format_ok = check_format(in_file)
    if not format_ok:
        print("Error: the format of the input vcf is not correct, it must contain clustering information")
        sys.exit(2)
    output_newvcf(in_file, out_file, min_cluster_size, max_cluster_size, min_rank)
    

def check_format(vcf_file):
    ''' Checks if the vcf has the correct format, ie : the INFO field must contain clustering information, such as:
        Ty=SNP;Rk=1;UL=1;UR=2;CL=.;CR=.;Genome=.;Sd=.;Cluster=79466;ClSize=12
        '''
    
    filin = open(vcf_file, 'r')
    checked = False
    while not checked:
        line = filin.readline()
        if line.startswith("#"): continue
        INFO_split = line.split("\t")[7].split(";")
        checked = True
        if len(INFO_split) < 10: return False
        tmp_cluster = INFO_split[8].split("Cluster=")
        if len(tmp_cluster) < 2: return False
        if tmp_cluster[1] == ".": return False
        try:
            cl_id = int(tmp_cluster[1])
        except ValueError:
            return False
        return True
            
    filin.close()

if __name__ == "__main__":
    main()
