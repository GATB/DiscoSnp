
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import getopt


def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]," : filter vcf variants based on their cluster size")
    print("-----------------------------------------------------------------------------")
    print("usage: ",sys.argv[0]," -i vcf file")
    print("  -m: min cluster size (included)")
    print("  -M: max cluster size (included)")
    print("  -o: output vcf file path (default = stdout)")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)


def output_newvcf(in_file, out_file, min_cluster_size, max_cluster_size):

    filin = open(in_file, 'r')
    if out_file:
        filout=open(out_file,'w')
    else:
        filout = sys.stdout

    for line in filin.readlines():
        line=line.strip()
        if line[0]=='#': 
            print(line)
            continue
        
    #SNP_higher_path_3       199     3       C       G       .       .       Ty=SNP;Rk=1.0;UL=86;UR=261;CL=169;CR=764;Genome=.;Sd=.;Cluster=0;ClSize=3  ...
        cluster_size = int(line.split("ClSize=")[1].split()[0])
        if cluster_size >= min_cluster_size and cluster_size <= max_cluster_size: 
            print(line)

    filin.close()
    filout.close()



def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:m:M:o:")
    except getopt.GetoptError as err:
        # print help information and exit:
        usage()
        sys.exit(2)
    
    # Default parameters
    min_cluster_size = 0
    max_cluster_size = sys.maxsize
    in_file = None
    out_file = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i"):
            in_file = arg
        elif opt in ("-m"):
            min_cluster_size = float(arg)
        elif opt in ("-M"):
            max_cluster_size = float(arg)
        elif opt in ("-o", "--out"):
            out_file = arg
        else:
            assert False, "unhandled option"

    if in_file==None:
        usage()
        sys.exit(2)
    output_newvcf(in_file, out_file, min_cluster_size, max_cluster_size)
    
    
if __name__ == "__main__":
    main()
