#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import getopt


#
#input: disco.vcf
#output : nouveau disco.vcf avec un SNP par cluster, celui ayant le moins de génotypes manquants
#
#utilisation : python filter_paralogs.py disco.vcf
#
    
    
def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]," : select one variant from a vcf with clusters")
    print("-----------------------------------------------------------------------------")
    print("usage: ",sys.argv[0])
    print("  -v vcf file [mandatory]")
    print("  -o: output vcf file path (default = stdout)")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)
    
def store_info(vcf_file):

    dict_ = {}  ## {num_cluster : [num_SNP ayant le moins de géno manquants, nb_missgeno]}
    filin = open(vcf_file, 'r')
    nb_SNP_tot = 0

    while True:
        """#SNP_higher_path_14643	30	14643	C	T	.	.	Ty=SNP;Rk=0.55424;UL=0;UR=0;CL=0;CR=0;Genome=.;Sd=.;Cluster=1285;ClSize=4	GT:DP:PL:AD:HQ	0/1:38:554,48,75:7,31:71,71	0/1:20:63,23,263:15,5:71,71"""

        line = filin.readline()
        if not line: break
        if line.startswith("#"): continue

        nb_SNP_tot += 1

        try:
            num_cluster = int(line.split("Cluster=")[1].split(";")[0])
        except ValueError:
                print ("No cluster size information stored in the vcf, exit")
                sys.exit(1)
          
        if num_cluster == -1: continue
        id_SNP = line.split("\t")[2]

        if num_cluster not in dict_:
            dict_[num_cluster] = ["x", sys.maxsize]                             # ghost best variant for this new cluster

        genotypes = [i.split(":")[0] for i in line.split("\t")[9:]]
        nb_missing = 0

        for ind in genotypes:
            if ind[0] == "." : 
                nb_missing += 1
            
       
        if nb_missing >= dict_[num_cluster][1]: continue
        dict_[num_cluster] = [id_SNP, nb_missing]

    filin.close()
    return dict_, nb_SNP_tot, len(dict_)



def output_newvcf(vcf_file, out_file, dict_) :
    
    filin = open(vcf_file, 'r')
    if out_file:
        filout=open(out_file,'w')
    else:
        filout = sys.stdout
        

    while True:

        line = filin.readline() #SNP_higher_path_14643	30	14643	C	T	.	.	Ty=SNP;Rk=0.55424;UL=0;UR=0;CL=0;CR=0;Genome=.;Sd=.;Cluster=1285;ClSize=4	GT:DP:PL:AD:HQ	0/1:38:554,48,75:7,31:71,71	0/1:20:63,23,263:15,5:71,71
        if not line: break
        if line.startswith("#"): 
             filout.write(line)        
             continue

        try:
            cluster = int(line.split("Cluster=")[1].split(";")[0])
        except ValueError:
                print ("No cluster size information stored in the vcf, exit")
                sys.exit(1)
        if cluster == -1 :  continue
 
        id_SNP = line.split("\t")[2]       
        if id_SNP != dict_[cluster][0]: continue
       
        filout.write(line)

    filin.close()
    filout.close()


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv:o:")
    except getopt.GetoptError as err:
        # print help information and exit:
        usage()
        sys.exit(2)
    
    # Default parameters
    vcf_file =       None
    out_file =      None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-v"):
            vcf_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        else:
            assert False, "unhandled option"

    if vcf_file==None:
        print ("-v missing")
        usage()
        sys.exit(2)



    dict_, nb_SNP_before, nb_SNP_after = store_info(vcf_file)
    output_newvcf(vcf_file, out_file, dict_)

    
    
if __name__ == "__main__":
    main()


