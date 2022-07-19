#!/usr/bin/env python3
# -*- coding: utf-8 -*-


''' ***********************************************

Script to select only one variant per locus (cluster) in a discoSnp vcf output file (.vcf)
Author - Claire Lemaitre, Pierre Peterlongo, Inria

Usage:
python3 1SNP_per_cluster.py -i vcf_file [-o new_vcf_file]

Details:
The selected variant is not chosen at random, for a given cluster, it is the one with the less missing genotypes (if ties, it is the first read in the input file)

*********************************************** '''

import sys
import getopt


    
def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]+" : selects one variant per cluster")
    print("-----------------------------------------------------------------------------")
    print("usage: "+sys.argv[0]+" -i vcf_file [-o output_file]")
    print("  -i: vcf file [mandatory]")
    print("  -o: output vcf file (default = stdout)")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)
    

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))
        usage()
        sys.exit(2)

    # Default parameters
    vcf_file =       None
    out_file =      None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i"):
            vcf_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        else:
            assert False, "unhandled option"

    if vcf_file==None:
        print ("Error: option -i is mandatory")
        usage()
        sys.exit(2)

    format_ok = check_format(vcf_file)
    if not format_ok:
        print("Error: the format of the input vcf is not correct, it must contain clustering information")
        sys.exit(2)
    
    dict_, nb_SNP_before, nb_SNP_after = store_info(vcf_file)
    output_newvcf(vcf_file, out_file, dict_)




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

        num_cluster = int(line.split("Cluster=")[1].split(";")[0])
          
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


    
    
if __name__ == "__main__":
    main()


