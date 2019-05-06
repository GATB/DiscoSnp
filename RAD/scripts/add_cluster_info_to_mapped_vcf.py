
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import getopt


def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]," : merge 2 vcf files, add cluster info present in the unmmaped vcf in the mapped vcf (ie. with mapping info on a ref genome)")
    print("-----------------------------------------------------------------------------")
    print("usage: ",sys.argv[0])
    print("  -u unmapped vcf file [mandatory]")
    print("  -m mapped vcf file [mandatory]")
    print("  -o: output vcf file path (default = stdout)")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)


def output_newvcf(unmapped_file, mapped_file, out_file):

    filin = open(unmapped_file, 'r')
    ## store cluster_info for each variant (identified by its unique id, 3rd column of the vcf line)
    id_to_cluster_info={}
    for line in filin.readlines():
        line=line.strip()
        if line[0]=="#": continue
    #SNP_higher_path_3       199     3       C       G       .       .       Ty=SNP;Rk=1.0;UL=86;UR=261;CL=169;CR=764;Genome=.;Sd=.;Cluster=0;ClSize=3  ...
        splitted = line.split("\t")
        cluster_info = ";".join(splitted[7].split(";")[8:10])
        id = splitted[2]
        id_to_cluster_info[id] = cluster_info
        #print(f"{id}-{cluster_info};")
    
    filin.close()
    
    ## Print vcf
    if out_file:
        filout=open(out_file,'w')
    else:
        filout = sys.stdout
        
    filin = open(mapped_file,"r")
    for line in filin.readlines():
        if line[0]=="#":
            filout.write (line)
        else:
            ##chr3R    26778135     3       C       G       .       PASS       Ty=SNP;Rk=1.0;UL=86;UR=261;CL=169;CR=764;Genome=.;Sd=1  ...
            splitted = line.split("\t")
            id = splitted[2]
            cluster_info = id_to_cluster_info[id]
            INFO = splitted[7] + ";" + cluster_info
            tojoin = splitted[:7] + [INFO] + splitted[8:]
            filout.write ("\t".join(tojoin))

    filin.close()
    filout.close()

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hu:m:o:")
    except getopt.GetoptError as err:
        # print help information and exit:
        usage()
        sys.exit(2)
    
    # Default parameters
    unmapped_file = None
    mapped_file = None
    out_file = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-u"):
            unmapped_file = arg
        elif opt in ("-m"):
            mapped_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        else:
            assert False, "unhandled option"

    if mapped_file==None:
        print ("-m missing")
        usage()
        sys.exit(2)

    if unmapped_file==None:
        print ("-u missing")
        usage()
        sys.exit(2)
    output_newvcf(unmapped_file, mapped_file, out_file)
    
    
if __name__ == "__main__":
    main()

