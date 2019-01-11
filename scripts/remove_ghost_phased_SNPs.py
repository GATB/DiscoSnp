import sys

# From a phased file produced by kissreads, this script removes SNPs for which only one of the two allele is phased, and thus exists in the file. 
# It updates the corresponding lines, and remove lines containing zero or one SNP (including pairend facts). 

def index_SNPs(phased_snps_file_name):
    phased_snps_file = open(phased_snps_file_name)
    snps2allele={}
    for line in phased_snps_file:
        # -1000058h_0;2436173h_8;3631958h_202;  -1000183h_0;3224567l_287; => 3
        line=line.replace(" ", "").rstrip().split("=")[0].rstrip(';').split(";") # keep only a list of ids
        for id_snp in line: 
            if id_snp[0]==('-'): id_snp=id_snp[1:]
            id_snp=id_snp.split('_')[0]
            id_snp_only=id_snp[:-1]
            allele=id_snp[-1]
            if id_snp_only not in snps2allele:
                snps2allele[id_snp_only]=set()
            snps2allele[id_snp_only].add(allele)
    phased_snps_file.close()
        
    biallelic_snps=set()
    
    for id_snp in snps2allele:
        if id_snp=="2791159":
            print (snps2allele[id_snp])
        if len(snps2allele[id_snp])==2:
            biallelic_snps.add(id_snp)
    
    return biallelic_snps

def print_phased_snps(phased_snps, biallelic_snps):
    """ Takes a string of phased snps (eg '-1000058h_0;2436173h_8;3631958h_202;') and 
    print an updated fact removing non biallelic ones.
    For instance, if 2436173 is not biallelic, we'd print '-1000058h_0;3631958h_210;'
    For instance, if 1000058 is not biallelic, we'd print '2436173h_0;3631958h_202;'
    """
    # print("print_phased_snps ",phased_snps)
    correction_distance_to_previous=0
    first_to_print=True
    for snp_info in phased_snps.split(';')[:-1]:
        snp_id_only=snp_info.split('_')[0][:-1]
        if snp_id_only[0]=='-': snp_id_only=snp_id_only[1:]
        distance_to_previous=int(snp_info.split('_')[-1])
        if snp_id_only not in biallelic_snps:
            correction_distance_to_previous+=distance_to_previous
        else:
            if first_to_print:
                print (snp_info.split('_')[0]+'_0;',end='')
                first_to_print=False
            else:
                print (snp_info.split('_')[0]+'_'+str(correction_distance_to_previous+distance_to_previous)+';',end='')
            correction_distance_to_previous=0

                

def print_biallelic_snps(phased_snps_file_name,biallelic_snps):
    # print ("# SNPs having only one of their allele present in the original file were removed ")
    phased_snps_file = open(phased_snps_file_name)
    
    for line in phased_snps_file:
        # print("line is", line)
        if line[0]=="#":
            print (line.rstrip())
            continue
        # case1: non paired: 
        # -1000058h_0;2436173h_8;3631958h_202; => 3
        # case2: paired:
        # -1000058h_0;2436173h_8;3631958h_202;  -1000183h_0;3224567l_287; => 3
        
        # first we check that there exists more than one biallelict SNPs in this fact: 
        snp_ids=line.replace(" ", "").rstrip().split("=")[0].rstrip(';').split(";") # keep only a list of ids
        biallelics=set()
        for id_snp in snp_ids: 
            if id_snp[0]==('-'): id_snp=id_snp[1:]
            id_snp=id_snp.split('_')[0]
            id_snp_only=id_snp[:-1]
            if id_snp_only in biallelic_snps:
                biallelics.add(id_snp_only)
        #print(biallelics)
        if len(biallelics)<2: continue 
        
        
        #print first set of phased SNPs:
        print_phased_snps(line.split()[0],biallelic_snps)
        if len(line.split())==4: 
            print(" ",end='')
            print_phased_snps(line.split()[1],biallelic_snps)
        print(" => "+line.split()[-1])
    phased_snps_file.close()
    
biallelic_snps = index_SNPs(sys.argv[1])
print_biallelic_snps(sys.argv[1],biallelic_snps)
