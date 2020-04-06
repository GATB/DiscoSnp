import sys
"""
From a vcf file (arg1) and a valid dataset_id (arg2) this script extracts the max, min and average DP for the give dataset. 
Returns an error if the dataset_id does not exists in the vcf

Input format (with one unique dataset, else lines would have one or more additional : GT:DP:PL:AD:HQ fields (eg. 0/1:1727:5117,821,20742:1255,472:73,73).
# Comments
ecoli2kb    599    5    A    C    .    PASS    Ty=SNP;Rk=0;UL=19;UR=19;CL=.;CR=.;Genome=A;Sd=-1    GT:DP:PL:AD:HQ    0/1:1727:5117,821,20742:1255,472:73,73
ecoli2kb    749    4    A    C    .    PASS    Ty=SNP;Rk=0;UL=19;UR=1721;CL=.;CR=.;Genome=G;Sd=1    GT:DP:PL:AD:HQ    0/1:1717:23527,1401,3352:353,1364:73,73
ecoli2kb    649    3    A    C    .    PASS    Ty=SNP;Rk=0;UL=19;UR=19;CL=.;CR=.;Genome=T;Sd=1    GT:DP:PL:AD:HQ    0/0:1770:1466,2566,28486:1562,208:73,73
ecoli2kb    699    2    A    C    .    PASS    Ty=SNP;Rk=0;UL=19;UR=19;CL=.;CR=.;Genome=A;Sd=1    GT:DP:PL:AD:HQ    1/1:1712:30704,3553,410:97,1615:73,73
ecoli2kb    549    1    A    C    .    PASS    Ty=SNP;Rk=0;UL=19;UR=518;CL=.;CR=.;Genome=A;Sd=-1    GT:DP:PL:AD:HQ    0/1:1657:22184,1236,3525:361,1296:73,73


Output format (stdout): 
param DP_avg := 1716.6
param DP_max := 1770;
param DP_min := 1657;
"""

def extract_DP(vcf_file_name, dataset_id):
    if dataset_id < 1:
        sys.stderr.write("Error: dataset id must be > 0\n")
        return
    min=sys.maxsize
    max=0
    sum=0
    nb=0
    
    with open(vcf_file_name) as vcf_file:
        # SNP_higher_path_9	22	9	C	G	.	.	Ty=SNP;Rk=0;UL=2;UR=3;CL=.;CR=.;Genome=.;Sd=.	GT:DP:PL:AD:HQ	1/1:540:9155,939,294:48,492:50,50
        for line in vcf_file:
            line = line.strip()
            if line[0] =="#": continue
            s_line = line.split()
            field_id = 8+dataset_id
            try:
                GT_DP = s_line[field_id]
                DP = int(GT_DP.split(":")[1])
                if DP > max: max=DP
                if DP < min: min=DP
                sum+=DP
                nb+=1
            except IndexError: 
                sys.stderr.write("Error: no field nb "+str(dataset_id)+" in file "+vcf_file_name+"\n")
                return
        print(f"param DP_avg := {int(sum/float(nb))};")
        print(f"param DP_max := {max};")
        print(f"param DP_min := {min};")
            
            

def main(vcf_file_name, dataset_id):
    '''
    Extraction of DP min, max and average from a VCf file
    Usage: 
        python ~/workspace/gatb-discosnp/scripts/k3000/extract_DP_from_vcf.py /discoRes_k_31_c_2_D_0_P_3_b_2_coherent.vcf >> graph_5haplotypes_filtered.dat #once  graph_5haplotypes_filtered.dat was created using K3000_gfa_to_dat.py
    '''
    extract_DP(vcf_file_name, dataset_id)



if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
