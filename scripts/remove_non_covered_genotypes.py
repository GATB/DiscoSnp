import sys
import gzip

if len(sys.argv)<3:
    print ("This tool replaces discoSnp VCF genotypes with DP lower or equal to a threshold to \"./.\"")
    print ("python remove_non_covered_genotypes.py \".vcf from discoSnp\" \"DP threshold\"")
    sys.exit()


if "gz" in sys.argv[1]:
    coherent_file=gzip.open(sys.argv[1],"r")
else: 
    coherent_file=open(sys.argv[1],"r")
dp_threshold = float(sys.argv[2])

while True:
    
    
    line = coherent_file.readline().rstrip()
    if not line: break
    if line[0]=='#': 
        print (line)
        continue
    
#    SNP_higher_path_3       196     3       C       G       .       .       Ty=SNP;Rk=1;UL=86;UR=261;CL=166;CR=761;Genome=.;Sd=.    GT:DP:PL:AD:HQ  0/0:124:10,378,2484:124,0:0,0   1/1:134:2684,408,10:0,134:0,0
    splitted_line = line.split()
    
    toprint=""
    for i in range(9):
        toprint+=splitted_line[i]+'\t'
    for i in range(9,len(splitted_line)):
        splitted_geno = splitted_line[i].split(':')
        if int(splitted_geno[1])<=dp_threshold: toprint+= "./.:"
        else: toprint+= splitted_geno[0]+":"
        
        for j in range(1,len(splitted_geno)):
            toprint+= splitted_geno[j]
            if j<len(splitted_geno)-1: toprint+= ":"
        if i<len(splitted_line): toprint+='\t'
    print (toprint)