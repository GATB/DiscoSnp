import sys
import gzip

if len(sys.argv)<2:
    print ("This tool removes from discoSnp sorted VCF the locus which are tri-allelic or more")
    print ("python remove_non_diploids.py \".vcf from discoSnp\" ")
    sys.exit()


if "gz" in sys.argv[1]:
    coherent_file=gzip.open(sys.argv[1],"r")
else: 
    coherent_file=open(sys.argv[1],"r")

previous_line=""
previous_chr=""
previous_pos=""
previous_triploid=False

while True:    
    line = coherent_file.readline().rstrip()
    if not line: break
    if line[0]=='#': 
        print (line)
        continue
#    SNP_higher_path_3       196     3       C       G       .       .       Ty=SNP;Rk=1;UL=86;UR=261;CL=166;CR=761;Genome=.;Sd=.    GT:DP:PL:AD:HQ  0/0:124:10,378,2484:124,0:0,0   1/1:134:2684,408,10:0,134:0,0
    splitted_line = line.split()
    chr=splitted_line[0]
    pos=splitted_line[1]
    if chr==previous_chr and pos==previous_pos:
        previous_triploid=True
        continue
    else:
        if not previous_triploid and len(previous_line)>0:
            print (previous_line)
    previous_chr=chr
    previous_pos=pos
    previous_triploid=False
    previous_line=line
