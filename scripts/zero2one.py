import argparse
import fileinput

def zero2one(input_vcf_file_name: str):
    with fileinput.FileInput(input_vcf_file_name, inplace = True) as vcf_file:
        for vcf_line in vcf_file:
            if vcf_line[0] == "#":
                print(vcf_line, end='')
                continue
            # From
            # SNP_higher_path_3       196     3       C       G       .       .       Ty=SNP;Rk=1;UL=86;UR=261;CL=166;CR=761;Genome=.;Sd=.    GT:DP:PL:AD:HQ  0/0:124:10,378,2484:124,0:0,0   1/1:134:2684,408,10:0,134:0,0
            # To 
            # SNP_higher_path_3       197     3       C       G       .       .       Ty=SNP;Rk=1;UL=86;UR=261;CL=166;CR=761;Genome=.;Sd=.    GT:DP:PL:AD:HQ  0/0:124:10,378,2484:124,0:0,0   1/1:134:2684,408,10:0,134:0,0
            vcf_line = vcf_line.strip().split()
            new_POS = int(vcf_line[1]) + 1
            end_line = '\t'.join(vcf_line[2:])
            new_vcf_line = f"{vcf_line[0]}\t{new_POS}\t{end_line}"
            print(new_vcf_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="name of input zero-based vcf file (will be overwritten!)", dest='input_vcf_file', type=str, required=True)
    args = parser.parse_args()
    zero2one(args.input_vcf_file)
    print(f"Transformed zero-based {args.input_vcf_file} to one-based")
