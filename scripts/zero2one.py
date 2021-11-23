import argparse
import fileinput

def zero2one(input_vcf_file_name: str):
    with fileinput.FileInput(input_vcf_file_name, inplace = True) as vcf_file:
        for vcf_line in vcf_file:
            if vcf_line[0] == "#":
                print(vcf_line, end='')
                continue
            # From
            #121_cassette_SSV9_CAG-eGFP-ISce-PASV40      3163     8622   A       C       .       MULTIPLE        Ty=SNP;Rk=8.0508e-05;UL=0;UR=0;CL=0;CR=0;Genome=A;Sd=1;XA=121_cassette_SSV9_CAG-eGFP-ISce-PASV40_17,121_cassette_SSV9_CAG-eGFP-ISce-PASV40_3081,121_cassette_SSV9_CAG-eGFP-ISce-PASV40_79  GT:DP:PL:AD:HQ  ./.:2:.,.,.:2,0:0,0     0/0:291967:8727,869383,5824307:291691,276:71,48
            # To 
            #121_cassette_SSV9_CAG-eGFP-ISce-PASV40      3163     8622   A       C       .       MULTIPLE        Ty=SNP;Rk=8.0508e-05;UL=0;UR=0;CL=0;CR=0;Genome=A;Sd=1;XA=121_cassette_SSV9_CAG-eGFP-ISce-PASV40_18,121_cassette_SSV9_CAG-eGFP-ISce-PASV40_3082,121_cassette_SSV9_CAG-eGFP-ISce-PASV40_80  GT:DP:PL:AD:HQ  ./.:2:.,.,.:2,0:0,0     0/0:291967:8727,869383,5824307:291691,276:71,48
            vcf_line = vcf_line.strip().split()
            new_POS = int(vcf_line[1]) + 1
            ID_REF_ALT_QUAL_FILTER = '\t'.join(vcf_line[2:7])
            INFO = vcf_line[7]
            # zero to one for XA positions
            s_INFO = INFO.split(';')
            if s_INFO[-1].startswith("XA="):
                XA="XA="
                all_XAs = s_INFO[-1].split("=")[-1].split(',')
                for current_XA in all_XAs:
                    pos = int(current_XA.split("_")[-1]) + 1
                    XA += "_".join(current_XA.split("_")[:-1])+"_"+str(pos)+","
                XA = XA[:-1] # remove additional ','
                INFO = ';'.join(s_INFO[:-1]) + ";" + XA
                    
                    
            end_line = '\t'.join(vcf_line[8:])
            new_vcf_line = f"{vcf_line[0]}\t{new_POS}\t{ID_REF_ALT_QUAL_FILTER}\t{INFO}\t{end_line}"
            print(new_vcf_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="name of input zero-based vcf file (will be overwritten!)", dest='input_vcf_file', type=str, required=True)
    args = parser.parse_args()
    zero2one(args.input_vcf_file)
    print(f"Transformed zero-based {args.input_vcf_file} to one-based")
