#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' ***********************************************

Script to filter out variants in a discoSnp vcf output file (.vcf), according to various features
Author - Claire Lemaitre, Pierre Peterlongo, Inria

Usage:
python3 filter_vcf_by_indiv_cov_max_missing_and_maf.py -i vcf_file -o new_vcf_file [-x 0.1 -y 0.5]

Details:
filter a vcf file by keeping only variants such that :
  - each individual genotype with DP < min_cov is replaced by a missing genotype
  - the proportion of missing genotype is < max_missing
  - the minor allele freq is >= min_maf
outputs a vcf
Note : NA are indicated as "./." genotypes.
*********************************************** '''


# Last modified : Oct. 2019

import sys
import getopt


def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]+" : vcf filter")
    print("-----------------------------------------------------------------------------")
    print("usage: "+sys.argv[0]," -i vcf_file -o output_file [-c min_cov -m max_missing -f min_maf -s]")
    print("  -i: input vcf file [mandatory]")
    print("  -o: output vcf file [mandatory]")
    print("  -m: max missing genotype proportion to keep a variant (between 0 and 1, def = 1)")
    print("  -c: min coverage to call a genotype (int, def=0)")
    print("  -f: min minor allele frequency (maf) (between 0 and 1, def=0)")
    print("  -s: snp only (def= all variants)")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:c:m:f:s", ["in=", "out=","cov=","miss=", "freq=","snp-only"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    # Default parameters
    output_file = 0
    input_file = 0
    min_cov = 0
    max_missing_prop = 1
    snp_only = 0
    min_maf = 0
    for opt, arg in opts:
        if opt in ("-i", "--in"):
            input_file = arg
        elif opt in ("-o", "--out"):
            output_file = arg
        elif opt in ("-c", "--cov"):
            min_cov = int(arg)
        elif opt in ("-m", "--miss"):
            max_missing_prop = float(arg)
        elif opt in ("-f", "--freq"):
            min_maf = float(arg)
        elif opt in ("-s", "--snp-only"):
            snp_only=1
        else:
            assert False, "unhandled option"

    if input_file == 0 or output_file == 0:
        print("Error: options -i and -o are mandatory")
        usage()
        sys.exit(2)
    else:

        filin = open(input_file, 'r')
        filout = open(output_file, 'w')
    
        n_geno = 0
        max_missing = 0 # integer
        var_count = 0
        var_written = 0
        na_count = 0
        new_na_count = 0
        selected_na_count = 0
    
        while True:
            #random_genome	49694	9997	G	T	.	MULTIPLE	Ty=SNP;Rk=1;UL=83;UR=4;CL=83;CR=4;Genome=G;Sd=1	GT:DP:PL:AD:HQ
            line = filin.readline()
            if not line: break
            if line.startswith("#"):
                filout.write(line)
                if line.startswith("#CHROM"):  # Used to record the total number of samples
                    n_geno=len(line.split("\t")) - 9
                    max_missing = int(max_missing_prop*n_geno)  # int -> floor(), ensures proportion of non na is greater or equal to 1-max_missing_prop (when using miss <= max_missing)
                continue
 
            var_count += 1
 
            # Filtering on the type of variant (SNP and INDELs or SNPs only)
            # WARNING : designed for discoSNP only, should remain compatible with Stacks (no error, but will not filter out INDEL)
            thistype = line.split("\t")[7].split(";")[0].split("=")[1].strip()
            #print(thistype)
            if snp_only and thistype == "INDEL": continue
     
            splitted = line.split("\t")
            start_geno = 9

            line_towrite = "\t".join(splitted[:start_geno])
            missing_count = 0
            
            ref_count = 0
            alt_count = 0

            for geno in splitted[start_geno:]:
                geno_info = geno.split(":")
                genotype = geno_info[0]
                if genotype == "./." or genotype == ".|.":
                    missing_count += 1
                    na_count += 1
                else:
                    DP = int(geno_info[1])
                    if DP < min_cov:
                        genotype = "./."
                        missing_count += 1
                        new_na_count += 1
                    else:
                        #count ref and alt alleles
                        allele1 = genotype[0]
                        allele2 = genotype[2]
                        if allele1 == "0":
                            ref_count += 1
                        else:
                            alt_count += 1
                        if allele2 == "0":
                            ref_count += 1
                        else:
                            alt_count += 1
                        
                #rewrite geno info
                geno_new = genotype+":"+":".join(geno_info[1:])
                line_towrite += "\t"+geno_new

            maf = 0
            # Should we also remove variants that are no longer variable if min_maf = 0, ie alt_count == 0 and ref_count ==0 ??
            if ref_count+alt_count>0:
                maf = (float) (min(ref_count,alt_count))/(ref_count+alt_count)
                ## Filter on missing count and maf
                if missing_count <= max_missing and maf >= min_maf:
                    filout.write(line_towrite)
                    var_written += 1
                    selected_na_count += missing_count

        filin.close()
        filout.close()
        print("# "+sys.argv[0])
        print("# input_file : "+input_file)
        print("# output_file : "+output_file)
        print("# filter parameters : indiv_DP>="+str(min_cov)+", missing<="+str(max_missing)+" (prop*n_geno = "+str(max_missing_prop)+" * "+str(n_geno)+"), snp-only="+str(snp_only))
        print("# "+str(var_count)+" seen variants, "+str(var_written)+" variants after filtering")
        print("# initial missing count (on snps only if snp-only=1) = "+str(na_count)+", genotypes changed to NA = "+str(new_na_count)+", final missing count in new vcf = "+str(selected_na_count))
        print("# initial missing percent (on snps only if snp-only=1) = {0:.2f} %, final missing percent in new vcf = {1:.2f} %".format(100*float(na_count)/float(var_count*n_geno),100*float(selected_na_count)/float(var_written*n_geno)))



if __name__ == "__main__":
    main()



