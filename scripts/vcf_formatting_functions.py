#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 


''' ***********************************************
    
    Functions to format discoSnp raw output file (.fa) into a vcf format file (.vcf) (ghost vcf)
    Author - Claire Lemaitre
    
    Used by create_filtered_vcf.py and fasta_and_cluster_to_filtered_vcf.py
    
    *********************************************** '''

import re

def vcf_header(source, date, fasta_file, nb_samples):
    '''
        Returns the first lines of the vcf file, composed of the COMMENTS and the HEADER
        source:     the name of the python program that writes this vcf
        date:       the date to appear in the vcf file
        fasta_file: the input fasta file
        nb_samples: the nb of samples
        '''

    VCF_COMMENTS = f'''##fileformat=VCFv4.1
##filedate={date}
##source={source}
##SAMPLE=file://{fasta_file}
##REF=<ID=REF,Number=1,Type=String,Description="Allele of the path Disco aligned with the least mismatches">
##FILTER=<ID=MULTIPLE,Description="Mapping type : PASS or MULTIPLE or .">
##INFO=<ID=Ty,Number=1,Type=String,Description="SNP, INS, DEL or .">
##INFO=<ID=Rk,Number=1,Type=Float,Description="SNP rank">
##INFO=<ID=UL,Number=1,Type=Integer,Description="length of the unitig left">
##INFO=<ID=UR,Number=1,Type=Integer,Description="length of the unitig right">
##INFO=<ID=CL,Number=1,Type=Integer,Description="length of the contig left">
##INFO=<ID=CR,Number=1,Type=Integer,Description="length of the contig right">
##INFO=<ID=Genome,Number=1,Type=String,Description="Allele of the reference;for indel reference is . ">
##INFO=<ID=Sd,Number=1,Type=Integer,Description="Reverse (-1) or Forward (1) Alignement">
##INFO=<ID=Cluster,Number=1,Type=Integer,Description="Cluster ID when variants have been clustered (RAD-seq mode)">
##INFO=<ID=ClSize,Number=1,Type=Integer,Description="Cluster size (RAD-seq mode)">
##INFO=<ID=XA,Number=.,Type=String,Description="Other mapping positions (chromosome_position). Position is negative in case of Reverse alignment. The position designs the starting position of the alignment, not the position of the variant itself.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled Genotype Likelihoods">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Depth of each allele by sample">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">'''

    HEADER = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]+["G"+str(i+1) for i in range(nb_samples)])
    
    header_str = VCF_COMMENTS + "\n" + HEADER + "\n"
    return header_str


def format_vcf(splitted_1, splitted_2, nb_samples, rank, sequence, cluster_id, cluster_size, tig_type = 0):
    '''
        Format information extracted from the fasta headers of a single bubble into one or several vcf lines

        splitted_1:     a list of values contained in the higher path fasta header (splitted by "|")
        splitted_2:     idem of lower path
        nb_samples:     nb of samples
        rank:           rank of the variant
        sequence :      DNA sequence of the lower path : usefull for INDELs (longest sequnce)
        cluster_id:     cluster_id, if not clustered "."
        cluster_size:   size of the cluster, if not clustered "."
        tig_type :      0 no extension, 1 unitig, 2 contig (ie. unitig and contig length are present)

        '''
    
    '''
        >SNP_higher_path_1487|P_1:30_A/G|low|nb_pol_1|left_unitig_length_
        SNP_higher_path_1487    34    1487    A    G    .    .    Ty=SNP;Rk=0.00072476;UL=4;UR=60;CL=.;CR=.;Genome=.;Sd=.;Cluster=.;ClSize=.    GT:DP:PL:AD:HQ    0/1:53:194,32,613:37,16:71,71    0/1:83:296,44,955:58,25:71,71
        '''
    
    ''' INDEL :
        >INDEL_higher_path_3205|P_1:30_27_3|high|nb_pol_1|left_unitig_length_14|right_unitig_length_6|C1_14|C2_17|Q1_71|Q2_71|G1_0/1:243,13,203|G2_0/1:268,13,248|rank_0.019011
        aaggcagcggccagTCCAGGATGTCCAAGAATTCAACCAATTCGAACAATTCTAAAGGATCGTTTAATTCAagcggc
        >INDEL_lower_path_3205|P_1:30_27_3|high|nb_pol_1|left_unitig_length_14|right_unitig_length_6|C1_16|C2_18|Q1_71|Q2_71|G1_0/1:243,13,203|G2_0/1:268,13,248|rank_0.019011
        aaggcagcggccagTCCAGGATGTCCAAGAATTCAACCAATTCG GGACAGTCCAGATAGTCGTATAACTCG AACAATTCTAAAGGATCGTTTAATTCAagcggc
        aaggcagcggccagTCCAGGATGTCCAAGAATTCAACCAAT TCGGGACAGTCCAGATAGTCGTATAAC TCGAACAATTCTAAAGGATCGTTTAATTCAagcggc
        INDEL_higher_path_3205    41    3205    T    TTCGGGACAGTCCAGATAGTCGTATAAC    .    .    Ty=INS;Rk=0.019011;UL=14;UR=6;CL=.;CR=.;Genome=.;Sd=.;Cluster=.;ClSize=.    GT:DP:PL:AD:HQ    0/1:30:243,13,203:14,16:71,71    0/1:35:268,13,248:17,18:71,71
        
        POS = UL + POS  : 1-based   (note: no longer left-normalized)
        '''
    
    FORMAT = "GT:DP:PL:AD:HQ"
    nb_fixed_fields = 4 + 2*tig_type

    vcf_line = ""
    nb_pol = int(splitted_1[3].split("nb_pol_")[1])
    path_name = splitted_1[0].lstrip(">")
    CHROM = path_name
    id = path_name.split("_")[3]
    ty = re.findall("(\w+)_higher_path",path_name)[0]
    
    if ty != "SNP":
        ty = "INS"
    
    INFO = f"Ty={ty};Rk={rank};"
    position_offset = 0
    if tig_type >0:
        unitig_len = re.findall("_unitig_length_(\d+)","|".join(splitted_1[4:6]))
        position_offset = int(unitig_len[0])
        INFO += f"UL={unitig_len[0]};UR={unitig_len[1]};";
        if tig_type == 2:
            contig_len = re.findall("_contig_length_(\d+)","|".join(splitted_1[6:8]))
            position_offset = int(contig_len[0]) # if contig POS = POS + left_contig_len
            INFO += f"CL={contig_len[0]};CR={contig_len[1]};"
        else:
            INFO += "CL=.;CR=.;"
    else:
        INFO += "UL=.;UR=.;CL=.;CR=.;"
    INFO += f"Genome=.;Sd=.;Cluster={cluster_id};ClSize={cluster_size}"

    # Genotype info (same for all polymorphisms)
    GENO = ""
    for indiv in range(nb_samples):
        ad_1 = splitted_1[nb_fixed_fields + indiv].split("_")[1]
        ad_2 = splitted_2[nb_fixed_fields + indiv].split("_")[1]
        dp = int(ad_1) + int(ad_2)
        qual_1 = splitted_1[nb_fixed_fields + nb_samples + indiv].split("_")[1]
        qual_2 = splitted_2[nb_fixed_fields + nb_samples + indiv].split("_")[1] #necessary : is qual_2 always == qual_1 ???
        
        geno_fields = splitted_1[nb_fixed_fields + 2*nb_samples + indiv].split("_")[1].split(":")
        genotype = geno_fields[0]
        likelihood = geno_fields[1]
        GENO += "\t" + ":".join([genotype,str(dp),likelihood,",".join([ad_1,ad_2]),",".join([qual_1,qual_2])])

    GENO = GENO.strip()

    cigar = splitted_1[1].split(",")
    isolated = False
    if len(cigar) == 1:
        isolated = True

    i = 1
    for pol in cigar:
        if ty == "SNP":
            POS, REF, ALT = re.findall("P_\d+:(\d+)_(\w)/(\w)",pol)[0]
            POS = int(POS) + position_offset  # POS is 1-based
        else:  # INDEL
            POS, indel_size = re.findall("P_\d+:(\d+)_(\d+)",pol)[0]
            POS = int(POS) + position_offset # 1-based
            ALT = sequence[(POS-1):(POS+int(indel_size))]
            REF = ALT[0]
        ID = id
        if not isolated:
            ID += f"_{i}"
        my_line = "\t".join([CHROM, str(POS), ID, REF, ALT, ".", ".", INFO, FORMAT, GENO])
        vcf_line += my_line + "\n"
        i += 1

    return vcf_line

    
