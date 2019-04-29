#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 


''' ***********************************************
    
    Script to filter and format discoSnp raw output file (.fa) into a vcf format file (.vcf)
    Author - Claire Lemaitre
    
    Usage:
    python3 create_filtered_vcf.py -i disco_bubbles_coherent.fa [-o disco_bubbles_coherent.vcf -m 0.95 -r 0.4]
    
    *********************************************** '''


import sys
import getopt
import random
import re #regular expressions
import time

''' Usages in discoSnp pipeline scripts (RAD/scripts/discoRAD_finalization.sh):

    If no clustering:
    python3 create_filtered_vcf.py -i disco_bubbles_coherent.fa -o disco_bubbles_coherent.vcf -m 0.95 -r 0.4
    end of the pipeline
    
    If clustering:
    python3 create_filtered_vcf.py -i disco_bubbles_coherent.fa -f -o disco_bubbles_coherent_filtered.fa_removemeplease -m 0.95 -r 0.4
    clustering with file disco_bubbles_coherent_filtered.fa_removemeplease
    python3 create_filtered_vcf.py -i disco_bubbles_coherent_filtered_clustered.fa_removemeplease -o disco_bubbles_coherent_clustered.vcf
    rm -f disco_bubbles_coherent_filtered.fa_removemeplease
    end of pipeline
    '''

def store_clusters(cluster_file):
    if cluster_file==None: return None, None
    clusters=open(cluster_file,"r")
    read_id_to_cluster_id={}
    cluster_id_to_cluster_size={}
    cluster_id=-1
    for cluster in clusters:
        # a line is "70166 70345 70409 70222 70406 70167 70223 69786 70407 69787 70408 70611 70610 70344 "
        cluster_id+=1
        cluster_id_to_cluster_size[cluster_id]=int(len(cluster.rstrip().split())/2)
        for read_id in cluster.rstrip().split():
            read_id_to_cluster_id[int(read_id.split('-')[0])]=cluster_id # A line can be formated as 70166 70345-info_about_similarity
    clusters.close()
    return read_id_to_cluster_id, cluster_id_to_cluster_size
    
def get_cluster_id_and_size(sequence_id, read_id_to_cluster_id, cluster_id_to_cluster_size):
    if not read_id_to_cluster_id: return ".", "."
    if sequence_id not in read_id_to_cluster_id:
        print("Warning, sequence id "+str(sequence_id)+" not in clusters",file=sys.stderr)
        return ".", "."
    return read_id_to_cluster_id[sequence_id], cluster_id_to_cluster_size[read_id_to_cluster_id[sequence_id]]
    

def usage():
    '''Usage'''
    print("-----------------------------------------------------------------------------")
    print(sys.argv[0]," : discoSnp output filtering and formatting in vcf")
    print("-----------------------------------------------------------------------------")
    print("usage: ",sys.argv[0]," -i disco_bubbles.fa")
    print("  -r: min rank value filter (default = 0)")
    print("  -m: max missing value filter (default = 1)")
    print("  -o: output vcf file path (default = stdout)")
    print("  -f: output a filtered fasta file instead of a vcf file")
    print("  -c: considers a cluster input file. In this situation, prints the cluster_id in the output")
    print("  -h: help")
    print("-----------------------------------------------------------------------------")
    sys.exit(2)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:r:m:o:fc:", ["help", "in=", "rank=", "miss=", "out=", "fastaout"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    # Default parameters
    fasta_file = 0
    fasta_only = 0
    min_rank = 0
    max_miss = 1
    k = 31
    out_file =      None
    cluster_file =  None
    with_cluster = False
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--in"):
            fasta_file = arg
        elif opt in ("-r", "--rank"):
            min_rank = float(arg)
        elif opt in ("-m", "--miss"):
            max_miss = float(arg)
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-f", "--fastaout"):
            fasta_only = 1
        elif opt in ("-c"):
            cluster_file = arg
            with_cluster = True
        else:
            assert False, "unhandled option"

    if fasta_file == 0:
        print("option -i (--in) is mandatory")
        usage()
        sys.exit(2)
    else:
        today = time.localtime()
        date = str(today.tm_year) + str(today.tm_mon) + str(today.tm_mday)
        
        
        # First identifying what kind of fasta we have with the first line
        tig_type = 0 # 0 : no extension, 1 unitig, 2 contig (ie. unitig and contig length are output)
        nb_samples = 0
        nb_fixed_fields = 0 #nb fields before C1_X|C2_Y|... depends if unitig and/or contig lengths have been output
        
        ## LOAD clusters
        read_id_to_cluster_id, cluster_id_to_cluster_size = store_clusters(cluster_file)
        
        with open(fasta_file, 'r') as filin:
            for line in filin:
                splitted_1 = line.split("|")

                #nb_samples:
                tig_type = len(re.findall("left_\w+_length",line))
                nb_fixed_fields = 4 + 2*tig_type
                nb_samples = (len(splitted_1) - (nb_fixed_fields + 1))/3
                if nb_samples % 1 != 0:
                    print(f"Warning: could not detect the correct nb of samples : {nb_samples}")
                    sys.exit(2)
                nb_samples = int(nb_samples)
                break

        sys.stdout.close = lambda: None  #make stdout unclosable, to use with and handle both `with open(…)` and `sys.stdout` nicely. cf. https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
        # Now going through all lines
        with open(fasta_file, 'r') as filin, (open(out_file,'w') if out_file else sys.stdout) as filout:
            if not fasta_only:
                # Write vcf comment lines
                VCF_COMMENTS = f'''##fileformat=VCFv4.1
##filedate={date}
##source=create_filtered_vcf.py
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
                filout.write(VCF_COMMENTS + "\n")
                HEADER = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]+["G"+str(i+1) for i in range(nb_samples)])
                filout.write(HEADER + "\n")

            nb_kept_variants = 0
            nb_analyzed_variants = 0
            line_count = 0
            sequence_id=-1
            fasta_4lines = ""  ## Remembering 4 consecutive lines for kept variants to write in a fasta file if fasta_only mode
            keep_variant = False
            cluster_id=-1
            cluster_size=0
            for line in filin:
                if line_count%2 == 0:   sequence_id+=1  # first sequence is 1, second (lower path of first variant) is 2, ...
                if line_count ==0:      
                    fasta_4lines=""                     #back to empty for incoming variant 
                    cluster_id,cluster_size = get_cluster_id_and_size(sequence_id, read_id_to_cluster_id, cluster_id_to_cluster_size)
                line_count += 1
                fasta_4lines += line
                if line_count == 1:
                    keep_variant = False
                    nb_analyzed_variants += 1

                    # Header higher path
                    line = line.strip()
                    splitted_1 = line.split("|")
                    #fasta_4lines = splitted_1[0] + "\n"  #simplified headers for fasta_only and src

                    ## FILTERING
                    #filter rank
                    rank = float(splitted_1[-1].split("rank_")[1])
                    if rank < min_rank:
                        continue
                    # filter missing genotype ratio
                    if max_miss !=1:
                        nb_missing = len(re.findall(r"G\d+_\./\.",line))
                        missing_ratio = nb_missing / nb_samples
                        if missing_ratio >= max_miss:
                            continue
                    
                    keep_variant = True  # for fasta_only mode
                    nb_kept_variants += 1
                
                if keep_variant and line_count == 3:
                    #Header lower path
                    line = line.strip()
                    splitted_2 = line.split("|")
                    #fasta_4lines += splitted_2[0] + "\n"  #simplified headers for fasta_only and src

                if line_count == 4:
                    line_count = 0
                    if keep_variant:
                        if fasta_only:
                            #fasta_4lines += line   #simplified headers for fasta_only and src
                            filout.write(fasta_4lines) # TODO: do we writte fasta variants if not in a cluster and a cluster file is provided?
                        else:
                            #now format in vcf format
                            line = line.strip()
                            filout.write(format_vcf(splitted_1, splitted_2, nb_samples, nb_fixed_fields, rank, line, cluster_id, cluster_size, tig_type))
                    


        #print(f"{nb_lost_variants} variant bubbles filtered out")
        #print(f"{nb_kept_variants} variant bubbles output out of {nb_tot_variants} ({nb_analyzed_variants} analyzed)")
        sys.stderr.write(f"{nb_kept_variants} variant bubbles output out of {nb_analyzed_variants}\n")

def format_vcf(splitted_1, splitted_2, nb_samples, nb_fixed_fields, rank, sequence, cluster_id, cluster_size, tig_type = 0):
    '''
        Format information extracted from the fasta headers of a single bubble into one or several vcf lines
        
        tig_type :  0 no extension, 1 unitig, 2 contig (ie. unitig and contig length are present)
        sequence : DNA sequence of the lower path : usefull for INDELs (longest sequnce)
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
        
        POS = UL + POS - fuzziness (here = 3 = 3 possible positions for the 27 bp insertion) : 1-based, left-normalized
        '''
    
    # INFO = Ty=SNP;Rk=0.44073;UL=0;UR=0;CL=0;CR=0;Genome=.;Sd=.
    # FORMAT = GT:DP:PL:AD:HQ
    FORMAT = "GT:DP:PL:AD:HQ"
    
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
            POS, indel_size, fuzziness = re.findall("P_\d+:(\d+)_(\d+)_(\d+)",pol)[0]
            POS = int(POS) + position_offset - int(fuzziness)  # left-normalization, 1-based
            ALT = sequence[(POS-1):(POS+int(indel_size))]
            REF = ALT[0]
        ID = id
        if not isolated:
            ID += f"_{i}"
        my_line = "\t".join([CHROM, str(POS), ID, REF, ALT, ".", ".", INFO, FORMAT, GENO])
        vcf_line += my_line + "\n"
        i += 1

    return vcf_line


if __name__ == "__main__":
    main()
                      


