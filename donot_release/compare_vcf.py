#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt

# homozygous = 1
# hetero = 0
def transform_geno(geno):
    sep="/"
    if geno[1]=="|": sep="|"
    if geno.split(sep)[0] == geno.split(sep)[1]: return 1
    return 0

def index_reference_vcf(ref_vcf_file):
    index={}
    nb_SNP=0
    nb_indel=0
    filin = open(ref_vcf_file, 'r')
    
    while True:
        # VCF 1000 humans: 
        # ['1', '46588503', 'rs10890383', 'A', 'G', '100', 'PASS', 'ERATE=0.0004;AN=2184;AC=877;THETA=0.0005;VT=SNP;AA=A;LDAF=0.4026;SNPSOURCE=LOWCOV;RSQ=0.9841;AVGPOST=0.9885;AF=0.40;ASN_AF=0.59;AMR_AF=0.47;AFR_AF=0.07;EUR_AF=0.44', 'GT:DS:GL'] 0|1:1.000:-4.70,-0.00,-3.47 0|0:0.000:-0.10,-0.70,-4.10
        line = filin.readline() 
        if not line: break
        if line.startswith("#") : continue
        # only position for now:
        # print line
        
        
        
        pos = int(line.split("'")[3])
        ref = line.split("'")[7]
        alt = line.split("'")[9]
        if len(ref)==len(alt): nb_SNP+=1
        else: nb_indel+=1
        geno1 = line.split("]")[1].split(" ")[1].split(":")[0].strip()
        geno2 = line.split("]")[1].split(" ")[2].split(":")[0].strip()
        
        
        index[pos]=[]
        index[pos].append(transform_geno(geno1))
        index[pos].append(transform_geno(geno2))
    filin.close()
    return index,nb_SNP,nb_indel
    

def get_index_entry (pos, index, wantedtype):
    
    if wantedtype=="SNP":
        if pos-1 in index: return index[pos-1]
        if pos in index: return index[pos]
        if pos+1 in index: return index[pos+1]
        return None
    # this is an indel we check for previous and next entries:
    span=20
    for  i in range(pos,pos-span,-1):
        if i in index: return index[i]
        
    for i in range(pos+1,pos+span):
        if i in index: return index[i]
    
    return None
        
    
    
    
def analyse_predicted_vcf(prediction_vcf_file,index,wantedtype,nb_ref, roc_file_name):

    nb_seen=0
    nb_TP=0
    nb_FP=0
    
    nb_genotype_TP=0
    
    nb_PASS=0
    nb_TP_PASS=0
    nb_FP_PASS=0
    nb_genotype_TP_PASS=0

    prev_out = sys.stdout
    f = file(roc_file_name, 'w')
    sys.stdout = f
    
    filin = open(prediction_vcf_file, 'r')
    
    while True:
        #gi|224384768|gb|CM000663.1|     108430814       99995   T       C       .       PASS    Ty=SNP;Rk=1.00000;MULTI=.;DT=0;UL=.;UR=.;CL=.;CR=.;C1=42,0;C2=0,38;Genome=T;Sd=1        GT:DP:PL        0/0:42:0        1/1:38:0
        line = filin.readline() 
        if not line: break
        if line.startswith("#") or line.startswith("0"): continue
        
        thistype=line.split("\t")[7].split(";")[0].split("=")[1].strip()
        if thistype=="INS" or thistype=="DEL": thistype="INDEL"
        if wantedtype!=thistype: continue
        
        nb_seen+=1
        # only position for now:
        pos = int(line.split("\t")[1].strip())-1 # -1 as our vcf is 1 based
        PASS=True
        if not line.split("\t")[6].strip()=="PASS": PASS=False
        rank = float(line.split("\t")[7].split(";")[1].split("=")[1].strip())
        geno1=transform_geno(line.split("\t")[9].split(":")[0].strip())
        geno2=transform_geno(line.split("\t")[10].split(":")[0].strip())
        # print pos,PASS,rank,geno1,geno2
        
        
        index_line = get_index_entry(pos,index,wantedtype)
        if index_line:
            nb_TP+=1
            if geno1==index_line[0]: nb_genotype_TP+=1
            if geno2==index_line[1]: nb_genotype_TP+=1
            
        else:
            nb_FP+=1
            
        print 100*nb_TP/float(nb_ref), 100*nb_TP/float(nb_seen), rank
            
        if PASS:
            nb_PASS+=1
            if index_line:
                nb_TP_PASS+=1

                if geno1==index_line[0]: nb_genotype_TP_PASS+=1
                if geno2==index_line[1]: nb_genotype_TP_PASS+=1
            else:
                nb_FP_PASS+=1
                

    f.close()
    sys.stdout = prev_out      
    print "#############################################"
    print "Results for", wantedtype
    print "Total simulated", nb_ref
    print "Total predicted", nb_TP+nb_FP
    print "Total TP", nb_TP
    print "Precision =", 100*nb_TP/float(nb_FP+nb_TP)
    print "recall =", 100*nb_TP/float(nb_ref)
    
    
    print "good Genotypes (of TP)", nb_genotype_TP
    
    print "Among mapped predictions (PASS in the predicted VCF):"
    print "Total TP PASS", nb_TP_PASS
    print "Total FP PASS", nb_FP_PASS
    
    
    print "good Genotypes (of TP & PASS)", nb_genotype_TP_PASS
    

index=None
print "INDEXING"
index,nb_SNP,nb_indel=index_reference_vcf(sys.argv[1])
print "ANALYSING"
analyse_predicted_vcf(sys.argv[2],index,"SNP", nb_SNP, "roc_SNP_from_vcf")
analyse_predicted_vcf(sys.argv[2],index,"INDEL", nb_indel, "roc_INDEL_from_vcf")