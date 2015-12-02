#!/usr/bin/env python
# -*- coding: utf-8 -*-
#V24102015 update : check the ref and the alt for the snp !
#Compare VCF to VCF or CSV to VCF
#Needs the number of chromosome in the studied genome
#You have to check if the format of the chromosome is the same in both files (chr1 ? or 1 ? or gi...)
#usage  python compare_vcf_V03112015.py <ref_file> <vcf_file> <boolean if chrom will be used or not (1 or 0)>

import sys
import getopt
import re
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
# homozygous = 1
# hetero = 0
def transform_geno(geno):
    sep="/"
    if geno[1]=="|": sep="|"
    if geno.split(sep)[0] == geno.split(sep)[1]: return 1
    return 0
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
def index_reference_csv(ref_vcf_file):
    index={}
    nb_SNP=0
    nb_indel=0
    filin = open(ref_vcf_file, 'r')
    typeFile=0
    while True:
        indel=False
        line = filin.readline().rstrip()
        if not line: break
        if line.startswith("run") : 
                typeFile=1
                #csv
                #   run;patient;function;gene;exonicFunction;chromosome;exonAnnotation;start;stop;size;reference;alteration;NM;nucleotideChange;aminoAcideChange;dpTot;dp_a1-R1;dp_a1-R2;dp_a2-R1;dp_a2-R2;dpTotAlt;dpAlt_a1-R1;dpAlt_a1-R2;dpAlt_a2-R1;dpAlt_a2-R2;faAltMean;faAlt_a1-R1;faAlt_a1-R2;faAlt_a2-R1;faAlt_a2-R2;pqAltMean;pqAlt_a1-R1;pqAlt_a1-R2;pqAlt_a2-R1;pqAlt_a2-R2;amplicon;exonAln;pcrName
                continue
        if line.startswith("PC") :
               #csv
               #PC;GENE;CHR;POS;EXON;AA;NT;AF;PCR_NAME
               typeFile=2
               continue  
        #if line.split()[6] != "PASS": continue
        # only position for now:
        # print line
        
        
        if typeFile==1:
                chro ="chr"+str(line.split(";")[5])
                pos = int(line.split(";")[7])
                
                ref = line.split(";")[10]
                #if ',' in ref : continue
                alt = line.split(";")[11]
                #if ',' in alt : continue
        if typeFile==2:
                #PC;GENE;CHR;POS;EXON;AA;NT;AF;PCR_NAME
                #c.1799T>A
                #ref = line.split(";")[6]
                #listref=re.findall('(\d+|[A-Za-z]|>|\.)',ref)
                ref=line.split(";")[4]
                alt=line.split(";")[5]
                print line.split(";")
                chro ="chr"+str(line.split(";")[2])
                pos = int(line.split(";")[3])
                #parsingCigarCode=re.findall('(\d+|[A-Za-z])',cigarcode)
                # if ">" not in listref :
#                       indel==True
#                       ref="AA"
#                       alt="A"
#                 else:
#                         #ref=listref[3]
#                 #if ',' in ref : continue
#                         #alt=listref[5]
#                         ref="A"
#                         alt
#                 #if ',' in alt : continue
        if not chro in index: index[chro]={}
        if not pos in index[chro]: 
            index[chro][pos]=[] 
            if len(ref)==len(alt): 
                nb_SNP+=1
                thistype="SNP"
            else:
                nb_indel+=1
                thistype="INDEL"
            index[chro][pos].append(thistype)
            thistype=""
            indel=False
            index[chro][pos].append(False)
            index[chro][pos].append(ref)
            index[chro][pos].append(alt)
            index[chro][pos].append(ReverseComplement(ref))
            index[chro][pos].append(ReverseComplement(alt))
            
    filin.close()
    return index,nb_SNP,nb_indel

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
def ReverseComplement(nucleotide):
	if nucleotide=="A": return "T"
	if nucleotide=="T": return "A"
	if nucleotide=="C": return "G"
	return "C"

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
def index_reference_vcf(ref_vcf_file,boolChrom):
    index={}
    nb_SNP=0
    nb_indel=0
    filin = open(ref_vcf_file, 'r')
    thistype=""
    while True:
        # VCF 1000 humans: 
        # ['1', '46588503', 'rs10890383', 'A', 'G', '100', 'PASS', 'ERATE=0.0004;AN=2184;AC=877;THETA=0.0005;VT=SNP;AA=A;LDAF=0.4026;SNPSOURCE=LOWCOV;RSQ=0.9841;AVGPOST=0.9885;AF=0.40;ASN_AF=0.59;AMR_AF=0.47;AFR_AF=0.07;EUR_AF=0.44', 'GT:DS:GL'] 0|1:1.000:-4.70,-0.00,-3.47 0|0:0.000:-0.10,-0.70,-4.10
        line = filin.readline()
        if not line: break
        if line.startswith("#") : continue
        if line.split("\t")[6] != "PASS": continue
        # only position for now:
        # print line
        if boolChrom==0:
                pos = int(line.split("\t")[1])
                ref = line.split("\t")[3]
                alt = line.split("\t")[4]
                if not pos in index:
                        index[pos]=[] 
                        if len(ref)==len(alt):
                                nb_SNP+=1 
                                thistype="SNP"
                        else: 
                                nb_indel+=1
                                thistype="INDEL"
                        index[pos].append(thistype)
                        index[pos].append(False)
        
        if boolChrom==1:
                chro ="chr"+str(line.split("\t")[0])
                pos = int(line.split("\t")[1])
                #re.findall('(\d+|[A-Za-z])',cigarcode)
                ref = line.split("\t")[3]
                #if ',' in ref : continue
                alt = line.split("\t")[3]
                #if ',' in alt : continue
                
                if not chro in index: index[chro]={}
                if not pos in index[chro]: 
                    index[chro][pos]=[] 
                    if len(ref)==len(alt): 
                        nb_SNP+=1
                        thistype="SNP"
                    else: 
                        nb_indel+=1
                        thistype="INDEL"
                    index[chro][pos].append(thistype)
                    index[chro][pos].append(False)                     
    filin.close()
    return index,nb_SNP,nb_indel
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------   
def get_index_entry (chro,pos, index,wantedtype,nb_chro,ref,alt):
    span=20
    if len(index)<=nb_chro: 
        if chro not in index : return index,False
        if wantedtype=="SNP":  
            if pos in index[chro] and index[chro][pos][0]=="SNP":
                    #if (index[chro][pos][2]==ref or index[chro][pos][4]==ref) and (index[chro][pos][3]==alt or index[chro][pos][5]==alt):
                            index[chro][pos][1]=True
                            return (index,True)
            #if pos-1 in index[chro] and index[chro][pos-1][0]=="SNP" :
                    #if (index[chro][pos-1][2]==ref or index[chro][pos-1][4]==ref) and (index[chro][pos-1][3]==alt or index[chro][pos-1][5]==alt):
             #               index[chro][pos-1][1]=True
             #               return (index,True)
            #if pos+1 in index[chro] and index[chro][pos+1][0]=="SNP" :
                    #if (index[chro][pos+1][2]==ref or index[chro][pos+1][4]==ref) and (index[chro][pos+1][3]==alt or index[chro][pos+1][5]==alt):
             #               index[chro][pos+1][1]=True
             #               return (index,True)
            else: return index,False      
        # this is an indel we check for previous and next entries:
        elif wantedtype=="INDEL":            
                for i in range(pos,pos-span,-1):
                    if i in index[chro] and index[chro][i][0]=="INDEL":
                            index[chro][i][1]=True
                            return (index,True)
                    
                for i in range(pos+1,pos+span):
                    if i in index[chro] and index[chro][i][0]=="INDEL":
                            index[chro][i][1]=True
                            return (index,True)  
        return index,False            
            
            
            
    else:
        if wantedtype=="SNP":
            if pos in index and index[pos][0]=="SNP":
                    index[pos][1]=True
                    return (index,True)
            elif pos-1 in index and index[pos-1][0]=="SNP":
                    index[pos-1][1]=True
                    return (index,True)
            elif pos+1 in index and index[pos+1][0]=="SNP":
                    index[pos+1][1]=True
                    return (index,True)
                  
        # this is an indel we check for previous and next entries:
        elif wantedtype=="INDEL":
                for i in range(pos,pos-span,-1):
                    if i in index and index[i][0]=="INDEL":
                            index[i][1]=True
                            return (index,True)
                    
                for i in range(pos+1,pos+span):
                    if i in index and index[i][0]=="INDEL":
                            index[i][1]=True
                            return (index,True)  
        return index,False

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------    
def countDict(dico,thistype,nb_chro):
        compt=0
        if len(dico)>=nb_chro:        
                for variantType,val in dico.values():
                        if val==True and variantType==thistype:
                                compt+=1
        if len(dico)<=nb_chro:
                for sousdico in dico.values():
                        for key,(variantType,val,ref,alt,refR,altR) in sousdico.items():
                                if val==True and variantType==thistype:
                                        compt+=1
        return compt
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------    
def analyse_predicted_vcf(prediction_vcf_file,index,wantedtype,nb_ref, roc_file_name,nb_chro):

    nb_seen=0
    nb_TP=0
    nb_FP=0
    nb_seenPass=0
    nbIndel=0
    nbSnp=0
    nb_genotype_TP=0
    nb_TP_seen=0
    
    nb_PASS=0
    nb_TP_PASS=0
    nb_FP_PASS=0
    nb_genotype_TP_PASS=0
    filin = open(prediction_vcf_file, 'r')
    
    while True:
        #gi|224384768|gb|CM000663.1|	10583	.	G	A	401.19	.	AC=1;AF=0.250;AN=4;BaseQRankSum=-0.557;ClippingRankSum=-1.021;DP=48;FS=1.248;MLEAC=1;MLEAF=0.250;MQ=60.00;MQ0=0;MQRankSum=1.508;QD=20.06;ReadPosRankSum=0.974;SOR=0.474	GT:AD:DP:GQ:PL	0/0:28,0:28:84:0,84,1058	0/1:7,13:20:99:431,0,204
        line = filin.readline() 
        if not line: break
        if line.startswith("#") or line.startswith("0"): continue
        listLine=line.split("\t")
        pos=int(listLine[1])
        ref=listLine[3]
        alt=listLine[4]
        chro=listLine[0]
        if len(ref)==len(alt): 
                thistype="SNP"
                nbSnp+=1
        else:
                thistype="INDEL"
                nbIndel+=1        
        if wantedtype!=thistype: continue
        # only position for now: 
        index,index_line = get_index_entry(chro,pos,index,wantedtype,nb_chro,ref,alt)
        if index_line:
                print line.rstrip()
                nb_TP+=1
        else:
                nb_FP+=1
    
    nb_TP_seen=countDict(index,wantedtype,nb_chro) 
    print "#############################################"
    print "Results for", wantedtype
    print "Total simulated", nb_ref
    print "Total predicted", nb_TP+nb_FP
    print "Total TP", nb_TP
    print "Precision =", 100*nb_TP/float(nb_FP+nb_TP)
    print "recall =", 100*nb_TP_seen/float(nb_ref)
    print "Number of INDEL",nbIndel
    print "Number of SNP", nbSnp
    #print "good Genotypes (of TP)", nb_genotype_TP
    
    
    #print "good Genotypes (of TP & PASS)", nb_genotype_TP_PASS   
index=None
nb_chro=46
boolChrom=sys.argv[3]
print "INDEXING"
if ".vcf" in sys.argv[1]:
    index,nb_SNP,nb_indel=index_reference_vcf(sys.argv[1],boolChrom)
elif ".csv" in sys.argv[1]:
    index,nb_SNP,nb_indel=index_reference_csv(sys.argv[1])
print "ANALYSING"
fileName=sys.argv[2]
if ".vcf" in fileName:
        analyse_predicted_vcf(sys.argv[2],index,"SNP", nb_SNP, "roc_SNP_from_vcf",nb_chro)
        analyse_predicted_vcf(sys.argv[2],index,"INDEL", nb_indel, "roc_INDEL_from_vcf",nb_chro)

else :
        print "!!! WRONG FILE FORMAT !!!"
        
        
        
            
