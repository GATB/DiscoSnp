#!/usr/bin/env python
import sys


def get_central_nuc(s1,s2):
    s_maj=""
    for l in s1:
        if l.isupper():
            s_maj=s_maj+l
    central1=s_maj[len(s_maj)/2]
    s_maj=""
    for l in s2:
        if l.isupper():
            s_maj=s_maj+l
    central2=s_maj[len(s_maj)/2]
    
    return s_maj[0:len(s_maj)/2]+central1+"/"+central2+s_maj[(len(s_maj)/2)+1:]
    
    


if len(sys.argv) !=3:
    print "Mandatory: python discoSnp_to_genotypes.py prefix_coherent_k_kval_c_cval.fa threshold_value"
    print "This program extracts the \genotypes\" of each SNP:"
    print "\t for each SNP and each input data set it indicates if the SNP is:"
    print "\t\t -homozygous ALT1 path (coverage ALT1 >= threshold and ALT2 < threshold): 1"
    print "\t\t -homozygous ALT2 path (coverage ALT1 < threshold and ALT2 >= threshold): -1"
    print "\t\t -heterozygous (coverage ALT1 >= threshold and ALT2 >= threshold): 2"
    print "\t\t -absent (coverage ALT1 < threshold and ALT2 < threshold): 0"
    print "\t then it outputs the central sequence of length 2k-1 replacing the central position by ALT1/ALT2"
    sys.exit(1)

f=open(sys.argv[1], "r")

t=int(sys.argv[2])


while 1:
    com1_1=f.readline()
    if not com1_1:
        break
    data1_1=f.readline()
    if not data1_1:
        break
    com1_2=f.readline()
    if not com1_2:
        break
    data1_2=f.readline()
    if not data1_2:
        break
    

    sequence=get_central_nuc(data1_1, data1_2)
    com1_tab=com1_1.split("|")
    com2_tab=com1_2.split("|")
    
    
    
    # get coverages for computing genotypings
    sys.stdout.write("GENOTYPES_SNP_")
    print com1_tab[0].split("_")[-1],
    sys.stdout.write("_THRESHOLD_")
    print t,
    i=4
    while com1_tab[i][0:1]!="C" and i<len(com1_tab): 
        i+=1
    while com1_tab[i][0:1]=="C":
        up=int(com1_tab[i].split("_")[1])
        lo=int(com2_tab[i].split("_")[1])
        if up>=t and lo<t: print 1,
        if up<t and lo>=t: print -1,
        if up<t and lo<t: print 0,
        if up>=t and lo>=t: print 2,
        i+=1
    print sequence
    




