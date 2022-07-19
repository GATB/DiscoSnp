#!/usr/bin/env python3
import sys
if len(sys.argv) !=2:
    sys.stdout.write("Mandatory: python discoSnp_to_csv.py prefix_coherent_k_kval_c_cval.fa\n")
    sys.stdout.write("This program formats the .fa to .csv format by puting each couple of .fa sequence (4 lines = 2 comments + 2 nucleotide sequences) into one line, replacing the '|' character by spaces and removing the CX_ formating")
    sys.exit(1)

f=open(sys.argv[1], "r")


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
    
    com1_tab=com1_1.split("|")
    
    # prints all before coverages
    for i in range(0,4):
        sys.stdout.write( com1_tab[i]+",")
        
    # prints coverages
    i=4
    while com1_tab[i][0:1]!="C" and i<len(com1_tab): 
        i+=1
    while com1_tab[i][0:1]=="C":
        sys.stdout.write( com1_tab[i].split("_")[1]+",")
        i+=1

    # prints all remaining
    while i<len(com1_tab)-1:
        sys.stdout.write( com1_tab[i]+",")
        i+=1
    sys.stdout.write( com1_tab[i][:-1]+",")
    
    sys.stdout.write(data1_1[:-1]+",")
    
    com2_tab=com1_2.split("|")
    
    # prints all before coverages
    for i in range(0,4):
        sys.stdout.write( com2_tab[i]+",")
        
    # prints coverages
    i=4
    while com1_tab[i][0:1]!="C" and i<len(com1_tab): 
        i+=1
    while com2_tab[i][0:1]=="C":
        sys.stdout.write( com2_tab[i].split("_")[1]+",")
        i+=1

    # prints all remaining
    while i<len(com2_tab)-1:
        sys.stdout.write( com2_tab[i]+",")
        i+=1
    sys.stdout.write( com2_tab[i][:-1]+",")
    
    print (data1_2,)
    
    
    




