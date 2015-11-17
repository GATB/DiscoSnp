#!/bin/python
# -*- coding: utf-8 -*-
###################################
#
#
import sys
file=open(sys.argv[1],"r")
minimal_coverage=int(sys.argv[2])


for line in file:
    if "#" == line[0]: 
        print line, 
        continue
    line=line.rstrip('\n')
    #SNP_higher_path_2201025 45      2201025_1       C       T       .       .       Ty=SNP;Rk=0.99973;UL=14;UR=1;CL=14;CR=1;Genome=.;Sd=.   GT:DP:PL:AD     0|1:0:4,4,4:0,0 1|1:159:3184,483,11:0,159
    splitted_line = line.split()
    for i in range(9): print splitted_line[i],
    # decrease the position value if it is a number:
    for i in range(9,len(splitted_line)):
        covs = splitted_line[i].split(":")[-1].split(",")
        covA=covs[0]
        covB=covs[1]
        
        if int(covA)>minimal_coverage or int(covB)>minimal_coverage: print splitted_line[i],
        else: print "-9"+splitted_line[i][3:],
    print

