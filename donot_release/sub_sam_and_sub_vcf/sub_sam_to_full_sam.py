#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt


# Given a SAM formated as 
# INDEL_lower_path_390|P_1:30_3_9|high|nb_pol_1|C1_2148|C2_1714|Q1_69|Q2_68|G1_1/1:41976,5624,67|G2_1/1:33373,4389,86|rank_0.01396        16      chr11:530242..537550    4033    37      54M     *       0       0       TCTTGCCCACACCGCCGGCGCCCACCACCACCAGCTTATATTCCGTCATCGCTC  *       XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:54
# Transform it as:
# INDEL_lower_path_390|P_1:30_3_9|high|nb_pol_1|C1_2148|C2_1714|Q1_69|Q2_68|G1_1/1:41976,5624,67|G2_1/1:33373,4389,86|rank_0.01396        16      chr11    534285    37      54M     *       0       0       TCTTGCCCACACCGCCGGCGCCCACCACCACCAGCTTATATTCCGTCATCGCTC  *       XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:54
# that is to say: put back the right position: 4043+530242 (beginning of the portion of mapped sequence, starting position 530242 on the chromosome 11


#A difficulity stands in the fact that for each line, maybe it does not contain any position information. In this case we don't modify the line


filein=open(sys.argv[1],'r')

for line in filein: 
    if line[0]=='@':
        print line,
        continue
    values=line.split()

    if not ':' in values[2]: # no position information, we continue
        print line, 
        continue
        
    start=int(values[2].split(':')[1].split('.')[0])
    position=int(values[3])+start-1
    for i in xrange(2):
        print values[i]+"\t",
    print values[2].split(':')[0]+"\t"+str(position)+"\t",
    for i in xrange(4,len(values)):
        print values[i]+"\t",
    print