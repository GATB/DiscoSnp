#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt


# Given a VCF formated as 
# chr11:530242..537550    4043    390     ACCA    A       .       PASS    Ty=INS;Rk=0.01396;UL=.;UR=.;CL=.;CR=.;Genome=.;Sd=-1    GT:DP:PL:AD     0/0:2196:67,5624,41976:2148,48  0/0:1760:86,4389,33373:1714,46
# Transform it as:
# chr11    534285    390     ACCA    A       .       PASS    Ty=INS;Rk=0.01396;UL=.;UR=.;CL=.;CR=.;Genome=.;Sd=-1    GT:DP:PL:AD     0/0:2196:67,5624,41976:2148,48  0/0:1760:86,4389,33373:1714,46
# that is to say: put back the right position: 4043+530242 (beginning of the portion of mapped sequence, starting position 530242 on the chromosome 11


filein=open(sys.argv[1],'r')

for line in filein:
    if line[0]=='#':
        print line,
        continue
    values=line.split()
    start=int(values[0].split(':')[1].split('.')[0])
    position=int(values[1])+start-1
    print values[0].split(':')[0]+"\t"+str(position)+"\t",
    for i in xrange(2,len(values)):
        print values[i]+"\t",
    print