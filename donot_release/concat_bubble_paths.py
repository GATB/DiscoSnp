#!/usr/bin/env python
import sys


f=open(sys.argv[1], "r")

#USE FOR ASSESSING THE UNICITY OF EACH PREDICTED BUBBLE.
# python concat_bubble_paths.py foo_coherent.fa > trash_me.txt
# wc -l trash_me.txt
# sort -u trash_me.txt | wc -l
# If results differ, the bubble detection is redundant.


while 1:
    com1_1=f.readline().rstrip()
    if not com1_1:
        break
    data1_1=f.readline().rstrip()
    if not data1_1:
        break
    com1_2=f.readline().rstrip()
    if not com1_2:
        break
    data1_2=f.readline().rstrip()
    if not data1_2:
        break
    

    print data1_1+data1_2
    




