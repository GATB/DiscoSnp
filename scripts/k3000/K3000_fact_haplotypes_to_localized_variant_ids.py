#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Given a graph (gfa format) in which nodes are facts and given some paths in this graph:
Output for each path the corresponding facts with their relative order and inner distances. 

Fragment of the gfa file: 
S       2       aCTTTACTTGCTTTACATAGAAGTTATTTGACTCCTGGTGATtcttcttcaGGTTGGACAgCTGGTGCTGCAgCTTATTATGTggNNNNNNNNNNNNctAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACagatgctg    AS:-77l;513l;118l;-1059h;1001h; SP:0_51;12_61;22_73;41_85;97_1
48;       BP:0_41;-23_41;-29_41;-30_41;16_41;     FC:i:270        min:319 max:491 mean:377.4      AC:364;319;335;378;491;
S       211     gcacaGTCTACAGCATCTGTAATGGTTCCATTTTCATTATATTTTAatagAaaagtcCTagNNNNNNNNNNNNccACATAATAAGCTGCAGCACCAGCTGTCCAACCTGAAGAAGAaTCACCAGGAGtcaaataacTTCtatgtaa      AS:768l;310h;-1001h;1059h;-118l;-513l;  SP:0_50;7_57;10_61;73_117;85_136;97_146;        BP:0_41;-36_41;-33_41;16_41;-30_41;-29_41;      FC:i:532        min:319 max:555 mean:431.1666666666667  AC:555;509;491;378;335;319;
L       2       +       211     -       123M

'''

__author__ = "Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr"
