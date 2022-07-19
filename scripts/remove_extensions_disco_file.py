#!/bin/python3
# -*- coding: utf-8 -*-
###################################
#Removes extensions in lowercase :
# usage : remove_extensions_disco_file.py <disco_snps_file>.fa
import os
import sys
import getopt
import time


fichier=sys.argv[1]
output=sys.argv[2]
nomFichier=(str(fichier).split('.'))
##Delete extension (lowercase) from disco file
snpfile=open(fichier,"r")
snpDiscobis=open(output, "w")
for line in snpfile:
	if ">"==line[0]:
		snpDiscobis.write(line)
	else:
		line=str(line)
		line=line.replace('a','')
		line=line.replace('g','')
		line=line.replace('t','')
		line=line.replace('c','')
		snpDiscobis.write(line)
snpDiscobis.close()

