#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Index Functions with the python dictionnary
'''
import string

def indexGenome(genome,k):
    index={}
    for i in range(0,len(genome)-k+1):
        key=genome[i:i+k]
        if key in index:
            index[key].append(i)
        else:
            index[key]=[i]
    return index

def queryIndex(index,seed):
    if str(seed) in index:
        return index[seed]
    else:
        return -1
        
