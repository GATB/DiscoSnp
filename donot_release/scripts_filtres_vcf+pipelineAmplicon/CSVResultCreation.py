#!/usr/bin/env python
# -*- coding: utf-8 -*-

#to Run it on every log file:
#for log in $( ls |grep "log"); do python ../Script/CSVResultCreation.py $log ; done

#It will create : 
#ResultatsAmpliconSNP.csv
#ResultatsAmpliconINDEL.csv



#log parsing : 
#Found and parse in log file of pipeline+MAF.sh or log file of compare:
#INDEXING
#ANALYSING
#############################################
#Results for SNP
#Total simulated 2688
#Total predicted 1496
#Total TP 932
#Precision = 62.2994652406
#recall = 23.6235119048
#Number of INDEL 67
#Number of SNP 1496
#############################################
#Results for INDEL
#Total simulated 4
#Total predicted 67
#Total TP 8
#Precision = 11.9402985075
#recall = 100.0
#Number of INDEL 67
#Number of SNP 1496


### In a csv file
#Name of File;Total simulated;Total predicted;Total TP;Precision;recall
#logTPCnb5aWithRef002.txt;2688;1496;932;932;62.2994652406;23.6235119048
#logTPCnb5aWithRef002.txt;2688;1496;932;932;62.2994652406;23.6235119048
#logTPCnb5aWithRef002.txt;2688;1496;932;932;62.2994652406;23.6235119048
#logTPCnb5aWithRef002.txt;2688;1496;932;932;62.2994652406;23.6235119048
#logTPCnb5aWithRef002.txt;2688;1496;932;932;62.2994652406;23.6235119048





import sys
import re
import os

fileName=sys.argv[1]
outputCSVSNP=str("ResultatsAmpliconSNP.csv")
outputCSVINDEL=str("ResultatsAmpliconINDEL.csv")
logFile=open(fileName,'r')


if os.path.isfile(outputCSVSNP):
        outputCSVSNP = open(outputCSVSNP,'a')
        
else:
        outputCSVSNP = open(outputCSVSNP,'a')
        outputCSVSNP.write("Name of File"+str(";")+"Total simulated"+str(";")+"Total predicted"+str(";")+"Total TP"+str(";")+"Precision"+str(";")+"recall\n")

if os.path.isfile(outputCSVINDEL):
        outputCSVINDEL = open(outputCSVINDEL,'a') 
else:
        outputCSVINDEL = open(outputCSVINDEL,'a')
        outputCSVINDEL.write("Name of File"+str(";")+"Total simulated"+str(";")+"Total predicted"+str(";")+"Total TP"+str(";")+"Precision"+str(";")+"recall\n")

while True:
        line = logFile.readline() 
        if not line: break
        if "INDEXING" in line :
                line = logFile.readline()
                if "ANALYSING" in line:
                          line=logFile.readline()
                          line=logFile.readline()
                          total_simulated_snp=logFile.readline().split()[2]
                          total_predicted_snp=logFile.readline().split()[2]
                          total_tp_snp=logFile.readline().split()[2]
                          precision_snp=logFile.readline().split()[2]
                          recall_snp=logFile.readline().split()[2]
                          outputCSVSNP.write(str(fileName)+str(";")+str(total_simulated_snp)+str(";")+str(total_predicted_snp)+str(";")+str(total_tp_snp)+str(";")+str(precision_snp)+str(";")+str(recall_snp)+str("\n"))
                          logFile.readline()
                          logFile.readline()
                          logFile.readline()
                          logFile.readline()
                          total_simulated_indel=logFile.readline().split()[2]
                          total_predicted_indel=logFile.readline().split()[2]
                          total_tp_indel=logFile.readline().split()[2]
                          precision_indel=logFile.readline().split()[2]
                          recall_indel=logFile.readline().split()[2]
                          outputCSVINDEL.write(str(fileName)+str(";")+str(total_simulated_indel)+str(";")+str(total_predicted_indel)+str(";")+str(total_tp_indel)+str(";")+str(precision_indel)+str(";")+str(recall_indel)+str("\n"))
                          
outputCSVINDEL.close()
outputCSVSNP.close()
logFile.close()
                          
                          
                          
                          
