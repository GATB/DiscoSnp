import sys

filin = open(sys.argv[1], 'r')


while True:    
    line1=filin.readline()[:-1]
    if not line1: break
    line2=filin.readline()[:-1]
    line3=filin.readline()[:-1]
    line4=filin.readline()[:-1]
    
    print float(line1.split("|")[1]),line1+";"+line2+";"+line3+";"+line4