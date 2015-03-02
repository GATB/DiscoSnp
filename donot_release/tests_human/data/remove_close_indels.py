import sys

filin = open(sys.argv[1], 'r')

k=31
#SNP 736288 gi T A
#INDEL 739140 gi AT A
#SNP 739209 gi A G

line1=filin.readline()
line2=filin.readline()
line3=filin.readline()
print line1,
while True:
    pos1=int(line1.split(" ")[1])
    pos2=int(line2.split(" ")[1])
    pos3=int(line3.split(" ")[1])
    
    if line2.startswith("SNP"):
        print line2,
    if line2.startswith("INDEL") and pos2-pos1>k and pos3-pos1>k:
        print line2,
    
    line1=line2
    line2=line3
    line3=filin.readline()
    if not line3: break

print line2,