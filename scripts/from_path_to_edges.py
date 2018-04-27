import sys
if len(sys.argv)>1:
    file = open(sys.argv[1])
else:
    file=sys.stdin
# print("source\t target\t type")
def printLine(line, type='r'):
    for i in range(len(line)-1):
        if int(line[i][:-1])<int(line[i+1][:-1]):
            print (line[i]+'\t'+line[i+1]+"\t"+type)
        else:
            print (line[i+1]+'\t'+line[i]+"\t"+type)


# 1000547h;2286435h; 1330792h;1152525l;
for line in file: 
    line=line.strip().rstrip().split(' ')
    line1=line[0].rstrip().split(';')[:-1]
    printLine(line1)
    if len(line)==2: # pairend links:
        line2=line[1].rstrip().split(';')[:-1]
        printLine(line2)
        # print the pairend link:
        printLine([line1[-1],line2[0]],'p')
        
    
    

    
    
    
