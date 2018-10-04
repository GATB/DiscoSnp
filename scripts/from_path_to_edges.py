import sys
if len(sys.argv)>1:
    file = open(sys.argv[1])
else:
    file=sys.stdin
# print("source\t target\t type")
def printLine(line, type='r'):
    for i in range(len(line)-1):
        if int(line[i].split('_')[0][:-1])<int(line[i+1].split('_')[0][:-1]):
            print (line[i].split('_')[0]+'\t'+line[i+1].split('_')[0]+"\t"+type)
        else:
            print (line[i+1].split('_')[0]+'\t'+line[i].split('_')[0]+"\t"+type)


# 1000547h_0;2286435h_12; 1330792h_0;1152525l_24;
for line in file: 
    line=line.strip().rstrip().split(' ')
    line1=line[0].rstrip().split(';')[:-1]
    printLine(line1)
    if len(line)==2: # pairend links:
        line2=line[1].rstrip().split(';')[:-1]
        printLine(line2)
        # print the pairend link:
        printLine([line1[-1],line2[0]],'p')
        
    
    

    
    
    
