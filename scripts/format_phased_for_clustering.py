import sys

if len(sys.argv)>1:
    file = open(sys.argv[1])
else:
    file=sys.stdin
    
for line in file: 
    line=line.strip().rstrip().split()
    print (line[0]+":",end='')
    for i in line:
        print(i+" ",end='')
    print ()

    
    
    