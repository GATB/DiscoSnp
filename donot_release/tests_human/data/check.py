import sys

filin = open(sys.argv[1], 'r')
letters={}

tot=0
for line in filin:
    
    if line[0] == ">": continue
    for i in range(len(line)-1):
        tot+=1
        if tot%1000000==0:
            print letters
        if line[i] not in letters:
            letters[line[i]]=1
        else:
            letters[line[i]]+=1

    
print letters