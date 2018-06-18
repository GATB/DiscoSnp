# cf https://www.evernote.com/l/ARUGTnFwaO9Kf7CO5fYJhYXh4T20R54V6eU

'''
Transforms the phased alleles from kissreads to a K2000 input format: (https://github.com/Malfoy/BWISE/blob/master/src/K2000/K2000.py)
#comments blabla toto
-1064h;-917l;1880l; => 2
...

becomes
#comments blabla toto
-2128;-1835;3761;


Moreover, K2000 is not able to deal with long range interactions (non overlaping). Thus pairend variants are forgotten:
-1064h;-917l;1880l; 24h;1l => 2

becomes: 
-2128;-1835;3761;
48;3;
'''

import sys

if len(sys.argv)>1:
    file = open(sys.argv[1])
else:
    file=sys.stdin
    
def f(variant):
    ''' 
    sVp 
    * s='-' or nothing
    * V=int value
    * p='h' or 'l' (path)

    f(sVp)=s2V+g(p) with g(p)=0 if p='h' and 1 if p='l'
    '''
    s=''
    if variant[0]=='-': 
        s='-'
        V=int(variant[1:-1])
    else: 
        V=int(variant[:-1])
    p=variant[-1]
    odd=0
    if p=='l':
        odd=1
    
    res=(V*2)+odd
    return s+str(res)

for line in file: 
    #-1064h;-917l;1880l; => 2
    # or
    #-1064h;-917l;1880l; 24h;1l => 2
    if line[0]=='#': 
        # print(line,end='')
        continue
    line=line.strip().rstrip().split(" ")
    if int(line[-1])<1: 
        continue
    for variant in line[0].split(';')[:-1]:
        print (f(variant)+';', end='')
    print ()
    if len(line)==4: # not paired
        for variant in line[1].split(';'):
            print (f(variant)+';', end='')
        print ()

    

    
    