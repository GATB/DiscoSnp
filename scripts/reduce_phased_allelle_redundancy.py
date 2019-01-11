'''
Merge equal facts: 
#comments blabla toto
-1064h;-917l;1880l; => 2
-1064h;-917l;1880l; => 3
...

becomes
#comments blabla toto
-1064h;-917l;1880l; => 5
'''

import sys

if len(sys.argv)>1:
    file = open(sys.argv[1])
else:
    file=sys.stdin

previous=file.readline().strip()
previous_value=""
for line in file:
    if line.split("=")[0].strip() != previous:
        if previous[0]=='#':
            print (previous) #print first comment only
        else:
            print (previous+"=> "+str(previous_value))
        previous=line.split("=")[0].strip()
        previous_value=int(line.split(">")[-1])
    else:
        previous_value+=int(line.split(">")[-1])
print (previous+" => "+str(previous_value))

    

    
    