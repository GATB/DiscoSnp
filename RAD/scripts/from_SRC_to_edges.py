import sys
if len(sys.argv)>1:
    file = open(sys.argv[1])
else:
    file=sys.stdin

#IN= "9976:4242 9976 9604 "
#OUT= 
# "9976 4242"
# "9976 9604"
def printLine(line):
    if line[0]=="#":
        print (line)
        return
    id_source=line.split(":")[0]
    ids_children=line.split(":")[1].split(" ")
    for id_child in ids_children:
        if id_child!=id_source:
            print (id_source,id_child)

for line in file:
    printLine(line.strip().rstrip())
        