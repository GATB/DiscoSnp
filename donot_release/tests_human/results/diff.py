filin = open("FN_indels", 'r')
# filin = open("tmp", 'r')

for line in filin:
    #INDEL 246808596 gi AG A
    print len(line.split(" ")[3])-len(line.split(" ")[4])+1