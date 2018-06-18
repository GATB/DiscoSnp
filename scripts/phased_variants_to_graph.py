# cf https://www.evernote.com/l/ARWYfQTvvoVFL4sOHvjcJ-AUksvCNp1zz-8
        
import sys

def printid(anid):
    if anid[0]=='-':    print(anid[1:],end='')
    else:               print(anid,end='')

def formatid(anid):
    if anid[0]=='-':    return anid[1:]
    else:               return anid


node2nodes={}
node2paired_nodes={}
nodes={}

def id1id2 (id1,id2,coverage,pairend):
    if pairend: struct=node2paired_nodes
    else:       struct=node2nodes
    if id2<id1:
        tmp=id1
        id1=id2
        id2=tmp
    if id1 not in struct:
        struct[id1]={}
    if id2 not in struct[id1]:
        struct[id1][id2]=0
    struct[id1][id2]+=coverage
    if id1 not in nodes: nodes[id1]=0
    if id2 not in nodes: nodes[id2]=0

def get_phasing_edges(file):
    
    # Print phasing edges
    for line in file: #-1064h;-917l;1880l; => 2
        if line[0]=='#':
            continue
        line=line.strip().rstrip().split(' ')
        pairend=False
        if len(line)==3: 
            coverage=int(line[2])
        else:
            coverage=int(line[3])
            pairend=True
            
        if coverage < 1: continue
        
        ids=line[0].split(';')[:-1]
        for i in range(len(ids)-1):
            id1=formatid(ids[i])
            id2=formatid(ids[i+1])
            id1id2 (id1,id2,coverage,False)
            
        if pairend:
            id1=line[0].split(';')[-2]
            id2=line[1].split(';')[0]
            id1id2 (id1,id2,coverage,True)
            ids=line[1].split(';')[:-1]
            for i in range(len(ids)-1):
                id1=formatid(ids[i])
                id2=formatid(ids[i+1])
                id1id2 (id1,id2,coverage,False)
            
            
            
        
def print_phasing_edges():
    for id1 in node2nodes:
        for id2 in node2nodes[id1]:
            print(str(id1)+"\t"+str(id2)+"\t"+str(node2nodes[id1][id2])+"\th")

def print_pairing_edges():
    for id1 in node2paired_nodes:
        for id2 in node2paired_nodes[id1]:
            print(str(id1)+"\t"+str(id2)+"\t"+str(node2paired_nodes[id1][id2])+"\tp")
            
def print_allele_edges():
    for id1 in nodes:
        if id1[-1]=='l': continue
        int_id=id1[:-1]
        id2=int_id+'l'
        if id2 in nodes: 
            print(id1+"\t"+id2+"\t0\ta")
        
if len(sys.argv)>1:
    file = open(sys.argv[1])
else:
    file=sys.stdin

get_phasing_edges(file)    
print("source\ttarget\tcoverage\ttype")
print_phasing_edges()
print_pairing_edges()
print_allele_edges()
            
#
# # Print phasing edges
# for line in file: #-1064h;-917l;1880l; => 2
#     if line[0]=='#':
#         continue
#     line=line.strip().rstrip().split(' ')
#     ids=line[0].split(';')[:-1]
#     coverage=line[2]
#     for i in range(len(ids)-1):
#         printid(ids[i])
#         print("\t",end='')
#         printid(ids[i+1])
#         print("\t",coverage, "\t p")
#
        
        
    
    
    
