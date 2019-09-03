import sys



"""
Dat file: 
AcorusCalamus_multigraph.dat

data;
param Title := AcorusCalamus;
param start := "3_0_F";
set V :=
3_0_F
1_0_R
0_0_R
0_0_F
1_0_F
1_1_F
3_0_R
1_1_R
;
param w :=
3_0_F	83691
1_0_R	26025
0_0_R	18476
0_0_F	18476
1_0_F	26025
1_1_F	26025
3_0_R	83691
1_1_R	26025
;
set reverse :=
3_0_F	3_0_R
1_0_R	1_0_F
0_0_R	0_0_F
0_0_F	0_0_R
1_0_F	1_0_R
1_1_F	1_1_R
3_0_R	3_0_F
1_1_R	1_1_F
;
set Edges :=
3_0_F	1_1_F	overlaps
3_0_F	1_1_F	links
3_0_F	1_0_F	overlaps
3_0_F	1_0_F	links
1_0_R	3_0_F	overlaps
1_0_R	3_0_F	links
1_0_R	3_0_R	overlaps
1_0_R	3_0_R	links
0_0_R	1_0_R	overlaps
0_0_R	1_0_R	links
0_0_R	1_1_R	overlaps
0_0_R	1_1_R	links
0_0_F	1_0_R	overlaps
0_0_F	1_0_R	links
0_0_F	1_1_R	overlaps
0_0_F	1_1_R	links
1_0_F	0_0_R	overlaps
1_0_F	0_0_R	links
1_0_F	0_0_F	overlaps
1_0_F	0_0_F	links
1_1_F	0_0_R	overlaps
1_1_F	0_0_R	links
1_1_F	0_0_F	overlaps
1_1_F	0_0_F	links
3_0_R	1_1_F	overlaps
3_0_R	1_1_F	links
3_0_R	1_0_F	overlaps
3_0_R	1_0_F	links
1_1_R	3_0_F	overlaps
1_1_R	3_0_F	links
1_1_R	3_0_R	overlaps
1_1_R	3_0_R	links
;
param l :=
3_0_F	1_1_F	overlaps	-99
3_0_F	1_1_F	links	0
3_0_F	1_0_F	overlaps	-99
3_0_F	1_0_F	links	0
1_0_R	3_0_F	overlaps	-99
1_0_R	3_0_F	links	0
1_0_R	3_0_R	overlaps	-99
1_0_R	3_0_R	links	0
0_0_R	1_0_R	overlaps	-99
0_0_R	1_0_R	links	0
0_0_R	1_1_R	overlaps	-99
0_0_R	1_1_R	links	0
0_0_F	1_0_R	overlaps	-99
0_0_F	1_0_R	links	0
0_0_F	1_1_R	overlaps	-99
0_0_F	1_1_R	links	0
1_0_F	0_0_R	overlaps	-99
1_0_F	0_0_R	links	0
1_0_F	0_0_F	overlaps	-99
1_0_F	0_0_F	links	0
1_1_F	0_0_R	overlaps	-99
1_1_F	0_0_R	links	0
1_1_F	0_0_F	overlaps	-99
1_1_F	0_0_F	links	0
3_0_R	1_1_F	overlaps	-99
3_0_R	1_1_F	links	0
3_0_R	1_0_F	overlaps	-99
3_0_R	1_0_F	links	0
1_1_R	3_0_F	overlaps	-99
1_1_R	3_0_F	links	0
1_1_R	3_0_R	overlaps	-99
1_1_R	3_0_R	links	0
;
"""

def print_header():
    print("data;")
    print("param Title := ???;")
    print("param start := \"???\";")


def print_nodes(gfa_file_name):
    print("set V :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if not line: break
        if line[0]=="S":
            "S       0       agaTAATATATGACTAAATGTTAAACTAAAATGAAAAAAAAAACATACATATGTAATGTATTAAGGTTGTAAGGTATAAATGACTGGAATTGCCAAAACTTTCCctgcaaaaAAATCTTGATAAAGTCTGTCTTGTAATGCAAAATGTCCACTTTTGCCCAGTGCGACATACTCAttat     FC:i:83 RC:i:36 SP:0_112;65_179;        BP:0_101;-36_107;       -100020h;212163h;"
            print("p"+line.split()[1])  # 'p' stands for "plus strand"
    print(";")
    gfa_file.close()
    
def print_nodes_weight(gfa_file_name):
    print("param w :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            "S       0       agaTAATATATGACTAAATGTTAAACTAAAATGAAAAAAAAAACATACATATGTAATGTATTAAGGTTGTAAGGTATAAATGACTGGAATTGCCAAAACTTTCCctgcaaaaAAATCTTGATAAAGTCTGTCTTGTAATGCAAAATGTCCACTTTTGCCCAGTGCGACATACTCAttat     FC:i:83 RC:i:36 SP:0_112;65_179;        BP:0_101;-36_107;       -100020h;212163h;"
            print("p"+line.split()[1]+"\t"+line.split()[4].split(":")[-1])
    print(";")
    gfa_file.close()


def print_reverse(gfa_file_name):
    print("set reverse :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="S":
            "S       0       agaTAATATATGACTAAATGTTAAACTAAAATGAAAAAAAAAACATACATATGTAATGTATTAAGGTTGTAAGGTATAAATGACTGGAATTGCCAAAACTTTCCctgcaaaaAAATCTTGATAAAGTCTGTCTTGTAATGCAAAATGTCCACTTTTGCCCAGTGCGACATACTCAttat     FC:i:83 RC:i:36 SP:0_112;65_179;        BP:0_101;-36_107;       -100020h;212163h;"
            print("p"+line.split()[1]+"\t"+"m"+line.split()[1])  # 'p' stands for "plus strand"
    print(";")
    gfa_file.close()
    
    
def print_edges(gfa_file_name):
    print("set Edges :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="L":
            "L      1       -       29384   +       8M"
            "to"
            "m1	p29384	overlaps"
            overlap_len = int(line.split()[5].rstrip("M"))
            if overlap_len>=0:
                sign_source="p"
                if line.split()[2]=='-': sign_source="m"
                sign_target="p"
                if line.split()[4]=='-': sign_target="m"
                type="overlaps"
                if overlap_len==0: type="links"
                print(sign_source+line.split()[1]+"\t"+sign_target+line.split()[3]+"\t"+type)  
            # do not print other edges (unitig linked and snp linked)
    print(";")
    gfa_file.close()
    
    
def print_edges_content(gfa_file_name):
    print("set Edges :=")
    gfa_file = open(gfa_file_name)
    for line in gfa_file.readlines():
        line=line.strip()
        if line[0]=="L":
            "L      1       -       29384   +       8M"
            "to"
            "m1	p29384	overlaps 8"
            
            overlap_len = int(line.split()[5].rstrip("M"))
            if overlap_len>=0:
                sign_source="p"
                if line.split()[2]=='-': sign_source="m"
                sign_target="p"
                if line.split()[4]=='-': sign_target="m"
                type="overlaps"
                if overlap_len==0: type="links"
                print(sign_source+line.split()[1]+"\t"+sign_target+line.split()[3]+"\t"+type+"\t"+str(overlap_len))  
            # do not print other edges (unitig linked and snp linked)
    print(";")
    gfa_file.close()
    
    

def main(gfa_file_name):
    '''
    Creation of a DAT file from a GFA file
    '''
    
    print_header()
    print_nodes(gfa_file_name)
    print_nodes_weight(gfa_file_name)
    print_reverse(gfa_file_name)
    print_edges(gfa_file_name)
    print_edges_content(gfa_file_name)



if __name__ == "__main__":
     main(sys.argv[1])
