import sys
import gzip

if len(sys.argv)<4:
    print ("This tool filters out discoSnp prediction whose number of read sets covering it is lower than a user defined threshold. A set covers a prediction if its coverage in at least one of the two alleles is higher than a user defined threshold")
    print ("python filter_out_using_ratio_of_covered_files.py \".fa from discoSnp\" \"number of sets threshold\" \"minimal coverage\"")
    sys.exit()


if "gz" in sys.argv[1]:
    coherent_file=gzip.open(sys.argv[1],"r")
else: 
    coherent_file=open(sys.argv[1],"r")
ratio_threshold = float(sys.argv[2])
minimal_coverage = int(sys.argv[3])

first_coverage_field_id=0
number_of_read_sets=0
while True:
    #>SNP_higher_path_1384|P_1:30_A/G|high|nb_pol_1|left_unitig_length_108|right_unitig_length_1156|C1_0|C2_0|C3_0|C4_0|C5_183|Q1_0|Q2_0|Q3_0|Q4_0|Q5_70|G1_1/1:19724,2972,47|G2_1/1:21024,3168,50|G3_1/1:17124,2581,42|G4_1/1:19564,2948,47|G5_0/1:16063,1163,1575|rank_0.36839
    
    comment1 = coherent_file.readline()
    if not comment1: break
    splitted_comment1 = comment1.split("|")
    path1 = coherent_file.readline()
    comment2 = coherent_file.readline()
    splitted_comment2 = comment2.split("|")
    path2 = coherent_file.readline()

    coverage_high=[]
    coverage_low= []    
    if first_coverage_field_id==0:
        while True:
            if splitted_comment1[first_coverage_field_id][0]=="C": break
            first_coverage_field_id+=1
            
    
    i=first_coverage_field_id
    while True:
        if splitted_comment1[i][0]!="C": break # no more a coverage
        coverage_high.append(int(splitted_comment1[i].split("_")[1]))
        coverage_low.append( int(splitted_comment2[i].split("_")[1]))
        i+=1
    
    if number_of_read_sets == 0:
        number_of_read_sets = i-first_coverage_field_id
    
    # print comment1,
    #print coverage_high
    #print coverage_low
    number_of_covered_sets=0
    for i in range(len(coverage_high)):
        if coverage_high[i]>=minimal_coverage or coverage_low[i]>=minimal_coverage: number_of_covered_sets+=1
    
    if 100*number_of_covered_sets/float(number_of_read_sets)>=ratio_threshold:
        print (comment1,path1,comment2,path2,)
    
      
            
        
   

