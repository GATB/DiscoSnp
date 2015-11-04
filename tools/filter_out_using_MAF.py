import sys



if len(sys.argv)<3:
    print "This tool filters out discoSnp prediction having a minor allele frequency lower than a provided threshold for ALL datasets."
    print "python filter_out_using_MAF.py \".fa from discoSnp\" \"MAF threshold\""
    sys.exit()

coherent_file=open(sys.argv[1],"r")
maf_threshold = float(sys.argv[2])


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
    i=6
    while True:
        if splitted_comment1[i][0]!="C": break # no more a coverage
        coverage_high.append(int(splitted_comment1[i].split("_")[1]))
        coverage_low.append( int(splitted_comment2[i].split("_")[1]))
        i+=1
    # print comment1,
    # print coverage_high
    
    to_output=False
    for i in range(len(coverage_high)):
        if coverage_high[i]==0 and coverage_low[i]==0: continue
        if (min(coverage_high[i],coverage_low[i]) / float(max(coverage_high[i],coverage_low[i]))) < maf_threshold:
            to_output=True
            break
    
    if to_output:
        print comment1,path1,comment2,path2,
    
      
            
        
   

