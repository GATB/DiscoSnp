import sys


def dostuff(fileVCF,fileFP,variant_type):
    print "index FP"
    # store the FP
    FP_IDs=[]#
    for line in fileFP:
        #SNP id 40
        #INDEL id 262195
        if not line.startswith(variant_type): continue
        id=int(line.split(" ")[2].strip())
        if id not in FP_IDs:
            FP_IDs.append(id)
        
    print "detects FP/VP VS MULTIPLE / SINGLE"
    
    if variant_type=="INDEL": variant_type="INS"
    # for each prediction: check if its a MULTIPLE OR NOT and IF ITS A FP or NOT
    SINGLE_TP=0
    SINGLE_FP=0
    UNMAPPED_TP=0
    UNMAPPED_FP=0
    MULTIPLE_TP=0
    MULTIPLE_FP=0
    i=0
    for line in fileVCF:
        i+=1
        if i%5000==0:
            print "For ",variant_type
            print "single = ",SINGLE_FP+SINGLE_TP," TP: ",SINGLE_TP, " FP : ", SINGLE_FP 
            print "multiple = ",MULTIPLE_FP+MULTIPLE_TP," TP: ",MULTIPLE_TP, " FP : ", MULTIPLE_FP 
            print "unmapped = ",UNMAPPED_TP+UNMAPPED_FP," TP: ",UNMAPPED_TP, " FP : ", UNMAPPED_FP,"\n"
        
        #gi|224384768|gb|CM000663.1|     105870192       200320  CC      C       .       PASS    Ty=INS;Rk=0.68919;MULTI=.;DT=0;UL=.;UR=.;CL=.;CR=.;C1=9,24;C2=33,2;Genome=.;Sd=-1       GT:DP:PL        0/1:33:0        0/0:35:0
        #gi|224384768|gb|CM000663.1|     145119644       54537   A       G       .       MULTIPLE        Ty=SNP;Rk=0.68912;MULTI=multi;DT=-1;UL=.;UR=.;CL=.;CR=.;C1=52,0;C2=21,42;Genome=A;Sd=-1 GT:DP:PL        0/0:52:0        0/1:63:0
        if line.startswith("#"): continue
        if line.startswith("0"): continue #TODO: BUG a indiquer a CHLOE
        if not line.split("\t")[7].split(";")[0].split("=")[1]==variant_type: continue
        TP=True
        Mapped=""
        id=int(line.split("\t")[2].split("_")[0].strip())
        if id in FP_IDs:
            TP=False
    
        if line.split("\t")[6].strip() == "MULTIPLE":
            Mapped="M" # Multiple
    
        if line.split("\t")[6].strip() == "PASS":
            Mapped="S" # Single
    
    
        if line.split("\t")[6].strip() == ".":
            Mapped="U" # Unmapped
    
        if Mapped=="S":
            if TP: SINGLE_TP+=1
            else: SINGLE_FP+=1
    
        if Mapped=="M":
            if TP: MULTIPLE_TP+=1
            else: MULTIPLE_FP+=1
        
        if Mapped=="U":
            if TP: UNMAPPED_TP+=1
            else: UNMAPPED_FP+=1
    

    print "Finally for ",variant_type
    print "single = ",SINGLE_FP+SINGLE_TP," TP: ",SINGLE_TP, " FP : ", SINGLE_FP 
    print "multiple = ",MULTIPLE_FP+MULTIPLE_TP," TP: ",MULTIPLE_TP, " FP : ", MULTIPLE_FP 
    print "unmapped = ",UNMAPPED_TP+UNMAPPED_FP," TP: ",UNMAPPED_TP, " FP : ", UNMAPPED_FP,"\n"



fileVCF = open(sys.argv[1], 'r') # open the VCF FILE
fileFP = open(sys.argv[2], 'r') # open the false positive file
dostuff(fileVCF,fileFP,"SNP")
fileVCF.close()
fileFP.close()
fileVCF = open(sys.argv[1], 'r') # open the VCF FILE
fileFP = open(sys.argv[2], 'r') # open the false positive file
dostuff(fileVCF,fileFP,"INDEL")
fileVCF.close()
fileFP.close()