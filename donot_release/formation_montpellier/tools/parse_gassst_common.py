

def homo_or_hereto(text_genotype):
    # 0/1: returns 0 (hetero) or 0|0
    # 1/1: returns 1 (homo)
    # 0/0: returns 1 (homo)
    sep="|" 
    if "/" in text_genotype: 
        sep="/"
    
    
    if text_genotype.split(sep)[0].strip()==text_genotype.split(sep)[1].strip():
        return 1
    return 0
  
def findNbMachedPolymorphism(list_reference_polymorphism, start, stop, predicted_positions, typePolymorphism, forward):
    nbMapped=0 # True positive prediction 
    if forward:
        rev_repredicted_positions=[]
        for pos in predicted_positions:
            rev_repredicted_positions.append(stop-start-pos-1) # stores also the reverse positions (in case the mapping ins in the reverse order)
        predicted_positions=rev_repredicted_positions
    for pos in range(start, stop+1):
        if pos in list_reference_polymorphism and list_reference_polymorphism[pos]==False:
            # we check that the genomic position matches the polymorphism found.
            # In case of SNP, this is a precise position (checked with 0 and 1 index)
            # In case of INDEL: we don't know preciselly where the breakpoint is. Thus, the simulated indel position matched are systematically considered as TP. 
            if pos-start in predicted_positions or pos-start-1 in predicted_positions or pos-start+1 in predicted_positions:
                list_reference_polymorphism[pos]=True
                nbMapped+=1
                # In case of indels: we stop the matching as soon as one indel was found (else we detect too much TP)
                # We maybe could improve this as we are not sure this indel is the one that was predicted. (It may affect up to 0.3% of the TP)
                if typePolymorphism == "INDEL": 
                    break

    if nbMapped > len(predicted_positions):
        print typePolymorphism,"WARNING: nb TP",nbMapped,"is bigger than the number of predictions",len(predicted_positions) 
    return nbMapped
    

def findNbMachedPolymorphism_with_genotype(list_reference_polymorphism, list_reference_genotypes, start, stop, predicted_positions, predicted_polymorpshimes, typePolymorphism, forward):
    nbMapped=0 # True positive prediction 
    this_nb_good_genotype_predictions=0
    this_nb_bad_genotype_predictions=0
    if forward:
        rev_repredicted_positions=[]
        for pos in predicted_positions:
            rev_repredicted_positions.append(stop-start-pos-1) # stores also the reverse positions (in case the mapping ins in the reverse order)
        predicted_positions=rev_repredicted_positions
    for pos in range(start, stop+1):
        if pos in list_reference_polymorphism and list_reference_polymorphism[pos]==False:
            # we check that the genomic position matches the polymorphism found.
            # In case of SNP, this is a precise position (checked with 0 and 1 index)
            # In case of INDEL: we don't know preciselly where the breakpoint is. Thus, the simulated indel position matched are systematically considered as TP. 
            if pos-start in predicted_positions or pos-start-1 in predicted_positions or pos-start+1 in predicted_positions:
                list_reference_polymorphism[pos]=True
                
                ############### VALID GENOTYPE                #
                # print list_reference_genotypes[pos] #0|0 1|1
                # print predicted_polymorpshimes #[1, 1]
                for i in range(len(predicted_polymorpshimes)):                    #
                    # print i, len(predicted_polymorpshimes),"ref", list_reference_genotypes[pos].split(" ")[i],
                    # print "predicted", predicted_polymorpshimes[i]
                    if homo_or_hereto(list_reference_genotypes[pos].split(" ")[i])==predicted_polymorpshimes[i]:
                        this_nb_good_genotype_predictions+=1
                    else:
                        this_nb_bad_genotype_predictions+=1
                # print nb_good_genotype_predictions,nb_bad_genotype_predictions
                ##############################
                
                
                nbMapped+=1
                # In case of indels: we stop the matching as soon as one indel was found (else we detect too much TP)
                # We maybe could improve this as we are not sure this indel is the one that was predicted. (It may affect up to 0.3% of the TP)
                if typePolymorphism == "INDEL": 
                    break

    if nbMapped > len(predicted_positions):
        print typePolymorphism,"WARNING: nb TP",nbMapped,"is bigger than the number of predictions",len(predicted_positions) 
    return nbMapped, this_nb_good_genotype_predictions, this_nb_bad_genotype_predictions
    

####################################################################
#                         PREDICTIONS
####################################################################
# Compute Number of predicted polymorphism
def getPredictedPolymorphism (discoResultsFile, threshold, typePolymorphism):
    list_predicted_polymorphism={}
    list_predicted_polymorphism_positions={} # for each polymorphism, contains the position of the polymorphism.
    filin = open(discoResultsFile, 'r') 
    nb_pol=0
    while 1: 
        line=filin.readline() #>INDEL_lower_path_99960|P_1:30_3_2|high|nb_pol_1|C1_24|C2_0|G1_1/1|G2_0/0|rank_1.00000 OR >SNP_higher_path_99995|P_1:30_C/T|high|nb_pol_1|C1_0|C2_38|G1_1/1|G2_0/0|rank_1.00000
                              #>SNP_higher_path_9979|P_1:30_A/G|high|nb_pol_1|left_unitig_length_48|right_unitig_length_4|left_contig_length_48|right_contig_length_1781|C1_33|C2_0|G1_0/0:6,104,664|G2_1/1:704,110,6|rank_1.00000
        if not line: break
        if line.startswith(">"+typePolymorphism):
            pol_id=int(line.split("|")[0].split("_")[-1])
            if "rank" in line :
                rank=float(line.split("|")[-1].split(" ")[0].split("_")[-1])                      
            else:
                rank=1
             
            upper_path=filin.readline() # read sequence upper, we don't care
            filin.readline() # read comment lower, we don't care
            lower_path=filin.readline() # read sequence lower, we don't care
            
        
            if rank>=threshold:
                if typePolymorphism=="INDEL":
                    nb_pol+=1
                    
                if typePolymorphism=="SNP":# We store the SNP positions for this bubble
                    list_predicted_polymorphism_positions[pol_id] = []
                    for i in range(len(upper_path)):
                        if upper_path[i] != lower_path[i]:
                            list_predicted_polymorphism_positions[pol_id].append(i)
                            nb_pol+=1
                
                list_predicted_polymorphism[pol_id] = 0 # This polymorphism has not been confirmed yet
            
            
            
    return list_predicted_polymorphism,nb_pol,list_predicted_polymorphism_positions
    
####################################################################
#                         IOs
####################################################################

def printRoc (discoResultsFile, list_reference, list_predicted_polymorphism, threshold, typePolymorphism):
    nb_simulated=0
    for chro in list_reference:
        nb_simulated+=len(list_reference[chro])
        
    nbTP=0
    nb_predicted=0
                
   
    filin = open(discoResultsFile, 'r') 
    
    while 1: 
        line=filin.readline() #>SNP_higher_path_95|high|nb_pol_1|C1_0|C2_40|rank_1.00000 OR >INDEL_9_higher_path_891|high|nb_pol_1|C1_0|C2_30|rank_1.00000
        if not line: break
        if line.startswith(">"+typePolymorphism):
            nb_pol=0
            pol_id=int(line.split("|")[0].split("_")[-1])
            if "rank" in line :
                rank=float(line.split("|")[-1].split("_")[-1])                      
            else:
                rank=1
             
            upper_path=filin.readline() # read sequence upper, we don't care
            filin.readline() # read comment lower, we don't care
            lower_path=filin.readline() # read sequence lower, we don't care
            
             
            if rank>=threshold:
                if typePolymorphism=="INDEL":
                    nb_pol+=1
                    
                if typePolymorphism=="SNP":# We store the SNP positions for this bubble
                    for i in range(len(upper_path)):
                        if upper_path[i] != lower_path[i]:
                            nb_pol+=1
                
                nbTP+=list_predicted_polymorphism[pol_id]
                nb_predicted+=nb_pol
                
                
                
                
                if nb_predicted>0 : print "%.2f "%  float(nbTP/float(nb_simulated)*100)+"%.2f "% float(nbTP/float(nb_predicted)*100) +"%.2f "% rank
            
                # if list_predicted_polymorphism[pol_id] != nb_pol:
               #      print "FP: ",pol_id
               #  else:
               #      print
            
    

def print_results(nb_predicted, list_reference, polymorphism, threshold):
   nbTP=0
   # if polymorphism=="INDEL": print list_reference
   for chro in list_reference:
       for i in list_reference[chro]:
           if list_reference[chro][i]==True:
               nbTP+=1

   nb_simulated=0
   for chro in list_reference:
       nb_simulated+=len(list_reference[chro])

               
        
   print "-------------------------------"
   print "             ",polymorphism
   if threshold>0:
       print " from",polymorphism,"predicted with a rank bigger or equal to",threshold,":"
   print nb_simulated, polymorphism,"in the reference.\t Among them", nbTP, "are correctly predicted"
   print nb_predicted, polymorphism,"were predicted.\t Among them", nbTP, "are correctly mapped"
   if nb_predicted>0: print polymorphism,"precision\t %.2f"% float(nbTP/float(nb_predicted)*100)
   if len(list_reference)>0: print polymorphism,"recall\t %.2f"% float(nbTP/float(nb_simulated)*100)
   print "-------------------------------"

# def print_few_results(nb_predicted, list_reference):
#     nbTP=0
#     # if polymorphism=="INDEL": print list_reference
#     for chro in list_reference:
#         for i in list_reference[chro]:
#             if list_reference[chro][i]==True:
#                 nbTP+=1
#
#     nb_simulated=0
#     for chro in list_reference:
#         nb_simulated+=len(list_reference[chro])
#
#     if nb_predicted>0 and len(list_reference)>0: print "%.2f "%  float(nbTP/float(nb_simulated)*100)+"%.2f "% float(nbTP/float(nb_predicted)*100)


def FN(list_reference, polymorphism):
   print "-----------------------------------------------"
   print polymorphism,"FALSE NEGATIVE GENOME POSITIONS "
   for chro in list_reference:
       for i in list_reference[chro]:
           if not list_reference[chro][i]:
               print "FN",polymorphism, "chr",chro," position", i
   print "-----------------------------------------------"


def FP(list_predicted, polymorphism):
   print "-----------------------------------------------"
   print polymorphism,"FALSE POSITIVES PREDICTION IDs "
   for i in list_predicted:
       if not list_predicted[i]:
           print "FP",polymorphism, "id", i
   print "-----------------------------------------------"