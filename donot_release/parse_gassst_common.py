    
def findNbMachedPolymorphism(list_reference_polymorphism, start, stop, predicted_positions):
    nbMapped=0 # True positive prediction 
    rev_repredicted_positions=[]
    for pos in predicted_positions:
        rev_repredicted_positions.append(stop-start-pos-1) # stores also the reverse positions (in case the mapping ins in the reverse order)
    for pos in range(start, stop+1):
        if pos in list_reference_polymorphism and list_reference_polymorphism[pos]==False:
            if pos-start-1 in predicted_positions or pos-start-1 in rev_repredicted_positions or pos-start in predicted_positions or pos-start in rev_repredicted_positions :
            #if pos-start in predicted_positions or pos-start in rev_repredicted_positions :
                list_reference_polymorphism[pos]=True
                nbMapped+=1

    return nbMapped
    

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
        line=filin.readline() #>SNP_higher_path_95|high|nb_pol_1|C1_0|C2_40|rank_1.00000 OR >INDEL_9_higher_path_891|high|nb_pol_1|C1_0|C2_30|rank_1.00000
        if not line: break
        if line.startswith(">"+typePolymorphism):
            pol_id=int(line.split("|")[0].split("_")[-1])
            rank=float(line.split("|")[-1].split("_")[-1])                        
            if rank>=threshold:
                nb_pol+=int(line.split("|")[2].split("_")[-1])
                list_predicted_polymorphism[pol_id] = False # This polymorphism has not been confirmed yet
            upper_path=filin.readline() # read sequence upper, we don't care
            filin.readline() # read comment lower, we don't care
            lower_path=filin.readline() # read sequence lower, we don't care
            
            # We store the SNP positions for this bubble
            if rank>=threshold and typePolymorphism!="INDEL":
                list_predicted_polymorphism_positions[pol_id] = []
                for i in range(len(upper_path)):
                    if upper_path[i] != lower_path[i]:
                        list_predicted_polymorphism_positions[pol_id].append(i)
            
    return list_predicted_polymorphism,nb_pol,list_predicted_polymorphism_positions
    
####################################################################
#                         IOs
####################################################################
def print_results(nb_predicted, list_reference, polymorphism, threshold):
   nbTP=0
   for chro in list_reference:
       for i in list_reference[chro]:
           if list_reference[chro][i]:
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