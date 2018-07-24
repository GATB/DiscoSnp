import sys
# coding: utf-8

# In[ ]:




# In[29]:

def store_clusters(cluster_file):
    clusters=open(cluster_file,"r")
    read_id_to_cluster_id={}
    cluster_id_to_cluster_size={}
    cluster_id=-1
    for cluster in clusters:
        # a line is "70166 70345 70409 70222 70406 70167 70223 69786 70407 69787 70408 70611 70610 70344 "
        cluster_id+=1
        cluster_id_to_cluster_size[cluster_id]=len(cluster.rstrip().split())
        for read_id in cluster.rstrip().split():
            read_id_to_cluster_id[int(read_id.split('-')[0])]=cluster_id # A line can be formated as 70166 70345-info_about_similarity
    clusters.close()
    return read_id_to_cluster_id, cluster_id_to_cluster_size
        
def assign_cluster_id_to_sequence_and_print(fasta_file, read_id_to_cluster_id,cluster_id_to_cluster_size):
    fasta=open(fasta_file,"r")
    read_id=-1
    while True: 
        read_id+=1
        comment=fasta.readline().rstrip()
        #>SNP_higher_path_9996|P_1:30_C/T|high|nb_pol_1|left_unitig_length_8|right_unitig_length_5|C1_55|C2_0|C3_0|C4_0|C5_0|Q1_0|Q2_0|Q3_0|Q4_0|Q5_0|G1_0/0:7,170,1104|G2_1/1:1084,167,7|G3_1/1:1064,164,7|G4_1/1:984,152,6|G5_1/1:1264,194,7|rank_1
        if not comment: break            
        sequence=fasta.readline().rstrip()
        #gcggaatgAATTAGTGGTATGTCAAGAGGGACTGCTATCAACACTTACGTAGTGCACATATTTCTTTGCatcgc
        SNP_id=comment.split("|")[0][1:] #SNP_higher_path_9996
        if read_id not in read_id_to_cluster_id:
            print("Warning, read id "+str(read_id)+" not in clusters",file=sys.stderr)
            print(comment,file=sys.stderr)
            print(sequence,file=sys.stderr)
            SNP_id=">cluster_-1_size_1_"+SNP_id
        else:
            cluster_id=read_id_to_cluster_id[read_id]
            SNP_id=">cluster_"+str(cluster_id)+"_size_"+str(cluster_id_to_cluster_size[cluster_id])+"_"+SNP_id
        comment_new=SNP_id
        for stuff in comment.split("|")[1:]:
            comment_new+="|"+stuff
        print (comment_new)
        print (sequence)
        
        
fasta_file=sys.argv[1]
cluster_file=sys.argv[2]
read_id_to_cluster_id,cluster_id_to_cluster_size = store_clusters(cluster_file)
assign_cluster_id_to_sequence_and_print(fasta_file,read_id_to_cluster_id,cluster_id_to_cluster_size)

