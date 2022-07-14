#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import getopt

# Calcule un recall locus-based : vcf disco filtré avec 1 SNP par cluster
# recall-locus = # locus avec un VP / # locus reference
# locus-dup = # locus de la ref qui sont présent au moins 2 fois dans le vcf disco
# compare SNP du génome simulé (argv[1] --> vcf de référence)
# avec ceux d'un vcf DiscoSnp (arg[2])
# 1 prediction  = chr_pos_allele_pair, VP si pos et allele pair identique

def index_reference(ref_vcf):

    ##construit un index avec le vcf de ref, sour la forme :  {"chr3R" : {"12957641" : ["A_T","A_G"]}}  + un index SNP_to_locus {"chr1_12957641" : "loci56681"} + un dico des locus {"loci56681" : 0}
    ## gere les sites multi-alleliques
    ## alleles sortés
    
    #chr3R    12957641    1    T    A    .    loci56681;128    .    ./.
    
    SNP_index = {}
    SNP_to_locus = {}  # for each distinct chr_pos, stores the locus, key is "chr_pos"
    all_loci = {}
    filin = open(ref_vcf, 'r')
    nb_total_loci = 0
   
    while True:
        line = filin.readline()
        if not line: break
        if line.startswith("#"): continue

        chr = line.split("\t")[0]
        pos = line.split("\t")[1]
        ref = line.split("\t")[3].upper()
        alt = line.split("\t")[4].upper()
        
        locus = line.split("\t")[6].split(";")[0]

        if ref<alt:
            bases = '_'.join([ref,alt])
        else:
            bases = '_'.join([alt,ref])
      
        ##initialise index si pas fait
        if chr not in SNP_index:
            SNP_index[chr] = {}
        if pos not in SNP_index[chr]:
            SNP_index[chr][pos] = []
            SNP_to_locus["{0}_{1}".format(chr, pos)] = locus
        if bases not in SNP_index[chr][pos]:
            SNP_index[chr][pos].append(bases)

        if locus not in all_loci:
            all_loci[locus] = 0
            nb_total_loci += 1

    filin.close()
    #print(SNP_index) 
    return SNP_index, SNP_to_locus, all_loci, nb_total_loci

def extract_multiple_pos(chr,pos,XA_field):
    dict_pos = {}
    #Ty=SNP;Rk=0.025231;UL=1;UR=1;CL=.;CR=.;Genome=.;Sd=1;XA=chr2L_+5176594,chr2L_+5177242,chr2L_+5177566

    ## ajouter position primaire
    dict_pos[chr]=[]
    dict_pos[chr].append(pos)
    for i in XA_field.split(","):
        pos = i.split("_")[-1]
        chr = i.split("_")[0]
        if chr not in dict_pos: dict_pos[chr] = []
        dict_pos[chr].append(pos)

    return dict_pos



def comp_disco_vcf(vcf_disco, index, SNP_to_locus, all_loci):

    cpt=0
    nb_total_prediction = 0 # pour vérifier : nb_TP_prediction+nb_FP=nb_total_prediction
    nb_locus_duplicated = 0 # compte le nb de vrais loci qui apparaissent + d'une fois
    nb_locus_recall = 0 # compte le nb de vrais locis trouvés (au moins 1 fois)
    
    filin = open(vcf_disco, 'r')
    while True:
        #random_genome	49694	9997	G	T	.	MULTIPLE	Ty=SNP;Rk=1;UL=83;UR=4;CL=83;CR=4;Genome=G;Sd=1	GT:DP:PL:AD:HQ
        line = filin.readline()
        if not line: break
        if line.startswith("#"): continue
 
        # WARNING : enleve ce test pour la compatibilité avec Stacks, dans nos tests on ne demande pas d'INDEL a disco, outes les lignes sont des ty=SNP.
        #thistype = line.split("\t")[7].split(";")[0].split("=")[1].strip()
        #print(thistype)
        #if thistype != "SNP": continue

        nb_total_prediction += 1
        
        chr = line.split("\t")[0]
        pos = line.split("\t")[1]
        id = line.split("\t")[2]
        ref = line.split("\t")[3].upper()
        alt = line.split("\t")[4].upper()

        if ref<alt:
            bases = '_'.join([ref,alt])
        else:
            bases = '_'.join([alt,ref])
    
        #print("chr"+str(chr)+" pos "+str(pos)+" alleles "+ref+"-"+alt+" geno"+str(pred_geno))
        
#        if line.startswith("SNP"): ##si non mappé
#            nb_FP += 1
#            filout.write("{0}\t{1}\n".format(id, "FP"))
#            continue

        #### cas des alignements multiples ####
        splitXA = line.split("\t")[7].split("XA=")
        if line.split("\t")[6].strip() == "MULTIPLE" and len(splitXA)>1  :

            SNP_trouve = 0
            XA_field = splitXA[1].split(";")[0]
            dict_multiple_pos = extract_multiple_pos(chr,pos,XA_field)
            for chr in dict_multiple_pos:
                if SNP_trouve:
                    break
                if chr not in index:
                    continue

                for pos_mult in dict_multiple_pos[chr]:
                    if SNP_trouve:
                        break # sort de la boucle sur les positions multiples
                    if pos_mult in index[chr]:
                        if len(index[chr][pos_mult])>1 or bases == index[chr][pos_mult][0]: # si multi-allelique compte comme TP sans checker les bases
                            SNP_trouve = 1
                            snp_found="{0}_{1}".format(chr, pos_mult)
                            locus_found = SNP_to_locus[snp_found]
                            if all_loci[locus_found] == 0:
                                nb_locus_recall +=1
                                all_loci[locus_found] = 1
                            elif all_loci[locus_found] == 1:
                                nb_locus_duplicated += 1
                                all_loci[locus_found] = 2

            continue ##on passe à la ligne suivante

        ### cas des alignements uniques ou multiples sans champ XA ###

        if chr in index:
            #print("SNP sur chr indexé :" + chr + "_" + pos)
            if not pos in index[chr]: 
                continue

            if len(index[chr][pos])>1 or bases == index[chr][pos][0]: # si multi-allelique compte comme TP sans checker les bases
                snp_found="{0}_{1}".format(chr, pos)
                locus_found = SNP_to_locus[snp_found]
                if all_loci[locus_found] == 0:
                    nb_locus_recall +=1
                    all_loci[locus_found] = 1
                elif all_loci[locus_found] == 1:
                    nb_locus_duplicated += 1
                    all_loci[locus_found] = 2

    filin.close()
    return nb_locus_recall, nb_locus_duplicated

def main(ref,disco):
    index, SNP_to_locus, all_loci, nb_total_loci = index_reference(ref)
    #print("indexing done\n")
    nb_locus_recall, nb_locus_duplicated = comp_disco_vcf(disco, index, SNP_to_locus, all_loci)

    locus_recall = 100*float(nb_locus_recall)/float(nb_total_loci)
    duplicated_ratio1 = 100*float(nb_locus_duplicated)/float(nb_total_loci)
    duplicated_ratio2 = 100*float(nb_locus_duplicated)/float(nb_locus_recall)
    
    print("#############################################")

    print("#############################################")
    print("Locus Recall on position-allele only. (All multiple predictions at the same position are counted)")
    print(f"Total nb of loci (true vcf) : {nb_total_loci}")
    print(f"Nb loci with a TP SNP : {nb_locus_recall}")
    print(f"Nb duplicated loci (with 2 or more TP SNPs) : {nb_locus_duplicated}")
    print(f"Locus recall : {locus_recall:.2f}")
    print(f"Locus duplication ratio (over all loci) : {duplicated_ratio1:.2f}")
    print(f"Locus duplication ratio (over found loci) : {duplicated_ratio2:.2f}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])



