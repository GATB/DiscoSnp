#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 

import sys
import getopt

# 
# compare SNP du génome simulé (argv[1] --> vcf de référence)
# avec ceux d'un vcf DiscoSnp (arg[2])
# 1 prediction  = chr_pos_allele_pair
# output : argv[2].eval : tabulated file with two columns : id_snp TP/FP
#

def index_reference(ref_vcf):

    ##construit un index avec le vcf de ref, sour la forme :  {"chr1" : {"54345677" : ["A_T","A_G"]}}
    ## gere les sites multi-alleliques
    ## alleles sortés
    
    SNP_index = {}
    info_index = {}  # for each distinct chr_pos, stores the locus and the position in locus : usefull to compute recall as a function of locus position, key is "chr_pos"
    filin = open(ref_vcf, 'r')
    nb_total_SNP_pos = 0  # 1 SNP (vcf line) :  1 position
    nb_total_SNP_distinct_pos = 0  # count distinct positions : multi-allelic SNPs counted only once
    
    while True:
        line = filin.readline()
        if not line: break
        if line.startswith("#"): continue

        chr = line.split("\t")[0]
        pos = line.split("\t")[1]
        ref = line.split("\t")[3].upper()
        alt = line.split("\t")[4].upper()
        
        info = line.split("\t")[6]

        if ref<alt:
            bases = '_'.join([ref,alt])
        else:
            bases = '_'.join([alt,ref])

        nb_total_SNP_pos += 1
        
        ##initialise index si pas fait
        if chr not in SNP_index:
            SNP_index[chr] = {}
        if pos not in SNP_index[chr]:
            SNP_index[chr][pos] = []
            nb_total_SNP_distinct_pos += 1
            info_index["{0}_{1}".format(chr, pos)] = info
        if bases not in SNP_index[chr][pos]:
            SNP_index[chr][pos].append(bases)

    filin.close()
    #print(SNP_index) 
    return SNP_index, nb_total_SNP_pos, nb_total_SNP_distinct_pos, info_index

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



def comp_disco_vcf(vcf_disco, index):
    cpt=0
    dict_TP = {} #créé une liste pour éviter la redondance (recall)
    nb_TP_recall = 0  # nb TP pour le calcul du recall = chaque SNP-geno de la référence est compté qu'une seule fois
    nb_TP_precision = 0  # nb TP pour le calcul de la precision = chaque prédiction SNP-geno (!= ./.) est comptée qu'une seule fois (en TP ou en FP)
    nb_FP = 0  # nb FP pour le calcul de la précision
    nb_total_prediction = 0 # pour vérifier : nb_TP_prediction+nb_FP=nb_total_prediction
    
    filin = open(vcf_disco, 'r')
    output_name = vcf_disco+".eval"
    filout = open(output_name, 'w')
    
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
        
        if line.startswith("SNP"): ##si non mappé
            nb_FP += 1
            filout.write("{0}\t{1}\n".format(id, "FP"))
            continue

        #### cas des alignements multiples ####
        splitXA = line.split("\t")[7].split("XA=")
        if line.split("\t")[6].strip() == "MULTIPLE" and len(splitXA)>1 :

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
                        if len(index[chr][pos_mult])>1: # si multi-allelique compte comme TP sans checker les bases
                            SNP_trouve = 1
                            nb_TP_precision += 1
                            filout.write("{0}\t{1}\n".format(id, "TP"))
                            snp_found="{0}_{1}".format(chr, pos_mult)
                            if snp_found not in dict_TP:
                                dict_TP[snp_found]=1
                                nb_TP_recall += 1
                        elif bases == index[chr][pos_mult][0]:
                            SNP_trouve = 1
                            nb_TP_precision += 1
                            filout.write("{0}\t{1}\n".format(id, "TP"))
                            snp_found="{0}_{1}".format(chr, pos_mult)
                            if snp_found not in dict_TP:
                                dict_TP[snp_found]=1
                                nb_TP_recall += 1

            if SNP_trouve == 0 : ##si aucun SNP n'a été associé aux différentes positions d'alignement, on compte en FP
                nb_FP += 1
                filout.write("{0}\t{1}\n".format(id, "FP"))


            continue ##on passe à la ligne suivante

        ### cas des alignements uniques ou multiples sans champ XA ###

        if chr in index:
            #print("SNP sur chr indexé :" + chr + "_" + pos)
            if not pos in index[chr]: 
                nb_FP += 1
                filout.write("{0}\t{1}\n".format(id, "FP"))
                continue

            if len(index[chr][pos])>1: # si multi-allelique compte comme TP sans checker les bases
                nb_TP_precision += 1
                filout.write("{0}\t{1}\n".format(id, "TP"))
                snp_found="{0}_{1}".format(chr, pos)
                if snp_found not in dict_TP:
                    dict_TP[snp_found]=1
                    nb_TP_recall += 1
            else:
                if bases == index[chr][pos][0]:
                    nb_TP_precision += 1
                    filout.write("{0}\t{1}\n".format(id, "TP"))
                    snp_found="{0}_{1}".format(chr, pos)
                    if snp_found not in dict_TP:
                        dict_TP[snp_found]=1
                        nb_TP_recall += 1
                else:
                    nb_FP += 1
                    filout.write("{0}\t{1}\n".format(id, "FP"))

        else:
            nb_FP += 1
            filout.write("{0}\t{1}\n".format(id, "FP"))
            #print("SNP sur chr non indexé :" + chr + "_" + pos)
       	    

    filin.close()
    filout.close()
    return nb_TP_recall, nb_TP_precision, nb_FP, nb_total_prediction, dict_TP

def main(ref,disco):
    index, nb_total_SNP_pos, nb_total_SNP_distinct_pos, info_index = index_reference(ref)
    nb_TP_recall, nb_TP_precision, nb_FP, nb_total_prediction, dict_TP = comp_disco_vcf(disco, index)


    ## Recall info :
    recall_file_name = disco+".recall"
    # WARNING : assuming loci of 150 nt !!!!
    nb_true_per_pos = [ 0 for i in range(150)]
    nb_found_per_pos = [ 0 for i in range(150)]
    filout = open(recall_file_name, 'w')
    for true_pos in info_index:
        pos_in_locus = int(info_index[true_pos].split(";")[1])
        nb_true_per_pos[pos_in_locus-1] += 1
        found = 0
        if true_pos in dict_TP:
            found = 1
            nb_found_per_pos[pos_in_locus-1] += 1
        filout.write("{0}\t{1}\t{2}\n".format(true_pos,info_index[true_pos], found))
    filout.close()


    stats_file_name = disco+".recall.stats"
    filout2 = open(stats_file_name, 'w')
    for i in range(len(nb_true_per_pos)):
        filout2.write("{0}\t{1}\t{2}\n".format(i+1,nb_true_per_pos[i],nb_found_per_pos[i]))
    filout2.close()

    #print(list_FP)
    #print(index)
    #print(nb_SNP_ref)
    #print(nb_TP)
    print("#############################################")

    print("#############################################")
    print("Precision/recall on position-allele only. All multiple predictions at the same position are counted \nPrecision pos_allele only: {0} %\nRecall pos_allele only: {1} %\nNb SNP_pos total : {2} \nNb SNP_distinct_pos total : {3} \nNb TP_rec : {4}\nNb TP_prec : {5}\nNb FP : {6}\nNb total predictions (vcf lines) : {7}".format((100*float(nb_TP_precision)/float((nb_TP_precision + nb_FP))), (100*float(nb_TP_recall)/float((nb_total_SNP_distinct_pos))), nb_total_SNP_pos, nb_total_SNP_distinct_pos, nb_TP_recall, nb_TP_precision, nb_FP, nb_total_prediction))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])



