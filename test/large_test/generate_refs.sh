 
#####################
# Default option run: 
#####################
../../run_discoSnp++.sh -r fof.txt -T -P 1
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi
mv discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa 
mv discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf


######################################################
# With reference run (using previously created graph): 
######################################################
../../run_discoSnp++.sh -r fof.txt -T -G humch1_first_5M.fasta -g discoRes_k_31_c_3.h5 -P 1
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi
mv discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf
mv discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf

#####################
# CLOSE SNPS  option run: 
#####################
../../run_discoSnp++.sh -r fof.txt -T -P 10
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi
mv discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa
mv discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf

######################################################
# With reference run (using previously created graph): 
######################################################
../../run_discoSnp++.sh -r fof.txt -T -G humch1_first_5M.fasta -g discoRes_k_31_c_3.h5 -P 10
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi
mv discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf 
mv discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf

# create archive
zip -r  data_test_disco.zip  humch1_00096_reads.fasta.gz humch1_00100_reads.fasta.gz humch1_first_5M.fasta ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf 
#(update
scp data_test_disco.zip ppeterlo@scm.gforge.inria.fr:/home/groups/gatb-discosnp/htdocs)

