# Get and decompress the data
rm -f data_test_disco.zip
wget http://gatb-discosnp.gforge.inria.fr/data_test_disco.zip 

# create: 
# zip -r  data_test_disco.zip  humch1_00096_reads.fasta.gz humch1_00100_reads.fasta.gz humch1_first_5M.fasta ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf 
#(update = scp data_test_disco.zip ppeterlo@scm.gforge.inria.fr:/home/groups/gatb-discosnp/htdocs)
unzip data_test_disco.zip

# Create the file of file: 
ls humch1_00* > fof.txt
 
sh local_large_test.sh


rm -f created ref discoRes*
rm -f data_test_disco.zip
rm -f humch1_*
rm -f fof.txt
rm -f ref*



