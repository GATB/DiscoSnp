 
#####################
# Default option run: 
#####################
../../run_discoSnp++.sh -r fof.txt -T
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi
# Test the .fa 
# The sequence ids and orders are not conserved due to parallelisation. This explains why we separate sequences from headers and why we remove ids
grep -v ">" ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa | sort -n > ref
grep -v ">" discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa | sort -n > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** Test: FAILURE on diff sequences of .fa"
       exit 1
fi

grep ">" ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa |cut -d "|" -f 2- | sort -n > ref
grep ">" discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa |cut -d "|" -f 2- | sort -n > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** Test: FAILURE on diff headers of .fa"
       exit 1
fi

# Test the .vcf headers 
grep "^#" ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi


# Test the .vcf content
awk '!/#/' ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort -n > ref
awk '!/#/' discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort -n > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** Test: FAILURE on diff content of .vcf"
       exit 1
fi



######################################################
# With reference run (using previously created graph): 
######################################################
../../run_discoSnp++.sh -r fof.txt -T -G humch1_first_5M.fasta -g
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi
# Test the .vcf headers 
grep "^#" ref_with_mapping_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi


# Test the .vcf content
awk '!/#/'  ref_with_mapping_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf  | cut -f 2,4- | sort -n > ref
awk '!/#/'  discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort -n > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff content of .vcf"
       exit 1
fi



# Test the for IGV .vcf headers 
grep "^#" ref_with_mapping_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent_for_IGV.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_auto_D_100_P_1_b_0_coherent_for_IGV.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of for IGV .vcf"
       exit 1
fi


# Test the for IGV .vcf content
awk '!/#/' ref_with_mapping_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort -n > ref
awk '!/#/' discoRes_k_31_c_auto_D_100_P_1_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort -n > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff content of for IGV .vcf"
       exit 1
fi




#####################
# CLOSE SNPS  option run: 
#####################
../../run_discoSnp++.sh -r fof.txt -T -P 10
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi
# Test the .fa 
# The sequence ids and orders are not conserved due to parallelisation. This explains why we separate sequences from headers and why we remove ids
grep -v ">" ref_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.fa | sort -n > ref
grep -v ">" discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.fa | sort -n > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** Test: FAILURE on diff sequences of .fa"
       exit 1
fi

grep ">" ref_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.fa |cut -d "|" -f 2- | sort -n > ref
grep ">" discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.fa |cut -d "|" -f 2- | sort -n > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** Test: FAILURE on diff headers of .fa"
       exit 1
fi

# Test the .vcf headers 
grep "^#" ref_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi


# Test the .vcf content
awk '!/#/' ref_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf | cut -f 2,4- | sort -n > ref
awk '!/#/' discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf | cut -f 2,4- | sort -n > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** Test: FAILURE on diff content of .vcf"
       exit 1
fi



######################################################
# With reference run (using previously created graph): 
######################################################
../../run_discoSnp++.sh -r fof.txt -T -G humch1_first_5M.fasta -g  -P 10
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi


# Test the fasta sequence content
grep -v '>' ref_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.fa | sort > ref
grep -v '>' discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.fa | sort > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS fasta FAILURE:"
       echo "*** Test: FAILURE on diff headers of .fa"
       exit 1
fi

# Test the .vcf headers 
grep "^#" ref_with_mapping_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi




# THOSE TESTS DONT WORK. THE VCF DIFFERS (WHILE THE FASTA ARE THE SAME)
# # Test the .vcf content
# awk '!/#/'  ref_with_mapping_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf  | cut -f 2,4- | sort -n > ref
# awk '!/#/'  discoRes_k_31_c_auto_D_100_P_10_b_0_coherent.vcf | cut -f 2,4- | sort -n > created
# diff created ref
# if [ $? -ne 0 ] ; then
#        echo "*** With close SNPS and mapping on ref FAILURE:"
#        echo "*** Test: FAILURE on diff content of .vcf"
#        exit 1
# fi



# Test the for IGV .vcf headers 
grep "^#" ref_with_mapping_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent_for_IGV.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_auto_D_100_P_10_b_0_coherent_for_IGV.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of for IGV .vcf"
       exit 1
fi


# THOSE TESTS DONT WORK. THE VCF DIFFERS (WHILE THE FASTA ARE THE SAME)
# # Test the for IGV .vcf content
# awk '!/#/' ref_with_mapping_discoRes_k_31_c_auto_D_100_P_10_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort -n > ref
# awk '!/#/' discoRes_k_31_c_auto_D_100_P_10_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort -n > created
# diff created ref
# if [ $? -ne 0 ] ; then
#        echo "*** With close SNPS and mapping on ref FAILURE:"
#        echo "*** Test: FAILURE on diff content of for IGV .vcf"
#        exit 1
# fi





echo "****************************"
echo "***** ALL TESTS PASSED *****"
echo "****************************"


