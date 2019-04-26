 
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
grep -v ">" ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa | sort > ref
grep -v ">" discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa | sort > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** Test: FAILURE on diff sequences of .fa"
       exit 1
fi

grep ">" ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa |cut -d "|" -f 2- | sort > ref
grep ">" discoRes_k_31_c_3_D_100_P_1_b_0_coherent.fa |cut -d "|" -f 2- | sort > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** Test: FAILURE on diff headers of .fa"
       exit 1
fi

# Test the .vcf headers 
grep "^#" ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi


# Test the .vcf content
awk '!/#/' ref_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort > ref
awk '!/#/' discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort > created
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
grep "^#" ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi


# Test the .vcf content
awk '!/#/'  ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf  | cut -f 2,4- | sort > ref
awk '!/#/'  discoRes_k_31_c_3_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff content of .vcf"
       exit 1
fi



# Test the for IGV .vcf headers 
grep "^#" ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of for IGV .vcf"
       exit 1
fi


# Test the for IGV .vcf content
awk '!/#/' ref_with_mapping_discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort > ref
awk '!/#/' discoRes_k_31_c_3_D_100_P_1_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort > created
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
grep -v ">" ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa | sort > ref
grep -v ">" discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa | sort > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** Test: FAILURE on diff sequences of .fa"
       exit 1
fi

grep ">" ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa |cut -d "|" -f 2- | sort > ref
grep ">" discoRes_k_31_c_3_D_100_P_10_b_0_coherent.fa |cut -d "|" -f 2- | sort > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** Test: FAILURE on diff headers of .fa"
       exit 1
fi

# Test the .vcf headers 
grep "^#" ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi


# Test the .vcf content
awk '!/#/' ref_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf | cut -f 2,4- | sort > ref
awk '!/#/' discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf | cut -f 2,4- | sort > created
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
# Test the .vcf headers 
grep "^#" ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of .vcf"
       exit 1
fi


# Test the .vcf content
awk '!/#/'  ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf  | cut -f 2,4- | sort > ref
awk '!/#/'  discoRes_k_31_c_3_D_100_P_10_b_0_coherent.vcf | cut -f 2,4- | sort > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff content of .vcf"
       exit 1
fi



# Test the for IGV .vcf headers 
grep "^#" ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf | grep -v filedate > ref
grep "^#" discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf | grep -v filedate > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff headers of for IGV .vcf"
       exit 1
fi


# Test the for IGV .vcf content
awk '!/#/' ref_with_mapping_discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort > ref
awk '!/#/' discoRes_k_31_c_3_D_100_P_10_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With close SNPS and mapping on ref FAILURE:"
       echo "*** Test: FAILURE on diff content of for IGV .vcf"
       exit 1
fi

######################################################
#           test the radseq option (-x)              #
######################################################

../../run_discoSnpRad.sh -r fof_rad.txt -k 31 -b 2 -D 0 -P 4 -p radtest 

if [ $? -ne 0 ] ; then
       echo "*** With truncated bubbles FAILURE:"
       echo "*** discoSnp failure"
       exit 1
fi

# Test the .fa
# The sequence ids and orders are not conserved due to parallelisation. This explains why we separate sequences from headers and why we remove ids
grep -v ">" rad_option_test/radref_k_31_c_3_D_0_P_4_b_2_coherent.fa | sort -n > radref
grep -v ">" radtest_k_31_c_3_D_0_P_4_b_2_coherent.fa | sort -n > radtest
diff radref radtest
if [ $? -ne 0 ] ; then
       echo "*** With truncated bubbles FAILURE:"
       echo "*** Test: FAILURE on diff sequences of .fa"
       exit 1
fi

grep ">" rad_option_test/radref_k_31_c_3_D_0_P_4_b_2_coherent.fa |cut -d "|" -f 2- | sort -n > radref
grep ">" radtest_k_31_c_3_D_0_P_4_b_2_coherent.fa |cut -d "|" -f 2- | sort -n > radtest
diff radref radtest
if [ $? -ne 0 ] ; then
       echo "*** With truncated bubbles FAILURE:"
       echo "*** Test: FAILURE on diff headers of .fa"
       exit 1
fi

rm radref radtest*



echo "****************************"
echo "***** ALL TESTS PASSED *****"
echo "****************************"


