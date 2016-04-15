# Get and decompress the data
wget http://www.irisa.fr/symbiose/people/ppeterlongo/data_test_disco.zip 
unzip data_test_disco.zip

# Create the file of file: 
ls humch1_00* > fof.txt

#####################
# Default option run: 
#####################
../../run_discoSnp++.sh -r fof.txt -T
# Test the .fa 
# The sequence ids and orders are not conserved due to parallelisation. This explains why we separate sequences from headers and why we remove ids
grep -v ">" ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa | sort > ref
grep -v ">" discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa | sort > created
diff ref created
if [ $? -ne 0 ] ; then
       echo "*** Default option FAILURE:"
  echo "*** Test: FAILURE on diff sequences of .fa"
  exit 1
fi

grep ">" ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa |cut -d "|" -f 2- | sort > ref
grep ">" discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa |cut -d "|" -f 2- | sort > created
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
awk '!/#/' ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort > ref
awk '!/#/' discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort > created
#diff created ref
#if [ $? -ne 0 ] ; then
#       echo "*** Default option FAILURE:"
#  echo "*** Test: FAILURE on diff content of .vcf"
#  exit 1
#fi



######################################################
# With reference run (using previously created graph): 
######################################################
../../run_discoSnp++.sh -r fof.txt -T -G humch1_first_5M.fasta -g

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
awk '!/#/'  ref_with_mapping_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf  | cut -f 2,4- | sort > ref
awk '!/#/'  discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf | cut -f 2,4- | sort > created
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
awk '!/#/' ref_with_mapping_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort > ref
awk '!/#/' discoRes_k_31_c_auto_D_100_P_1_b_0_coherent_for_IGV.vcf | cut -f 2,4- | sort > created
diff created ref
if [ $? -ne 0 ] ; then
       echo "*** With mapping on ref FAILURE:"
  echo "*** Test: FAILURE on diff content of for IGV .vcf"
  exit 1
fi


rm -f created ref discoRes*
rm -f data_test_disco.zip
rm -f humch1_*
rm -f fof.txt



