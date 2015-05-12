#! /bin/bash

################################################################################
# we download some banks
################################################################################
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR039/ERR039477/ERR039477.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR387/SRR387476/SRR387476.fastq.gz"


################################################################################
# we launch DiscoSNP
################################################################################
./run_discoSnp++.sh -r "./ERR039477.fastq.gz ./SRR387476.fastq.gz" -u 1


################################################################################
# we check the result
################################################################################
echo "eb49fdefddd8627eeecc427b9b4ce7c5  discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.fa"   > reference1
echo "3fd37576255360df6206debebb6a3f16  discoRes_k_31_c_4_D_0_P_1_b_0_withlow_uncoherent.fa" > reference2
echo "ead4f4d85b1596648a6e35364de62c69  discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.vcf"  > reference3

md5sum discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.fa   > check1
md5sum discoRes_k_31_c_4_D_0_P_1_b_0_withlow_uncoherent.fa > check2
md5sum discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.vcf  > check3

diff ./reference1 ./check1
if [ $? -eq 0 ]; then
   echo "TEST 1 OK"
else
   echo "TEST 1 KO"
fi

diff ./reference2 ./check2
if [ $? -eq 0 ]; then
   echo "TEST 2 OK"
else
   echo "TEST 2 KO"
fi

diff ./reference3 ./check3
if [ $? -eq 0 ]; then
   echo "TEST 3 OK"
else
   echo "TEST 3 KO"
fi


################################################################################
# clean up
################################################################################
rm -f  ERR039477.fastq.gz  SRR387476.fastq.gz
rm -f discoRes_k_31_c_4.h5  
rm -f discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.fa  discoRes_k_31_c_4_D_0_P_1_b_0_withlow_uncoherent.fa   discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.vcf
rm -f reference1 reference2 reference3  check1 check2 check3
