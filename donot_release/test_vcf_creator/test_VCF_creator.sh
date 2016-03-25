#!/bin/bash
# usage : ./test_VCF_creator.sh /home/Workhan_Data/Documents/Programmes/bwa-0.7.10/ >logTest.txt
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && cd .. && pwd )
echo -e " Test VCF_Creator will take : "
echo -e "\t#HMMM_with_indels_k_31_c_4_D_500_P_5_b_1_withlow_coherent.fa"
echo -e "\t#GenomePseudomonasNC_007492.fasta"
echo -e "\t#data_test_creator.sam"

pathBWA=$1
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
echo -e "****************Running Test 1****************"
$DIR/run_VCF_creator.sh -G GenomePseudomonasNC_007492.fasta -p HMMM_with_indels_k_31_c_4_D_500_P_5_b_1_withlow_coherent.fa -o afac.vcf -B $pathBWA -I
if ! diff afac.vcf GOLDstandardHMMPseudomonas.vcf;then
        echo "!!! Test 1 : Difference between the two files !!!"
fi

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
echo -e "****************Running Test 2****************"
$DIR/run_VCF_creator.sh -f data_test_creator.sam -o afac_VCF_data_test.vcf 
if  ! diff VCF_data_test.vcf GoldVCF_data_test.vcf ;then
        echo "!!! Test 2 : Difference between the two files!!!"
fi

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
echo -e "****************Running Test 3****************"
$DIR/run_VCF_creator.sh -p HMMM_with_indels_k_31_c_4_D_500_P_5_b_1_withlow_coherent.fa -o afac_withoutref.vcf
if ! diff afac_withoutref.vcf GOLDstandardHMMPseudomonasWithoutRef.vcf;then
        echo "!!! Test 1 : Difference between the two files !!!"
fi

