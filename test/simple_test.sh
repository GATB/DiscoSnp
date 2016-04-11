#!/bin/bash

../run_discoSnp++.sh -r fof.txt -T

diff discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.fa
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff .fa"
  exit 1
fi

# we do not want to compare fildate lines since they are de facto not similar
awk '!/##filedate=/' discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf > discoRes_a.vcf
awk '!/##filedate=/' ref_discoRes_k_31_c_auto_D_100_P_1_b_0_coherent.vcf > discoRes_b.vcf
diff discoRes_a.vcf discoRes_b.vcf
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff .vcf"
  exit 1
fi

rm -f discoRes_*  reference_genome.fa.*

echo "*** Test: OK"
