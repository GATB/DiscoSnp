#!/bin/bash



for fof in fof2.txt fof3.txt fof4.txt ; do 
       for b in 0 1 2 ; do
              fof_base=${fof%.*}
              echo testing fof=$fof_base b=${b}
              ../../run_discoSnp++.sh -r ${fof} -P 2 -c 1  -b ${b} -p test_${fof_base} > log 2>log2
       
              
              
              # TEST THE FASTA IF EXISTS
              A=golden_standard_results/ref_${fof_base}_k_31_c_1_D_100_P_2_b_${b}_coherent.fa
              B=test_${fof_base}_k_31_c_1_D_100_P_2_b_${b}_coherent.fa
              if [ -f $A ]
              then
                     DIFF=$(diff $A $B)
                     if [ "$DIFF" != "" ] 
                     then
                            echo "FASTA Results differ with $b on $fof"
                            echo "try diff $B $B"
                            exit
                     fi
              else
                     if [ -f $b ]
                     then 
                            echo "FASTA Result should not exist with $b on $fof"
                            exit
                     fi
              fi
              echo "FASTA OK with ${fof_base} ${b}"
              

              # TEST THE VCF IF EXISTS
              A=golden_standard_results/ref_${fof_base}_k_31_c_1_D_100_P_2_b_${b}_coherent.vcf
              B=test_${fof_base}_k_31_c_1_D_100_P_2_b_${b}_coherent.vcf
              if [ -f $A ]
              then
                     grep -v "#" $A > AFACREF
                     grep -v "#" $B > AFACTEST
                     DIFF=$(diff AFACREF AFACTEST)
                     if [ "$DIFF" != "" ] 
                     then
                            echo "VCF Results differ with $b on $fof"
                            echo "try diff AFACREF AFACTEST"
                            exit
                     fi
              else
                     if [ -f $b ]
                     then 
                            echo "VCF Result should not exist with $b on $fof"
                            exit
                     fi
              fi
              echo "VCF OK with $fof_base ${b}"
              echo 
              echo
              
              rm -f $B AFACTEST AFACREF test_*
       done
done



