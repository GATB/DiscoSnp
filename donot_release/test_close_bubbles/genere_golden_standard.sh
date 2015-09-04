#!/bin/bash



for fof in fof2.txt fof3.txt fof4.txt ; do 
       for b in 0 1 2 ; do
              fof_base=${fof%.*}
              ../../run_discoSnp++.sh -r ${fof} -P 2 -c 1  -b ${b} -p ref_${fof_base} > log 2>log2
       done
done
              
           
rm -f golden_standard_results/*
mv ref_* golden_standard_results

