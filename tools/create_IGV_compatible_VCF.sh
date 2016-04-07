#!/bin/bash


function help {
echo " ##############################"
echo "   Create a IGV compatible file"
echo " ##############################"
echo "Usage : ./create_IGV_compatible_VCF.sh VCF_file"

echo -e "\t VCF_file: the vcf created by run_VCF_creator"
echo -e "\t 1/ Remove from this file the unmapped variants"
echo -e "\t 2/ Make the vcf 0-based"
echo -e "\t 3/ sort variants by position"
echo -e "\t-h: print this message"
}
if test -z "$1" 
then
       help
       exit
fi
vcffile=$1
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


igvfiletemp=$(basename $vcffile .vcf)"_for_IGV.vcf.tmp"
igvfile=$(basename $vcffile .vcf)"_for_IGV.vcf"
cat $vcffile|grep "#">$igvfile
#cat $vcffile|grep -v  "#"|sort -k 2n,2n -n|grep -v "^SNP"|grep -v "^INDEL">>$igvfile
cat $vcffile|grep -v  "#"|sort -k 1,1 -k 2,2n |grep -v "^SNP"|grep -v "^INDEL">>$igvfile # from 2 2 6
#python $DIR/tools/one2zeroBased_vcf.py $igvfiletemp 
#cat VCFone2zeroBAsed.vcf >> $igvfile
#rm -f $igvfiletemp VCFone2zeroBAsed.vcf
echo -e "... Creation of the vcf file for IGV: done ...==> $igvfile"



