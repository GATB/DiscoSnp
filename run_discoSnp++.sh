#!/bin/bash
#*****************************************************************************
#   discoSnp++: discovering polymorphism from raw unassembled NGS reads
#   A tool from the GATB (Genome Assembly Tool Box)
#   Copyright (C) 2014  INRIA
#   Authors: P.Peterlongo, E.Drezen
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************

option_cores_gatb=""
option_cores_post_analysis=""

Ttot="$(date +%s)" 
#### constant #####
max_C=2147483647 #$((2**31-1))

###########################################################
#################### DEFAULT VALUES #######################
###########################################################
version="2.2.1"
read_sets="" # A file of file(s)
prefix="discoRes" # all intermediate and final files will be written will start with this prefix
k=31 # size of kmers
b=0 # smart branching approach: bubbles in which both paths are equaly branching are  discarded, all others are accepted
c=auto # minimal coverage
C=$max_C # maximal coverage
M=4
d=1 # estimated number of error per read (used by kissreads only)
D=100 # maximal size of searched deletions
max_ambigous_indel=20
P=1 # number of polymorphsim per bubble
l="-l"
extend=""
output_coverage_option=""
genotyping="-genotype"
paired=""
remove=1
EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DISCO_BUILD_PATH="$EDIR/build/"





useref=""
genome="" 
bwa_path_option=""
bwa_distance=4

#######################################################################
#################### END HEADER                 #######################
#######################################################################

function help {
echo "run_discoSnp++.sh, a pipelining kissnp2 and kissreads for calling SNPs and small indels from NGS reads without the need of a reference genome"
echo "Version "$version
echo "Usage: ./run_discoSnp++.sh -r read_file_of_files [OPTIONS]"
echo -e "\tMANDATORY:"
echo -e "\t\t -r read_file_of_files" 
echo -e "\t\t    Example: -r bank.fof with bank.fof containing the two lines \n\t\t\t data_sample/reads_sequence1.fasta\n\t\t\t data_sample/reads_sequence2.fasta.gz" 

echo -e "\tDISCOSNP++ OPTIONS:"
echo -e "\t\t -g: reuse a previously created graph (.h5 file) with same prefix and same k and c parameters."
echo -e "\t\t -b value. "
echo -e "\t\t\t 0: forbid variants for which any of the two paths is branching (high precision, lowers the recall in complex genomes). Default value"
echo -e "\t\t\t 1: (smart branching) forbid SNPs for which the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position"
echo -e "\t\t\t 2: No limitation on branching (lowers the precision, high recall)"
echo -e "\t\t -D value. discoSnp++ will search for deletions of size from 1 to D included. Default=100"
echo -e "\t\t -a value. Maximal size of ambiguity of INDELs. INDELS whose ambiguity is higher than this value are not output  [default '20']"
echo -e "\t\t -P value. discoSnp++ will search up to P SNPs in a unique bubble. Default=1"
echo -e "\t\t -p prefix. All out files will start with this prefix. Default=\"discoRes\""
echo -e "\t\t -l: remove low complexity bubbles"
echo -e "\t\t -k value. Set the length of used kmers. Must fit the compiled value. Default=31"
echo -e "\t\t -t: extend found polymorphisms with unitigs"
echo -e "\t\t -T: extend found polymorphisms with contigs"
echo -e "\t\t -c value. Set the minimal coverage per read set: Used by kissnp2 (don't use kmers with lower coverage) and kissreads (read coherency threshold). This coverage can be automatically detected per read set or specified per read set, see the documentation. Default=auto"
echo -e "\t\t -C value. Set the maximal coverage for each read set: Used by kissnp2 (don't use kmers with higher coverage). Default=2^31-1"
echo -e "\t\t -d value. Set the number of authorized substitutions used while mapping reads on found SNPs (kissreads). Default=1"
echo -e "\t\t -n: do not compute the genotypes"
echo -e "\t\t -u: max number of used threads"


echo -e "\t REFERENCE GENOME AND/OR VCF CREATION OPTIONS"

echo -e "\t\t -G: reference genome file (fasta, fastq, gzipped or nor). In absence of this file the VCF created by VCF_creator won't contain mapping related results."
echo -e "\t\t -R: use the reference file also in the variant calling, not only for mapping results"
echo -e "\t\t -B: bwa path. e.g. /home/me/my_programs/bwa-0.7.12/ (note that bwa must be pre-compiled)"
echo -e "\t\t\t Optional unless option -G used and bwa is not in the binary path."
echo -e "\t\t -M: Maximal number of mapping errors during BWA mapping phase."
echo -e "\t\t\t Useless unless mapping on reference genome is required (option -G). Default=4. "
echo -e "\t\t -h: Prints this message and exist"
echo "Any further question: read the readme file or contact us via the Biostar forum: https://www.biostars.org/t/discosnp/"
}


#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts ":r:p:k:c:C:d:D:b:P:htTlRmgnG:B:M:u:a:" opt; do
case $opt in
       R)
       useref="true"
       output_coverage_option="-dont_output_first_coverage"
       ;;
       
       a)
       max_ambigous_indel=$OPTARG
       ;;
       
	t)
	extend="-t"
	;;
	
	T)
	extend="-T"
	;;

	g)
	remove=0
	;;
	
	
	n)
	genotyping=""
	;;
	
	l)
	l=""
	;;
h)
help
exit
;;

r)
echo "use read set: $OPTARG" >&2
read_sets=$OPTARG
;;


b)
echo "use branching strategy: $OPTARG" >&2
b=$OPTARG
;;

p)
echo "use prefix=$OPTARG" >&2
prefix=$OPTARG
;;

k)
echo "use k=$OPTARG" >&2
k=$OPTARG
;;


P)
echo "use P=$OPTARG" >&2
P=$OPTARG
;;

c)
echo "use c=$OPTARG" >&2
c=$OPTARG
;;

C)
echo "use C=$OPTARG" >&2
C=$OPTARG
;;

d)
echo "use d=$OPTARG" >&2
d=$OPTARG
;;

D)
echo "use D=$OPTARG" >&2
D=$OPTARG
;;

	B)
	echo -e "BWA directory: $OPTARG" >&2
	bwa_path_option="-B "$OPTARG
	;;


	G)
	echo -e "use genome : $OPTARG" >&2
	genome=$OPTARG
	;;
       
       M)
       echo "use M=$OPTARG" >&2
       M=$OPTARG
       ;;
       
       u)
       echo "use at most $OPTARG cores" >&2
       option_cores_gatb="-nb-cores $OPTARG"
       option_cores_post_analysis="-t $OPTARG"
       ;;
       
\?)
echo "Invalid option: -$OPTARG" >&2
exit 1
;;

:)
echo "Option -$OPTARG requires an argument." >&2
exit 1
;;
esac
done
#######################################################################
#################### END GET OPTIONS            #######################
#######################################################################

if [ -z "$read_sets" ]; then
	echo "You must provide at least one read set (-r)"
help
exit
fi





if [ -d "$DISCO_BUILD_PATH" ] ; then
echo "Binaries in $DISCO_BUILD_PATH"
else
ls "$DISCO_BUILD_PATH"
echo "error; for some reason, the discoSnp++ build path ($DISCO_BUILD_PATH) is not accessible"
echo "Please indicate (option -B) the location of the discoSnp++ build directory"
exit 1  # fail
fi


######### CHECK THE k PARITY ##########
rest=$(( $k % 2 ))
if [ $rest -eq 0 ]
then
echo "k=$k is even number, to avoid palindromes, we set it to $(($k-1))"
k=$(($k-1))
fi


#######################################
if [ $C -ne $max_C ]
then
	h5prefix=${prefix}_k_${k}_c_${c}_C_${C}
else
	h5prefix=${prefix}_k_${k}_c_${c}
	
fi
kissprefix=${h5prefix}_D_${D}_P_${P}_b_${b}
readsFilesDump=${prefix}_read_files_correspondance.txt


#######################################
c_dbgh5=$c
if [[ "$useref" == "true" ]]; then

       if [ -z "$genome" ]; then
              echo "You can't use option -R without providing a reference genome (-G)"
              help
              exit
       fi
              
       echo $genome > ${read_sets}_${kissprefix}_removemeplease
       c_dbgh5="1,"$c
       echo $c_dbgh5
       cat $read_sets >> ${read_sets}_${kissprefix}_removemeplease
fi
cat $read_sets > ${read_sets}_${kissprefix}_removemeplease



#######################################################################
#################### OPTIONS SUMMARY            #######################
#######################################################################
MY_PATH="`( cd \"$MY_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$MY_PATH" ] ; then
# error; for some reason, the path is not accessible
# to the script (e.g. permissions re-evaled after suid)
exit 1  # fail
fi
echo -e "\tRunning discoSnp++ "$version", in directory "$MY_PATH" with following parameters:"
echo -e "\t\t read_sets="$read_sets
echo -e "\t\t prefix="$h5prefix
echo -e "\t\t c="$c
echo -e "\t\t C="$C
echo -e "\t\t k="$k
echo -e "\t\t b="$b
echo -e "\t\t d="$d
echo -e "\t\t D="$D

echo -e -n "\t starting date="
date
echo

#######################################################################
#################### END OPTIONS SUMMARY        #######################
#######################################################################

#############################################################
#################### DUMP READ FILES  #######################
#############################################################
$DISCO_BUILD_PATH/tools/read_file_names/read_file_names -in $read_sets > $readsFilesDump



############################################################
#################### GRAPH CREATION  #######################
############################################################
if [ $remove -eq 1 ]; then
	rm -f $h5prefix.h5
fi

if [ ! -e $h5prefix.h5 ]; then
       T="$(date +%s)"      
	echo -e "\t############################################################"
	echo -e "\t#################### GRAPH CREATION  #######################"
	echo -e "\t############################################################"
	#echo $DISCO_BUILD_PATH/ext/gatb-core/bin/dbgh5 -in ${read_sets}_${kissprefix}_removemeplease -out $h5prefix -kmer-size $k -abundance-min $c_dbgh5 -abundance-max $C -solidity-kind one $option_cores_gatb
	#$DISCO_BUILD_PATH/ext/gatb-core/bin/dbgh5 -in ${read_sets}_${kissprefix}_removemeplease -out $h5prefix -kmer-size $k -abundance-min $c_dbgh5 -abundance-max $C -solidity-kind one $option_cores_gatb
       
       graphCmd="$DISCO_BUILD_PATH/ext/gatb-core/bin/dbgh5 -in ${read_sets}_${kissprefix}_removemeplease -out $h5prefix -kmer-size $k -abundance-min ${c_dbgh5} -abundance-max $C -solidity-kind one ${option_cores_gatb}"
       echo ${graphCmd}
       ${graphCmd}
       
	if [ $? -ne 0 ]
	then
		echo "there was a problem with graph construction"
		exit
	fi

       T="$(($(date +%s)-T))"
       echo "Graph creation time in seconds: ${T}"

else
	echo -e "File $h5prefix.h5 exists. We use it as input graph"
fi
	

######################################################
#################### KISSNP2   #######################
######################################################
T="$(date +%s)"      
echo -e "\t############################################################"
echo -e "\t#################### KISSNP2 MODULE  #######################"
echo -e "\t############################################################"
kissnp2Cmd="$DISCO_BUILD_PATH/tools/kissnp2/kissnp2 -in $h5prefix.h5 -out $kissprefix  -b $b $l -P $P  -D $D $extend $option_cores_gatb $output_coverage_option -coverage_file ${h5prefix}_cov.h5 -max_ambigous_indel ${max_ambigous_indel}"
echo ${kissnp2Cmd}
${kissnp2Cmd}

if [ $? -ne 0 ]
then
    echo "there was a problem with kissnp2"
    exit
fi

T="$(($(date +%s)-T))"
echo "Bubble detection time in seconds: ${T}"

if [ ! -f $kissprefix.fa ]
then
 echo "No polymorphism predicted by discoSnp++"
echo -e -n "\t ending date="
date
echo -e "\t Thanks for using discoSnp++ - http://colibread.inria.fr/discoSnp/"
 exit
fi

#######################################################################
#################### KISSREADS                  #######################
#######################################################################

T="$(date +%s)" 
echo
echo -e "\t#############################################################"
echo -e "\t#################### KISSREADS MODULE #######################"
echo -e "\t#############################################################"


if (( $smallk>31 ))  ; then
  smallk=31
fi
i=5 #avoid modidy this (or increase this if memory needed by kissread is too high. Min 1. Large i (7-10) decreases memory and increases time).
index_stride=$(($i+1)); size_seed=$(($smallk-$i)) # DON'T modify this.

kissreadsCmd="$DISCO_BUILD_PATH/tools/kissreads2/kissreads2 -predictions $kissprefix.fa -reads  $read_sets -co ${kissprefix}_coherent -unco ${kissprefix}_uncoherent -k $k -size_seeds ${size_seed} -index_stride ${index_stride} -hamming $d  $genotyping -coverage_file ${h5prefix}_cov.h5 $option_cores_gatb"

echo $kissreadsCmd
$kissreadsCmd

if [ $? -ne 0 ]
then
echo "there was a problem with kissreads2": 
exit
fi

T="$(($(date +%s)-T))"
# echo "Kissreads (mapping reads on bubbles) time in seconds: ${T}"


#######################################################################
#################### SORT AND FORMAT  RESULTS #########################
#######################################################################

echo -e "\t###############################################################"
echo -e "\t#################### SORT AND FORMAT  RESULTS #################"
echo -e "\t###############################################################"
sort -rg ${kissprefix}_coherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_coherent.fa
if [ $? -ne 0 ]
then
echo "there was a problem with the result sorting."
exit
fi

sort -rg ${kissprefix}_uncoherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_uncoherent.fa
if [ $? -ne 0 ]
then
echo "there was a problem with the result sorting"
exit
fi

rm -f $kissprefix.fa ${kissprefix}_coherent ${kissprefix}_uncoherent

#######################################################################
#################### DISCOSNP FINISHED ###############################
#######################################################################






#######################################################################
#################### Deal with VCF ###############################
#######################################################################

T="$(date +%s)" 
echo -e "\t###############################################################"
echo -e "\t#################### CREATE VCF         #######################"
echo -e "\t###############################################################"

if [ -z "$genome" ]; then #  NO reference genome use, vcf creator mode 1
       vcfCreatorCmd="$EDIR/run_VCF_creator.sh -p ${kissprefix}_coherent.fa -o ${kissprefix}_coherent.vcf"
       echo $vcfCreatorCmd
       $vcfCreatorCmd
       if [ $? -ne 0 ]
       then
       echo "there was a problem with VCF creation. See how to use the \"run_VCF_creator.sh\" alone."
       fi
else # A Reference genome is provided, vcf creator mode 2
       vcfCreatorCmd="$EDIR/run_VCF_creator.sh $bwa_path_option -G $genome $bwa_path_option -p ${kissprefix}_coherent.fa -o ${kissprefix}_coherent.vcf  -I $option_cores_post_analysis"
       echo $vcfCreatorCmd
       $vcfCreatorCmd

       if [ $? -ne 0 ]
       then
       echo "there was a problem with VCF creation. See how to use the \"run_VCF_creator.sh\" alone."
       fi
fi

T="$(($(date +%s)-T))"
echo "Vcf creation time in seconds: ${T}"

echo -e "\t###############################################################"
echo -e "\t#################### DISCOSNP++ FINISHED ######################"
echo -e "\t###############################################################"
Ttot="$(($(date +%s)-Ttot))"
echo "DiscoSnp++ total time in seconds: ${Ttot}"
echo -e "\t################################################################################################################"
echo -e "\t fasta of predicted variant is \""${kissprefix}_coherent.fa"\""

if [ -z "$genome" ]; then
       echo -e "\t Ghost VCF file (1-based) is \""${kissprefix}_coherent.vcf"\""
else
       echo -e "\t VCF file (1-based) is \""${kissprefix}_coherent.vcf"\""
       echo -e "\t An IGV ready VCF file (sorted by position, only mapped variants, 0-based) is \""${kissprefix}_coherent_for_IGV.vcf"\""
fi
echo -e "\t Thanks for using discoSnp++ - http://colibread.inria.fr/discoSnp/ - Forum: http://www.biostars.org/t/discoSnp/"
echo -e "\t################################################################################################################"
