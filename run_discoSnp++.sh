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

#!/bin/bash

#### constant #####
max_C=$((2**31-1))

###########################################################
#################### DEFAULT VALUES #######################
###########################################################
version="2.0.7"
read_sets="" # FOR instance: "read_set1.fa.gz read_set2.fq.gz"
prefix="discoRes" # all intermediate and final files will be written will start with this prefix
k=31 # size of kmers
b=0 # smart branching approach: bubbles in which both paths are equaly branching are  discarded, all others are accepted
c=4 # minimal coverage
C=$max_C # maximal coverage
d=1 # estimated number of error per read (used by kissreads only)
D=0 # maximal size of searched deletions
P=1 # number of polymorphsim per bubble
l="-l"
extend=""
genotyping="-g"
paired=""
remove=1
EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DISCO_BUILD_PATH="$EDIR/build/"


#######################################################################
#################### END HEADER                 #######################
#######################################################################

function help {
echo "run_discoSnp++.sh, a pipelining kissnp2 and kissreads for calling SNPs and small indels from NGS reads without the need of a reference genome"
echo "Version "$version
echo "Usage: ./run_discoSnp++.sh OPT"
echo -e "\tMANDATORY:"
echo -e "\t\t -r list of reads separated by space, surrounded by the '\"' character. Note that reads may be in fasta or fastq format, gzipped or not."
echo -e "\t\t    Example: -r \"data_sample/reads_sequence1.fasta   data_sample/reads_sequence2.fasta.gz\"."

echo -e "\tOPT:"
echo -e "\t\t -m: indicates that read sets are paired. Each couple of read sets is considered as paired (set1_1.fa set1_2.fa set2_1.fa set2_2.fa ...) "
echo -e "\t\t -g: reuse a previously created graph (.h5 file) with same prefix and same k and c parameters."
echo -e "\t\t -b value. "
echo -e "\t\t\t 0: forbid variants for which any of the two paths is branching (high precision, lowers the recall in complex genomes). Default value"
echo -e "\t\t\t 1: (smart branching) forbid SNPs for which the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position"
echo -e "\t\t\t 2: No limitation on branching (lowers the precision, high recall)"
echo -e "\t\t -D value. discoSnp++ will search for deletions of size from 1 to D included. Default=0"
echo -e "\t\t -P value. discoSnp++ will search up to P SNPs in a unique bubble. Default=1"
echo -e "\t\t -p prefix. All out files will start with this prefix. Default=\"discoRes\""
echo -e "\t\t -l: remove low complexity bubbles"
echo -e "\t\t -k value. Set the length of used kmers. Must fit the compiled value. Default=31"
echo -e "\t\t -t: extend found polymorphisms with unitigs"
echo -e "\t\t -T: extend found polymorphisms with contigs"
echo -e "\t\t -c value. Set the minimal coverage per read set: Used by kissnp2 (don't use kmers with lower coverage) and kissreads (read coherency threshold). Default=4"
echo -e "\t\t -C value. Set the maximal coverage per read set: Used by kissnp2 (don't use kmers with higher coverage). Default=2^31-1"
echo -e "\t\t -d value. Set the number of authorized substitutions used while mapping reads on found SNPs (kissreads). Default=1"
echo -e "\t\t -n: do not compute the genotypes"
echo -e "\t\t -h: Prints this message and exist"
echo "Any further question: read the readme file or contact us: pierre.peterlongo@inria.fr"
}


#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts ":r:p:k:c:C:d:D:b:P:htTlmgn" opt; do
case $opt in
	t)
	extend="-t"
	;;
	
	T)
	extend="-T"
	;;

	g)
	remove=0
	;;
	
	m)
	paired="-P"
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

######### CHECK THE number of readsets parity if the -P option is requiered ###########
if [[ $paired == "-P" ]]
then
	nb_read_sets=`echo $read_sets | wc -w`
	rest=$(( $nb_read_sets % 2 ))
	if [ $rest -eq 1 ]
	then
	echo "[ERROR] You cannot ask for paired read sets with an odd number of read files ($nb_read_sets)"
	exit
	fi
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
	h5prefix=$prefix\_k_$k\_c_$c\_C_$C
else
	h5prefix=$prefix\_k_$k\_c_$c
	
fi
kissprefix=$h5prefix\_D_$D\_P_$P\_b_$b

if [[ $l == "-l" ]]
then
	kissprefix=$kissprefix\_withlow
fi


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
if [[ $paired == "-P" ]]
then
	echo -e "\t\t Reads are paired"
fi
echo -e -n "\t starting date="
date
echo

#######################################################################
#################### END OPTIONS SUMMARY        #######################
#######################################################################



if [ $remove -eq 1 ]; then
	rm -f $h5prefix.h5
fi

if [ ! -e $h5prefix.h5 ]; then
	echo -e "\t############################################################"
	echo -e "\t#################### GRAPH CREATION  #######################"
	echo -e "\t############################################################"
	echo "$DISCO_BUILD_PATH/ext/gatb-core/bin/dbgh5 -in `echo $read_sets | tr " " ","` -out $h5prefix -kmer-size $k -abundance-min $c -abundance-max $C -solidity-kind max"
	$DISCO_BUILD_PATH/ext/gatb-core/bin/dbgh5 -in `echo $read_sets | tr " " ","` -out $h5prefix -kmer-size $k -abundance-min $c -abundance-max $C -solidity-kind max
	if [ $? -ne 0 ]
	then
		echo "there was a problem with graph construction, command line: $DISCO_BUILD_PATH/ext/gatb-core/bin/dbgh5 -in `echo $read_sets | tr " " ","` -out $h5prefix -kmer-size $k -abundance-min $c -abundance-max $C -solidity-kind max"
		exit
	fi

else
	echo -e "File $h5prefix.h5 exists. We use it as input graph"
fi
	

echo -e "\t############################################################"
echo -e "\t#################### KISSNP2 MODULE  #######################"
echo -e "\t############################################################"
echo "$DISCO_BUILD_PATH/tools/kissnp2/kissnp2 -D $D -in $h5prefix.h5 -out $kissprefix  -b $b $l -P $P $extend"
$DISCO_BUILD_PATH/tools/kissnp2/kissnp2 -D $D -in $h5prefix.h5 -out $kissprefix  -b $b $l -P $P $extend

if [ $? -ne 0 ]
then
    echo "there was a problem with kissnp2, command line: $DISCO_BUILD_PATH/tools/kissnp2/kissnp2 -D $D -in $h5prefix.h5 -out $kissprefix  -b $b $l -P $P $extend"
    exit
fi


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
echo
echo -e "\t#############################################################"
echo -e "\t#################### KISSREADS MODULE #######################"
echo -e "\t#############################################################"

i=4 #avoid modidy this (or increase this if memory needed by kissread is too high. Min 1. Large i (7-10) decreases memory and increases time).
smallk=$(($k-$i-1)) # DON'T modify this.


echo "$DISCO_BUILD_PATH/tools/kissreads/kissreads $kissprefix.fa $read_sets -k $smallk -i $i -O $k -c $c -d $d -n $genotyping -o $kissprefix\_coherent -u $kissprefix\_uncoherent $paired"


$DISCO_BUILD_PATH/tools/kissreads/kissreads $kissprefix.fa $read_sets -k $smallk -i $i -O $k -c $c -d $d -n $genotyping -o $kissprefix\_coherent -u $kissprefix\_uncoherent $paired
if [ $? -ne 0 ]
then
echo "there was a problem with kissnp2, command line: $DISCO_BUILD_PATH/tools/kissreads/kissreads $kissprefix.fa $read_sets -k $smallk -i $i -O $k -c $c -d $d -n $genotyping -o $kissprefix\_coherent -u $kissprefix\_uncoherent $paired"
exit
fi

echo -e "\t###############################################################"
echo -e "\t#################### SORT AND FORMAT  RESULTS #################"
echo -e "\t###############################################################"
#######################################################################
#################### SORT AND FORMAT  RESULTS #########################
#######################################################################
echo "sort -rg $kissprefix\_coherent | cut -d " " -f 2 | tr ';' '\n' > $kissprefix\_coherent.fa"
sort -rg $kissprefix\_coherent | cut -d " " -f 2 | tr ';' '\n' > $kissprefix\_coherent.fa
if [ $? -ne 0 ]
then
echo "there was a problem with the result sorting, command line: sort -rg $kissprefix\_coherent | cut -d " " -f 2 | tr ';' '\n' > $kissprefix\_coherent.fa"
exit
fi

echo "sort -rg $kissprefix\_uncoherent | cut -d " " -f 2 | tr ';' '\n' > $kissprefix\_uncoherent.fa"
sort -rg $kissprefix\_uncoherent | cut -d " " -f 2 | tr ';' '\n' > $kissprefix\_uncoherent.fa
if [ $? -ne 0 ]
then
echo "there was a problem with the result sorting, command line: sort -rg $kissprefix\_uncoherent | cut -d " " -f 2 | tr ';' '\n' > $kissprefix\_uncoherent.fa"
exit
fi

rm -f $kissprefix.fa $kissprefix\_coherent $kissprefix\_uncoherent


#######################################################################
#################### DISCOSNP FINISHED ###############################
#######################################################################

echo -e "\t###############################################################"
echo -e "\t#################### DISCOSNP++ FINISHED #######################"
echo -e "\t###############################################################"
echo -e -n "\t ending date="
date
echo -e "\t SNPs and indels are stored in \""$kissprefix\_coherent.fa"\""
echo -e "\t Thanks for using discoSnp++ - http://colibread.inria.fr/discoSnp/ - Forum: http://www.biostars.org/t/discoSnp/"



