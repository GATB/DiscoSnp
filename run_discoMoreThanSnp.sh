#!/bin/bash

###########################################################
#################### DEFAULT VALUES #######################
###########################################################
version="2.0.0"
read_sets="" # FOR instance: "read_set1.fa.gz read_set2.fq.gz"
prefix="discoRes" # all temp and final files will be written will start with this prefix
k=31 # size of kmers
b=0 # smart branching approach: bubbles in which both paths are equaly branching are  discarded, all others are accepted
c=4 # minimal coverage
d=1 # estimated number of error per read (used by kissreads only)
D=0 # maximal size of searched deletions
#######################################################################
#################### END HEADER                 #######################
#######################################################################

function help {
echo "run_discoSnp.sh, a pipelining kissnp2 and kissreads for calling SNPs from NGS reads without the need of a reference genome"
echo "Version "$version
echo "Usage: ./run_discoSnp.sh OPT"
echo -e "\tMANDATORY:"
echo -e "\t \t -r list of reads separated by space, surrounded by the '\"' character. Note that reads may be in fasta or fastq format, gzipped or not."
echo -e "\t \t    Example: -r \"data_sample/reads_sequence1.fasta   data_sample/reads_sequence2.fasta.gz\"."
echo -e "\tOPT:"
echo -e "\t\t -b value. "
echo -e "\t\t\t 0: forbid SNPs for wich any of the two paths is branching (high precision, lowers the recal in complex genomes). Default value"
echo -e "\t\t\t 1: (smart branching) forbid SNPs for wich the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position"
echo -e "\t\t\t 2: No limitation on branching (lowers the precision, high recall)"
echo -e "\t\t -D value. If specified, discoXXX will search for deletions of size from 1 to D included. Default=0"
echo -e "\t\t -p prefix. All out files will start with this prefix. Default=\"discoRes\""
echo -e "\t\t -k value. Set the length of used kmers. Must fit the compiled value. Default=31"
echo -e "\t\t -c value. Set the minimal coverage: Used by kissnp2 (don't use kmers with lower coverage) and kissreads (read coherency threshold). Default=4"
echo -e "\t\t -d value. Set the number of authorized substitutions used while mapping reads on found SNPs (kissreads). Default=1"
echo -e "\t\t -h: Prints this message and exist"
echo "Any further question: read the readme file or contact us: pierre.peterlongo@inria.fr"
}


#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts ":r:p:k:c:d:D:b:h" opt; do
case $opt in

	w)
	remove=0
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

c)
echo "use c=$OPTARG" >&2
c=$OPTARG
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
	


######### CHECK THE k PARITY ##########
rest=$(( $k % 2 ))
if [ $rest -eq 0 ]
then
echo "k=$k is even number, to avoid palindromes, we set it to $(($k-1))"
k=$(($k-1))
fi
#######################################
prefix=$prefix\_k_$k\_c_$c

#######################################################################
#################### OPTIONS SUMMARY            #######################
#######################################################################
MY_PATH="`( cd \"$MY_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$MY_PATH" ] ; then
# error; for some reason, the path is not accessible
# to the script (e.g. permissions re-evaled after suid)
exit 1  # fail
fi
echo -e "\tRunning discoSnp "$version", in directory "$MY_PATH" with following parameters:"
echo -e "\t\t read_sets="$read_sets
echo -e "\t\t prefix="$prefix
echo -e "\t\t c="$c
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





if [ ! -e $prefix.h5 ]; then
	echo -e "\t############################################################"
	echo -e "\t#################### GRAPH CREATION  #######################"
	echo -e "\t############################################################"
	echo "./ext/gatb-core/bin/dbgh5 -in `echo $read_sets | tr " " ","` -out $prefix -kmer-size $k -abundance $c "
	./ext/gatb-core/bin/dbgh5 -in `echo $read_sets | tr " " ","` -out $prefix -kmer-size $k -abundance $c 
	if [ $? -ne 0 ]
	then
		echo "there was a problem with kissnp2, command line: ./ext/gatb-core/bin/dbgh5 -in `echo $read_sets | tr " " ","` -out $prefix -kmer-size $k -abundance $c  "
		exit
	fi

else
	echo -e "File $prefix.h5 exists. We use it as input graph"
fi
	

echo -e "\t############################################################"
echo -e "\t#################### KISSNP2 MODULE  #######################"
echo -e "\t############################################################"

./tools/kissnp2/kissnp2 -D $D -in $prefix.h5 -out $prefix\_D_$D -T -b $b 

if [ $? -ne 0 ]
then
    echo "there was a problem with kissnp2, command line: ./tools/kissnp2/kissnp2 -D $D -in $prefix.h5 -out $prefix -T -b $b "
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


./tools/kissreads/kissreads $prefix\_D_$D.fa $read_sets -k $smallk -i $i -O $k -c $c -d $d -n  -o $prefix\_D_$D\_coherent -u $prefix\_D_$D\_uncoherent 
if [ $? -ne 0 ]
then
echo "there was a problem with kissnp2, command line: ./tools/kissreads/kissreads $prefix\_D_$D.fa $read_sets -k $smallk -i $i -O $k -c $c -d $d -n  -o $prefix\_D_$D\_coherent -u $prefix\_D_$D\_uncoherent" 
exit
fi



echo -e "\t###############################################################"
echo -e "\t#################### DISCOSNPs FINISHED #######################"
echo -e "\t###############################################################"
#######################################################################
#################### SORT AND FORMAT COHERENT RESULTS #################
#######################################################################
sort -rg $prefix\_D_$D\_coherent | cut -d " " -f 2 | tr ';' '\n' > $prefix\_D_$D\_coherent.fa
if [ $? -ne 0 ]
then
echo "there was a problem with the result sorting, command line: sort -rg $prefix\_D_$D\_coherent | cut -d " " -f 2 | tr ';' '\n' > $prefix\_D_$D\_coherent.fa"
exit
fi



echo -e -n "\t ending date="
date
echo -e "\t SNPs are stored in \""$prefix\_D_$D\_coherent.fa"\""
echo -e "\t Thanks for using discoSnp - http://colibread.inria.fr/discoSnp/"



