#!/bin/bash
#*****************************************************************************
#   discoSnpRad: discovering polymorphism from raw unassembled RADSEQ NGS reads
#   A tool from the GATB (Genome Assembly Tool Box)
#   Copyright (C) 2017  INRIA
#   Authors: J. Gauthier, C. Mouden, C. Lemaitre, P. Peterlongo
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

# cmd="$EDIR/run_discoSnpRad.sh "$@" -x -t -e -c 3 -b 1"
# echo "I run discoSnpRad with following command line: " ${cmd}
# ${cmd}

#TODO: call automatically the clustering +radseq filters. 






EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
function myrealpath { echo $(cd $(dirname $1); pwd)/$(basename $1); }

option_cores_gatb=""
option_cores_post_analysis=""

Ttot="$(date +%s)"
#### constant #####
max_C=2147483647 #$((2**31-1))

###########################################################
#################### DEFAULT VALUES #######################
###########################################################
version="2.3.X"
read_sets="" # A file of file(s)
prefix="discoRad" # all intermediate and final files will be written will start with this prefix
k=31 # size of kmers
b=1 # smart branching approach: bubbles in which both paths are equaly branching are  discarded, all others are accepted
c=3 # minimal coverage
C=$max_C # maximal coverage
M=4
d=1 # estimated number of error per read (used by kissreads only)
D=100 # maximal size of searched deletions
max_ambigous_indel=20
P=5 # number of polymorphsim per bubble
option_max_symmetrical_crossroads=""
l="-l"
extend="-t"
x="-x"
e="-e"
output_coverage_option=""
genotyping="-genotype"
remove=1
verbose=1
short_read_connector_path=""
EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [ -d "$EDIR/build/" ] ; then # VERSION SOURCE COMPILED
    read_file_names_bin=$EDIR/build/bin/read_file_names
    dbgh5_bin=$EDIR/build/ext/gatb-core/bin/dbgh5
    kissnp2_bin=$EDIR/build/bin/kissnp2
    kissreads2_bin=$EDIR/build/bin/kissreads2
else # VERSION BINARY
    read_file_names_bin=$EDIR/bin/read_file_names
    dbgh5_bin=$EDIR/bin/dbgh5
    kissnp2_bin=$EDIR/bin/kissnp2
    kissreads2_bin=$EDIR/bin/kissreads2
fi


chmod +x $EDIR/scripts/*.sh $EDIR/scripts_RAD/*.sh $EDIR/run_discoSnpRad.sh 2>/dev/null # Usefull for binary distributions

useref=""
genome=""
bwa_path_option=""
bwa_distance=4

#######################################################################
#################### END HEADER                 #######################
#######################################################################

function help {
    echo ""
    echo ""
    echo ""
    echo "run_discoSnpRad.sh, pipelining kissnp2 and kissreads and clustering per locus for calling SNPs and small indels from RAD-seq data without the need of a reference genome"
    echo "Version "$version
    echo "Usage: ./run_discoSnpRad.sh -r read_file_of_files [OPTIONS]"
    echo -e "\tMANDATORY:"
    echo -e "\t\t -r read_file_of_files"
    echo -e "\t\t    Example: -r bank.fof with bank.fof containing the two lines \n\t\t\t data_sample/reads_sequence1.fasta\n\t\t\t data_sample/reads_sequence2.fasta.gz"
    echo -e "\t\t -S **absolute** path to short_read_connector directory, containing the ``short_read_connector.sh'' file. "
    echo -e "\t\t\t -Note1: short read connector must be compiled."
    echo -e "\t\t\t -Note2: if this option is missing, discoSnpRad will still however provide a fasta file containing SNPs and INDELS, that won't be clustered by locus" 
    echo -e "\tDISCOSNPRAD OPTIONS:"
    echo -e "\t\t -g: reuse a previously created graph (.h5 file) with same prefix and same k and c parameters."
    echo -e "\t\t -b value. "
    echo -e "\t\t\t 1: (smart branching) forbid SNPs for which the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position Default value"
    echo -e "\t\t\t 2: No limitation on branching (lowers the precision, high recall)"
    echo -e "\t\t -s value. In b2 mode only: maximal number of symmetrical croasroads traversed while trying to close a bubble. Default: no limit"
    echo -e "\t\t -D value. discoSnpRad will search for deletions of size from 1 to D included. Default=100"
    echo -e "\t\t -a value. Maximal size of ambiguity of INDELs. INDELS whose ambiguity is higher than this value are not output  [default '20']"
    echo -e "\t\t -P value. discoSnpRad will search up to P SNPs in a unique bubble. Default=5"
    echo -e "\t\t -p prefix. All out files will start with this prefix. Default=\"discoRad\""
    echo -e "\t\t -l: remove low complexity bubbles"
    echo -e "\t\t -k value. Set the length of used kmers. Must fit the compiled value. Default=31"
    echo -e "\t\t -c value. Set the minimal coverage per read set: Used by kissnp2 (don't use kmers with lower coverage) and kissreads (read coherency threshold). This coverage can be automatically detected per read set or specified per read set, see the documentation. Default=auto"
    echo -e "\t\t -C value. Set the maximal coverage for each read set: Used by kissnp2 (don't use kmers with higher coverage). Default=2^31-1"
    echo -e "\t\t -d value. Set the number of authorized substitutions used while mapping reads on found SNPs (kissreads). Default=1"
    echo -e "\t\t -u: max number of used threads"
    echo -e "\t\t -v: verbose 0 (avoids progress output) or 1 (enables progress output) -- default=1."
    # echo -e "\t\t -x: variant detection radseq optimization" #CHARLOTTE -- NOW CALLED BY DEFAULT USING SCRIP run_discoSnpRad
    # echo -e "\t\t -y: variant coverage radseq optimization" #CHARLOTTE -- NOW CALLED BY DEFAULT USING SCRIP run_discoSnpRad


    echo 
    echo -e "\t\t -h: Prints this message and exist"
    # echo -e "\t\t -e: map SNP predictions on reference genome with their extensions. - Forced usage when using discoSnpRad"
    echo "Any further question: read the readme file or contact us via the Biostar forum: https://www.biostars.org/t/discosnp/"
}


#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts ":r:p:k:c:C:d:D:b:s:P:S:htTlgu:a:v:" opt; do
    case $opt in
 
    S)
        short_read_connector_path=$OPTARG
        ;;
    a)
        max_ambigous_indel=$OPTARG
        ;;
       
    v)
        verbose=$OPTARG
        ;;

    s)
        option_max_symmetrical_crossroads="-max_symmetrical_crossroads "$OPTARG
        echo ${option_max_symmetrical_crossroads}
        ;;
    g)
        remove=0
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
        # TODO B0 FORBIDEN
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
    echo -e "\t\t\t**************************************************************************"
    echo -e "\t\t\t** ERROR: You must provide at least one read set (-r) "
    echo -e "\t\t\t**************************************************************************"
    exit 1
fi

src_file="$short_read_connector_path/short_read_connector.sh"


if [ -f "$src_file" ]; then
    echo "short_read_connector is $src_file"
else
    echo -e "\t\t\t**************************************************************************"
    echo -e "\t\t\t** WARNING: I cannot find short_read_connector (-S). "
    echo -e "\t\t\t** $src_file does not exist"
    echo -e "\t\t\t** I will not cluster variants per RAD locus"
    echo -e "\t\t\t**************************************************************************"
fi


######### CHECK THE k PARITY ##########
rest=$(( $k % 2 ))
if [ $rest -eq 0 ]
then
    echo "k=$k is even number, to avoid palindromes, we set it to $(($k-1))"
    k=$(($k-1))
fi


#######################################
c_filename=`echo ${c} | tr ',' '_'`
if [ $C -ne $max_C ]
then
    h5prefix=${prefix}_k_${k}_c_${c_filename}_C_${C}
else
    h5prefix=${prefix}_k_${k}_c_${c_filename}

fi
kissprefix=${h5prefix}_D_${D}_P_${P}_b_${b}
readsFilesDump=${prefix}_read_files_correspondance.txt


#######################################
c_dbgh5=$c
rm -f ${read_sets}_${kissprefix}_removemeplease
if [[ "$useref" == "true" ]]; then

    if [ -z "$genome" ]; then
        echo "You can't use option -R without providing a reference genome (-G)"
        help
        exit 1
    fi
       
    myrealpath $genome > ${read_sets}_${kissprefix}_removemeplease
    c_dbgh5="1,"$c

fi
cat $read_sets >> ${read_sets}_${kissprefix}_removemeplease



#######################################################################
#################### OPTIONS SUMMARY            #######################
#######################################################################

echo -e "\tRunning discoSnpRad "$version", in directory "$EDIR" with following parameters:"
echo -e "\t\t read_sets="$read_sets
echo -e "\t\t short_read_connector path="$short_read_connector_path
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
${read_file_names_bin} -in $read_sets > $readsFilesDump



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

    graphCmd="${dbgh5_bin} -in ${read_sets}_${kissprefix}_removemeplease -out $h5prefix -kmer-size $k -abundance-min ${c_dbgh5} -abundance-max $C -solidity-kind one ${option_cores_gatb} -verbose $verbose  -skip-bcalm -skip-bglue -no-mphf"
    echo ${graphCmd}
    ${graphCmd}

    if [ $? -ne 0 ]
    then
        echo "there was a problem with graph construction"
        exit 1
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
kissnp2Cmd="${kissnp2_bin} -in $h5prefix.h5 -out $kissprefix  -b $b $l $x -P $P  -D $D $extend $option_cores_gatb $output_coverage_option -coverage_file ${h5prefix}_cov.h5 -max_ambigous_indel ${max_ambigous_indel} ${option_max_symmetrical_crossroads}  -verbose $verbose"
echo ${kissnp2Cmd}
${kissnp2Cmd}

if [ $? -ne 0 ]
then
    echo "there was a problem with kissnp2"
    exit 1
fi

T="$(($(date +%s)-T))"
echo "Bubble detection time in seconds: ${T}"

if [ ! -f $kissprefix.fa ]
then
    echo "No polymorphism predicted by discoSnpRad"
    echo -e -n "\t ending date="
    date
    echo -e "\t Thanks for using discoSnpRad - http://colibread.inria.fr/discoSnp/"
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

smallk=$k
if (( $smallk>31 ))  ; then
    smallk=31
fi
i=5 #avoid modidy this (or increase this if memory needed by kissread is too high. Min 1. Large i (7-10) decreases memory and increases time).
index_stride=$(($i+1)); size_seed=$(($smallk-$i)) # DON'T modify this.

kissreadsCmd="${kissreads2_bin} -predictions $kissprefix.fa -reads  $read_sets -co ${kissprefix}_coherent -unco ${kissprefix}_uncoherent -k $k -size_seeds ${size_seed} -index_stride ${index_stride} -hamming $d  $genotyping -coverage_file ${h5prefix}_cov.h5 $option_cores_gatb  -verbose $verbose"

echo $kissreadsCmd
$kissreadsCmd

if [ $? -ne 0 ]
then
    echo "there was a problem with kissreads2":
    exit 1
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
    exit 1
fi

sort -rg ${kissprefix}_uncoherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_uncoherent.fa
if [ $? -ne 0 ]
then
    echo "there was a problem with the result sorting"
    exit 1
fi

#rm -f $kissprefix.fa ${kissprefix}_coherent ${kissprefix}_uncoherent
rm -f ${read_sets}_${kissprefix}_removemeplease

#######################################################################
#################### DISCOSNP FINISHED ###############################
#######################################################################






##################################################################################
#################### Deal with Downstream analyses ###############################
##################################################################################

echo -e "\t###############################################################"
echo -e "\t#################### CLUSTERING PER LOCUS #####################"
echo -e "\t###############################################################"

T="$(date +%s)"
if [ -f "$src_file" ]; then
    cmd="$EDIR/scripts_RAD/discoRAD_finalization.sh ${kissprefix}_coherent.fa $short_read_connector_path" 
    echo $cmd
    $cmd
    T="$(($(date +%s)-T))"
    echo "RAD clustering per locus time in seconds: ${T}"
else
    echo "IF YOU WANT TO CLUSTERIZE RESULTS, RUN: "
    echo "  $EDIR/scripts_RAD/discoRAD_finalization.sh ${kissprefix}_coherent.fa short_read_connector_path"
    echo "  With short_read_connector_path indicating the directory containing short_read_connector.sh command "
fi


echo -e "\t###############################################################"
echo -e "\t#################### DISCOSNPRAD FINISHED ######################"
echo -e "\t###############################################################"
Ttot="$(($(date +%s)-Ttot))"
echo "DiscoSnpRad total time in seconds: ${Ttot}"
echo -e "\t################################################################################################################"
echo -e "\t fasta of predicted variant is \""${kissprefix}_coherent.fa"\""
echo -e "\t Ghost VCF file (1-based) is \""${kissprefix}_coherent_sorted_with_clusters.vcf"\""
echo -e "\t Thanks for using discoSnpRad - http://colibread.inria.fr/discoSnp/ - Forum: http://www.biostars.org/t/discoSnp/"
echo -e "\t################################################################################################################"



