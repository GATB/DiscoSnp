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



die() {
    printf '%s\n' "$1" >&2
    exit 1
}




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
b=2 # all bubbles accepted"
c=3 # minimal coverage
C=$max_C # maximal coverage
M=4
d=10 # estimated number of error per read (used by kissreads only)
D=3 # maximal size of searched deletions
max_ambigous_indel=20
P=5 # number of polymorphsim per bubble
option_max_symmetrical_crossroads="5"
l="-l"
extend="-T"
x="-x"
e="-e"
output_coverage_option=""
genotyping="-genotype"
remove=1
verbose=1
short_read_connector_path=""
#EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
EDIR=$( python -c "import os.path; print(os.path.dirname(os.path.realpath(\"${BASH_SOURCE[0]}\")))" ) # as suggested by Philippe Bordron 

if [ -d "$EDIR/../build/" ] ; then # VERSION SOURCE COMPILED
    read_file_names_bin=$EDIR/../build/bin/read_file_names
    dbgh5_bin=$EDIR/../build/ext/gatb-core/bin/dbgh5
    kissnp2_bin=$EDIR/../build/bin/kissnp2
    kissreads2_bin=$EDIR/../build/bin/kissreads2
else # VERSION BINARY
    read_file_names_bin=$EDIR/../bin/read_file_names
    dbgh5_bin=$EDIR/../bin/dbgh5
    kissnp2_bin=$EDIR/../bin/kissnp2
    kissreads2_bin=$EDIR/../bin/kissreads2
fi


chmod u+x $EDIR/../scripts/*.sh $EDIR/scripts/*.sh $EDIR/run_discoSnpRad.sh 2>/dev/null # Usefull for binary distributions

wraith="false"
max_truncated_path_length_difference=0
option_phase_variants=""
#######################################################################
#################### END HEADER                 #######################
#######################################################################

function help {
    echo " ************"
    echo " *** HELP ***"
    echo " ************"
    echo "run_discoSnpRad.sh, pipelining kissnp2 and kissreads and clustering per locus for calling SNPs and small indels from RAD-seq data without the need of a reference genome"
    echo "Version "$version
    echo "Usage: ./run_discoSnpRad.sh --fof read_file_of_files --src_path <directory> [OPTIONS]"
    echo -e "MANDATORY"
    echo -e "\t -r|--fof <file name of a file of file(s)>"
    echo -e "\t\t The input read files indicated in a file of file(s)"
    echo -e "\t\t Example: -r bank.fof with bank.fof containing the two lines \n\t\t\t data_sample/reads_sequence1.fasta\n\t\t\t data_sample/reads_sequence2.fasta.gz"
    
    echo -e "\nOPTIONS"
    echo -e "\t -S|--src_path <directory>"
    echo -e "\t\t **absolute** path to short_read_connector directory, containing the \"short_read_connector.sh\" file. "
    echo -e "\t\t -Note1: short read connector must be compiled."
    echo -e "\t\t -Note2: with this option, discoSnpRad provide a vcf file containing SNPs and INDELS, clustered by locus" 

    echo -e "\t -k | --k_size value <int value>"
    echo -e "\t\t Set the length of used kmers. Must fit the compiled value."
    echo -e "\t\t Default=31"
    echo -e "\t -c | --min_coverage value <int value>"
    echo -e "\t\t Set the minimal coverage per read set: Used by kissnp2 (don't use kmers with lower coverage) and kissreads (read coherency threshold)." 
    echo -e "\t\t This coverage can be automatically detected per read set (in this case use \"auto\" or specified per read set, see the documentation."
    echo -e "\t\t Default=3"
    echo -e "\t -C | --max_coverage value <int value in 0, 1 or 2>"
    echo -e "\t\t Set the maximal coverage for each read set: Used by kissnp2 (don't use kmers with higher coverage)."
    echo -e "\t\t Default=2^31-1"
    
    echo -e "\t --high_precision | -R"
    echo -e "\t\t low recall / high precision mode. With this parameter no symmetrical crossroads may be traversed during bubble detection (by default up to 5 symmetrical crossroads may be traversed during bubble detection)."
    
    echo -e "\t -L | --max_diff_len <integer>"
    echo -e "\t\t Longest accepted difference length between two paths of a truncated bubble"
    echo -e "\t\t default 0"
    
    echo -e "\t -g | --graph <file name>"
    echo -e "\t\t reuse a previously created graph (.h5 file) with same prefix and same k and c parameters."
    echo -e "\t -D | --deletion_max_size <int>"
    echo -e "\t\t discoSnpRad will search for deletions of size from 1 to D included. Default=100"
    echo -e "\t -a | --ambiguity_max_size <int>"
    echo -e "\t\t Maximal size of ambiguity of INDELs. INDELS whose ambiguity is higher than this value are not output  [default '20']"
    echo -e "\t -P | --max_snp_per_bubble <int>"
    echo -e "\t\t discoSnpRad will search up to P SNPs in a unique bubble. Default=3"
    echo -e "\t -p | --prefix <string>"
    echo -e "\t\t All out files will start with this prefix. Default=\"discoRes\""
    echo -e "\t -l | --no_low_complexity"
    echo -e "\t\t Remove low complexity bubbles"
    echo -e "\t -d | --max_substitutions <int>"
    echo -e "\t\t Set the number of authorized substitutions used while mapping reads on found SNPs (kissreads). Default=1"
    echo -e "\t -u | --max_threads <int>"
    echo -e "\t\t Max number of used threads. 0 means all threads"
    echo -e "\t\t default 0"

    
    echo -e "\t -w\t Wraith mode: only show all discoSnpRad commands without running them"
    echo -e "\t -v <0 or 1>"
    echo -e "\t\t Verbose 0 (avoids progress output) or 1 (enables progress output) -- default=1."
    echo -e "\t -h | --help"
    echo -e "\t\t Prints this message and exist\n"
    
    echo "Any further question: read the readme file or contact us via the Biostar forum: https://www.biostars.org/t/discosnp/"
}




#######################################################################
#################### GET OPTIONS                #######################
#######################################################################



while :; do
    case $1 in
    -w)
        wraith="true"
        ;;
    -S|--src_path)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
            short_read_connector_path=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -a|--ambiguity_max_size)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
            max_ambigous_indel=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;
       
    -v)        
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            verbose=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -R|--high_precision)
        option_max_symmetrical_crossroads=0
        ;; 

    -L | --max_diff_len)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            max_truncated_path_length_difference=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -g|--graph)
        remove=0
        ;;

    -l|--no_low_complexity)
        l=""
        ;;
        
    -h|-\?|--help)
        help
        exit 
        ;;

    -r|--fof)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            read_sets=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

        
    -p|--prefix)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            prefix=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -k | --k_size)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            k=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;


    -P|--max_snp_per_bubble)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            P=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -c|--min_coverage)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            c=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -C|--max_coverage)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            C=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -d|--max_substitutions)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            d=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -D|--deletion_max_size)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            D=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -x)
        x="-x" ##CHARLOTTE
        ;;

    -y)
        y="-x" ##CHARLOTTE
        ;;

    -G|--reference_genome)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            genome=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;


    -e)
        e="-e"
        ;;

    -u|--max_threads)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            option_cores_gatb="-nb-cores $2"
            option_cores_post_analysis="-t $2"
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -?*)
        printf 'WARN: Unknown option (exit): %s\n' "$1" >&2
        exit 1
        ;;

    :)
        echo "Option $1 requires an argument." >&2
        exit 1
        ;;
        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac

    shift
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
if [[ "$wraith" == "false" ]]; then
    echo $src_file
fi

if [[ "$wraith" == "false" ]]; then
    if [ -f "$src_file" ]; then
        if [[ "$wraith" == "false" ]]; then
            echo "short_read_connector is $src_file"
        fi
    else
        if [[ "$wraith" == "false" ]]; then
            echo -e "\t\t\t**************************************************************************"
            echo -e "\t\t\t** WARNING: I cannot find short_read_connector (-S). "
            echo -e "\t\t\t** $src_file does not exist"
            echo -e "\t\t\t** I will not cluster variants per RAD locus"
            echo -e "\t\t\t**************************************************************************"
        fi
    fi
fi 


######### CHECK THE k PARITY ##########
rest=$(( $k % 2 ))
if [ $rest -eq 0 ]
then
    echo "# k=$k is even number, to avoid palindromes, we set it to $(($k-1))"
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
kissprefix=${h5prefix}_D_${D}_P_${P}_m_${option_max_symmetrical_crossroads}
readsFilesDump=${prefix}_read_files_correspondance.txt


#######################################
c_dbgh5=$c
rm -f ${read_sets}_${kissprefix}_removemeplease
cmdFofRemove="cat ${read_sets}" > ${read_sets}_${kissprefix}_removemeplease
echo $cmdFofRemove "> ${read_sets}_${kissprefix}_removemeplease"
if [[ "$wraith" == "false" ]]; then
    $cmdFofRemove > ${read_sets}_${kissprefix}_removemeplease
fi


#######################################################################
#################### OPTIONS SUMMARY            #######################
#######################################################################

if [[ "$wraith" == "false" ]]; then
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
    echo -e "\t\t max_truncated_path_length_difference="$max_truncated_path_length_difference
    echo -e -n "\t starting date="
    date
    echo
fi

#######################################################################
#################### END OPTIONS SUMMARY        #######################
#######################################################################

#############################################################
#################### DUMP READ FILES  #######################
#############################################################
dumpCmd="${read_file_names_bin} -in $read_sets"
    echo ${dumpCmd} "> $readsFilesDump"
    if [[ "$wraith" == "false" ]]; then
        ${dumpCmd} > $readsFilesDump
    fi

if [ $? -ne 0 ]
then
    echo "there was a problem with readFileName Dumping":
    exit 1
fi


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
    if [[ "$wraith" == "false" ]]; then
        ${graphCmd}
    fi
    
    if [ $? -ne 0 ]
    then
        echo "there was a problem with graph construction"
        exit 1
    fi

    T="$(($(date +%s)-T))"
    if [[ "$wraith" == "false" ]]; then
        echo "Graph creation time in seconds: ${T}"
    fi
else
    if [[ "$wraith" == "false" ]]; then
        echo -e "File $h5prefix.h5 exists. We use it as input graph"
    fi
fi

cleanCmd="rm -rf trashme_*"
echo ${cleanCmd}
if [[ "$wraith" == "false" ]]; then
    ${cleanCmd}
fi

######################################################
#################### KISSNP2   #######################
######################################################
T="$(date +%s)"
echo -e "\t############################################################"
echo -e "\t#################### KISSNP2 MODULE  #######################"
echo -e "\t############################################################"
kissnp2Cmd="${kissnp2_bin} -in $h5prefix.h5 -out ${kissprefix}_r  -b $b $l $x -P $P  -D $D $extend $option_cores_gatb $output_coverage_option -coverage_file ${h5prefix}_cov.h5 -max_ambigous_indel ${max_ambigous_indel} -max_symmetrical_crossroads ${option_max_symmetrical_crossroads}  -verbose $verbose -max_truncated_path_length_difference ${max_truncated_path_length_difference}"
echo ${kissnp2Cmd}
if [[ "$wraith" == "false" ]]; then
    ${kissnp2Cmd}
fi
if [ $? -ne 0 ]
then
    echo "there was a problem with kissnp2"
    exit 1
fi

T="$(($(date +%s)-T))"
if [[ "$wraith" == "false" ]]; then
    echo "Bubble detection time in seconds: ${T}"
fi
if [ ! -f ${kissprefix}_r.fa ]
then
        if [[ "$wraith" == "false" ]]; then
        echo "No polymorphism predicted by discoSnpRad"
        echo -e -n "\t ending date="
        date
        echo -e "\t Thanks for using discoSnpRad - http://colibread.inria.fr/discoSnp/"
        exit 
    fi
fi

#######################################################################
#################### REDUNDANCY REMOVAL         #######################
#######################################################################
echo -e "\t############################################################"
echo -e "\t#################### REDUNDANCY REMOVAL  ###################"
echo -e "\t############################################################"
redundancy_removal_cmd="python $EDIR/../scripts/redundancy_removal_discosnp.py ${kissprefix}_r.fa $k $kissprefix.fa"
echo ${redundancy_removal_cmd}
if [[ "$wraith" == "false" ]]; then
    ${redundancy_removal_cmd}
fi
if [ $? -ne 0 ]
then
    echo "there was a problem with redundancy removal":
    exit 1
fi

#######################################################################
#################### KISSREADS                  #######################
#######################################################################

T="$(date +%s)"
echo -e "\t#############################################################"
echo -e "\t#################### KISSREADS MODULE #######################"
echo -e "\t#############################################################"

smallk=$k
if (( $smallk>31 ))  ; then
    smallk=31
fi
i=5 #avoid modidy this (or increase this if memory needed by kissread is too high. Min 1. Large i (7-10) decreases memory and increases time).
index_stride=$(($i+1)); size_seed=$(($smallk-$i)) # DON'T modify this.

kissreadsCmd="${kissreads2_bin} -predictions $kissprefix.fa -reads  $read_sets -co ${kissprefix}_coherent -unco ${kissprefix}_uncoherent -k $k -size_seeds ${size_seed} -index_stride ${index_stride} -hamming $d  $genotyping -coverage_file ${h5prefix}_cov.h5 $option_cores_gatb  -verbose $verbose ${option_phase_variants}"

echo $kissreadsCmd
if [[ "$wraith" == "false" ]]; then
    eval $kissreadsCmd
fi
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
if [[ "$wraith" == "false" ]]; then
    sort -rg ${kissprefix}_coherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_coherent_raw.fa
fi
if [ $? -ne 0 ]
then
    echo "there was a problem with the result sorting."
    exit 1
fi
if [[ "$wraith" == "false" ]]; then
    sort -rg ${kissprefix}_uncoherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_uncoherent.fa
fi
if [ $? -ne 0 ]
then
    echo "there was a problem with the result sorting"
    exit 1
fi

rm -f ${read_sets}_${kissprefix}_removemeplease 
rm -f $kissprefix.fa ${kissprefix}_coherent ${kissprefix}_uncoherent

#######################################################################
#################### DISCOSNP FINISHED ###############################
#######################################################################






##################################################################################
#################### Deal with Downstream analyses ###############################
##################################################################################

echo -e "\t###############################################################"
echo -e "\t######## CLUSTERING PER LOCUS AND/OR FORMATTING ###############"
echo -e "\t###############################################################"

T="$(date +%s)"
if [ -f "$src_file" ]; then
    if [[ "$wraith" == "false" ]]; then
        echo "Clustering and vcf formmatting"
    fi
    final_output="${kissprefix}_coherent_clustered.vcf"
    cmd="$EDIR/scripts/discoRAD_clustering.sh -f ${kissprefix}_coherent_raw.fa -s $short_read_connector_path -o ${final_output}"
    echo $cmd
    if [[ "$wraith" == "false" ]]; then
        eval $cmd
    fi  
    T="$(($(date +%s)-T))"
    if [[ "$wraith" == "false" ]]; then
        echo "RAD clustering per locus time in seconds: ${T}"
    fi
else
    if [[ "$wraith" == "false" ]]; then
        echo "NO CLUSTERING (missing -s option)"
        echo "IF YOU WANT TO CLUSTERIZE RESULTS, RUN: "
        echo "  $EDIR/scripts/discoRAD_clustering.sh -f ${kissprefix}_coherent_raw.fa -s short_read_connector_path"
        #echo "  With short_read_connector_path indicating the directory containing short_read_connector.sh command "
        echo "Filtering and vcf formatting"
    fi
    final_output="${kissprefix}_coherent.vcf"
    cmd="python3 $EDIR/../scripts/create_filtered_vcf.py -i ${kissprefix}_coherent_raw.fa -o ${final_output} -m 0.95 -r 0.4"
    echo $cmd
    if [[ "$wraith" == "false" ]]; then
        eval $cmd
    fi
    T="$(($(date +%s)-T))"
    if [[ "$wraith" == "false" ]]; then
        echo "Filtering and vcf formatting time in seconds: ${T}"
    fi
fi

if [[ "$wraith" == "false" ]]; then
    echo -e "\t###############################################################"
    echo -e "\t#################### DISCOSNPRAD FINISHED ######################"
    echo -e "\t###############################################################"
    Ttot="$(($(date +%s)-Ttot))"
    echo "DiscoSnpRad total time in seconds: ${Ttot}"
    echo -e "\t################################################################################################################"
    echo -e "\t fasta of predicted variant is \""${kissprefix}_coherent_raw.fa"\""
    echo -e "\t Ghost VCF file (1-based) is \""${final_output}"\""
    echo -e "\t Thanks for using discoSnpRad - http://colibread.inria.fr/discoSnp/ - Forum: http://www.biostars.org/t/discoSnp/"
    echo -e "\t################################################################################################################"
fi


