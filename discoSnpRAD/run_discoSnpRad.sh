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
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
cyan=`tput setaf 6`
bold=`tput bold`
reset=`tput sgr0`


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
D=0 # maximal size of searched deletions
max_ambigous_indel=20
P=5 # number of polymorphsim per bubble
option_max_symmetrical_crossroads="5"
l="-l"
extend="-T"
x="-x"
e="-e"
output_coverage_option=""
genotyping="-genotype"
verbose=1
clustering="false"
short_read_connector_path=""
option_phase_variants=""
graph_reused="Egg62hdS7knSFvF3" # with -g option, we use a previously created graph. 


max_size_cluster=150
max_missing=0.95
min_rank=0.4

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
    echo $reset
    echo " ************"
    echo " *** HELP ***"
    echo " ************"
    echo "run_discoSnpRad.sh, pipelining kissnp2 and kissreads and clustering per locus for calling SNPs and small indels from RAD-seq data without the need of a reference genome"
    echo "Version: "$version
    echo "Cookbook: You may find a Cookbook here https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/COOKBOOK.md providing classical use cases."
    echo "Usage: ./run_discoSnpRad.sh --fof read_file_of_files --src [src_path] [OPTIONS]"
    echo "MANDATORY"
    echo "      -r|--fof <file name of a file of file(s)>"
    echo "           The input read files indicated in a file of file(s)"
    echo "           Example: -r bank.fof with bank.fof containing the two lines" 
    echo "             data_sample/reads_sequence1.fasta"
    echo "             data_sample/reads_sequence2.fasta.gz"
    echo "      Note: DiscoSnp-RAD uses files demultiplexed to samples. Each sample corresponds to exactly one line in the input file of files."
    echo ""
    echo "PARAMETERS"
    echo "      -k | --k_size value <int value>"
    echo "           Set the length of used kmers. Must fit the compiled value."
    echo "           Default=31"
    echo "      -c | --min_coverage value <int value>"
    echo "           Set the minimal coverage per read set: Used by kissnp2 (don't use kmers with lower coverage) and kissreads (read coherency threshold)." 
    echo "           This coverage can be automatically detected per read set (in this case use \"auto\" or specified per read set, see the documentation.)"
    echo "           Default=3"
    echo "      --high_precision | -R"
    echo "           lower recall / higher precision mode. With this parameter no symmetrical crossroads may be traversed during bubble detection (by default up to 5 symmetrical crossroads may be traversed during bubble detection)."

    echo ""
    echo "OPTIONS"
    echo "      -g | --graph <file name>"
    echo "           Reuse a previously created graph (.h5 file)"
    echo "      -p | --prefix <string>"
    echo "           All out files will start with this prefix. Default=\"discoRes\""
    echo "      -C | --max_coverage <int value>" 
    echo "           Set the maximal coverage for each read set: Used by kissnp2 (don't use kmers with higher coverage)."
    echo "           Default=2^31-1"
    echo "      -l | --no_low_complexity"
    echo "           Remove low complexity bubbles"
    echo "      -D | --deletion_max_size <int value>"
    echo "           discoSnpRad will search for deletions of size from 0 to D included. Default=0 (no deletion)"
    
    echo ""
    echo "CLUSTERING OPTION"
    echo "      -S|--src [src_path]"
    echo "           performs clustering of variants with short_read_connector"
    echo "           src_path: **absolute** path to short_read_connector directory, containing the \"short_read_connector.sh\" file. "
    echo "           -Note1: short read connector must be compiled."
    echo "           -Note2: if no value is given, it assumes short_read_connector.sh is in the PATH env variable."
    echo "           -Note3: with this option, discoSnpRad outputs a vcf file containing the variants clustered by locus" 
    echo "      --max_size_cluster <int value> "
    echo "           Discards cluster containing more than this number of variants. (Default 150)"
    echo "           Requires the -S or --src option"
    
    echo ""
    echo "FILTERING OPTION"
    echo "      --max_missing <float value> "
    echo "           Remove variants with more than max_missing % missing values. (Default 0.95, removes variants detected in 5% and less populations)"
    echo "      --min_rank <float value>"
    echo "           Remove variants whose rank is smaller than this threshold. (Default 0.4)"
    


    
    
#    echo "      -L | --max_diff_len <integer>" # Hidden - 
#    echo "           Longest accepted difference length between two paths of a truncated bubble"
#    echo "           default 0"
    

#    echo "      -a | --ambiguity_max_size <int>" # Hidden
#    echo "           Maximal size of ambiguity of INDELs. INDELS whose ambiguity is higher than this value are not output  [default '20']"
    echo ""
    echo "ADVANCED OPTIONS"
    echo "      -P | --max_snp_per_bubble <int>"
    echo "           discoSnpRad will search up to P SNPs in a unique bubble. Default=5"
    echo "      -d | --max_substitutions <int>"
    echo "           Set the number of authorized substitutions used while mapping reads on found SNPs (kissreads). Default=10"
    
    
    echo ""
    echo "MISC."
    echo "      -u | --max_threads <int>"
    echo "           Max number of used threads. 0 means all threads"
    echo "           default 0"
    echo "      -w      Wraith mode: only show all discoSnpRad commands without running them"
    echo "      -v <0 or 1>"
    echo "           Verbose 0 (avoids progress output) or 1 (enables progress output) -- default=1."
    echo "      -h | --help"
    echo "           Prints this message and exist"
    echo ""
    
    echo "Any further question: read the readme file or contact us via the Biostar forum: https://www.biostars.org/t/discosnp/"
}




#######################################################################
#################### GET OPTIONS                #######################
#######################################################################


echo "${yellow}"

while :; do
    case $1 in
    --max_missing)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
        max_missing=$2
        shift 1
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    --min_rank)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
        min_rank=$2
        shift 1
        echo "min_rank" ${min_rank}
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    --max_size_cluster)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
        max_size_cluster=$2
        shift 1
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;
    
    -A) 
        option_phase_variants="-phasing"
        echo "Will phase variants during kissreads process - WARNING this option is too experimental and thus not described in the help message"
        echo "You can obtain clusters using script : \"script/from_phased_alleles_to_clusters.sh file_name_of_phased_alleles\" (the filename(s) is/are given during kissreads process"
        ;;
    -w)
        wraith="true"
        ;;
    -S|--src)
    	clustering="true"
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
            short_read_connector_path=$2
            shift
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
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            graph_reused=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
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
echo $reset

#######################################################################
#################### END GET OPTIONS            #######################
#######################################################################

if [ -z "$read_sets" ]; then
    echo "${red}               **************************************************************************"
    echo "               ** ERROR: You must provide at least one read set (-r) "
    echo "               **************************************************************************"
    echo $reset
    exit 1
fi

#Checks if clustering can be performed

if [[ "$clustering" == "true" ]]; then
	# first tests the directory given by user if any
	if [ -n "$short_read_connector_path" ]; then
		src_file="$short_read_connector_path/short_read_connector.sh"
    	if [ -f "$src_file" ]; then
            echo "${yellow}short_read_connector path is $src_file$reset"
    	else
    		echo "${red}               **************************************************************************"
            echo "               ** WARNING: I cannot find short_read_connector (-S). "
            echo "               ** $src_file does not exist"
            echo "               ** I will not cluster variants per RAD locus"
            echo "               **************************************************************************"
            echo $reset
    		clustering="false"
    	fi
    else
    	#then tests if src is in the PATH env variable
    	src_file=$(command -v short_read_connector.sh)
    	if [ -n "$src_file" ]; then
    		echo "${yellow}short_read_connector path is $src_file$reset"
    	else
    		echo "${red}               **************************************************************************"
            echo "               ** WARNING: I cannot find short_read_connector in PATH. "
            echo "               ** Try giving the absolute path of short_read_connector directory with option -S"
            echo "               ** I will not cluster variants per RAD locus"
            echo "               **************************************************************************"
            echo $reset
    		clustering="false"
    	fi
    fi
fi
    		


######### CHECK THE k PARITY ##########
rest=$(( $k % 2 ))
if [ $rest -eq 0 ]
then
    echo "${red}# k=$k is even number, to avoid palindromes, we set it to $(($k-1))${reset}"
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
echo $green$cmdFofRemove "> ${read_sets}_${kissprefix}_removemeplease$cyan"
if [[ "$wraith" == "false" ]]; then
    $cmdFofRemove > ${read_sets}_${kissprefix}_removemeplease
fi


#######################################################################
#################### OPTIONS SUMMARY            #######################
#######################################################################

if [[ "$wraith" == "false" ]]; then
    echo "${yellow}     Running discoSnpRad "$version", in directory "$EDIR" with following parameters:"
    echo "           read_sets="$read_sets
    echo "           short_read_connector path="$short_read_connector_path
    echo "           prefix="$h5prefix
    if [ -f ${graph_reused} ]; then
        echo "           reuse graph="${graph_reused}
    fi
    echo "           c="$c
    echo "           C="$C
    echo "           k="$k
    echo "           b="$b
    echo "           d="$d
    echo "           D="$D
    echo "           max_truncated_path_length_difference="$max_truncated_path_length_difference
    echo -n "      starting date="
    date
    echo $reset
    echo
fi

#######################################################################
#################### END OPTIONS SUMMARY        #######################
#######################################################################

#############################################################
#################### DUMP READ FILES  #######################
#############################################################
dumpCmd="${read_file_names_bin} -in $read_sets"
    echo $green${dumpCmd} "> $readsFilesDump"$cyan
    if [[ "$wraith" == "false" ]]; then
        ${dumpCmd} > $readsFilesDump
    fi

if [ $? -ne 0 ]
then
    echo "${red}there was a problem with readFileName Dumping$reset"
    exit 1
fi


############################################################
#################### GRAPH CREATION  #######################
############################################################


if [ ! -f ${graph_reused} ]; then # no graph was given or the given graph was not a file. 
    T="$(date +%s)"
    echo "${yellow}     ############################################################"
    echo "     #################### GRAPH CREATION  #######################"
    echo "     ############################################################${reset}"

    graphCmd="${dbgh5_bin} -in ${read_sets}_${kissprefix}_removemeplease -out $h5prefix -kmer-size $k -abundance-min ${c_dbgh5} -abundance-max $C -solidity-kind one ${option_cores_gatb} -verbose $verbose  -skip-bcalm -skip-bglue -no-mphf"
    echo $green${graphCmd}$cyan
    if [[ "$wraith" == "false" ]]; then
        ${graphCmd}
    fi
    
    if [ $? -ne 0 ]
    then
        echo "${red}there was a problem with graph construction${reset}"
        exit 1
    fi

    T="$(($(date +%s)-T))"
    if [[ "$wraith" == "false" ]]; then
        echo "${yellow}Graph creation time in seconds: ${T}${reset}"
    fi
    graph_reused=$h5prefix.h5
else
    if [[ "$wraith" == "false" ]]; then
        echo "${yellow}File ${graph_reused} exists. We use it as input graph${reset}"
    fi
fi

cleanCmd="rm -rf trashme_*"
echo $green${cleanCmd}$cyan
if [[ "$wraith" == "false" ]]; then
    ${cleanCmd}
fi

######################################################
#################### KISSNP2   #######################
######################################################
T="$(date +%s)"
echo "${yellow}     ############################################################"
echo "     #################### KISSNP2 MODULE  #######################"
echo "     ############################################################${reset}"
kissnp2Cmd="${kissnp2_bin} -in ${graph_reused} -out ${kissprefix}_r  -b $b $l $x -P $P  -D $D $extend $option_cores_gatb $output_coverage_option -coverage_file ${h5prefix}_cov.h5 -max_ambigous_indel ${max_ambigous_indel} -max_symmetrical_crossroads ${option_max_symmetrical_crossroads}  -verbose $verbose -max_truncated_path_length_difference ${max_truncated_path_length_difference}"
echo $green${kissnp2Cmd}$cyan
if [[ "$wraith" == "false" ]]; then
    ${kissnp2Cmd}
fi
if [ $? -ne 0 ]
then
    echo "${red}there was a problem with kissnp2${reset}"
    exit 1
fi

T="$(($(date +%s)-T))"
if [[ "$wraith" == "false" ]]; then
    echo "${yellow}Bubble detection time in seconds: ${T}${reset}"
fi
if [ ! -f ${kissprefix}_r.fa ]
then
        if [[ "$wraith" == "false" ]]; then
        echo "${red}No polymorphism predicted by discoSnpRad"
        echo -n "      ending date="
        date
        echo "${yellow}      Thanks for using discoSnpRad - http://colibread.inria.fr/discoSnp/${reset}"
        exit 
    fi
fi

#######################################################################
#################### REDUNDANCY REMOVAL         #######################
#######################################################################
echo "${yellow}     ############################################################"
echo "     #################### REDUNDANCY REMOVAL  ###################"
echo "     ############################################################$reset"
redundancy_removal_cmd="python $EDIR/../scripts/redundancy_removal_discosnp.py ${kissprefix}_r.fa $k $kissprefix.fa"
echo $green${redundancy_removal_cmd}$cyan
if [[ "$wraith" == "false" ]]; then
   eval ${redundancy_removal_cmd}
fi
if [ $? -ne 0 ]
then
    echo "${red}there was a problem with redundancy removal$reset":
    exit 1
fi

#######################################################################
#################### KISSREADS                  #######################
#######################################################################

T="$(date +%s)"
echo "${yellow}     #############################################################"
echo "     #################### KISSREADS MODULE #######################"
echo "     #############################################################$reset"

smallk=$k
if (( $smallk>31 ))  ; then
    smallk=31
fi
i=5 #avoid modidy this (or increase this if memory needed by kissread is too high. Min 1. Large i (7-10) decreases memory and increases time).
index_stride=$(($i+1)); size_seed=$(($smallk-$i)) # DON'T modify this.

kissreadsCmd="${kissreads2_bin} -predictions $kissprefix.fa -reads  $read_sets -co ${kissprefix}_coherent -unco ${kissprefix}_uncoherent -k $k -size_seeds ${size_seed} -index_stride ${index_stride} -hamming $d  $genotyping -coverage_file ${h5prefix}_cov.h5 $option_cores_gatb  -verbose $verbose ${option_phase_variants}"

echo $green$kissreadsCmd$cyan
if [[ "$wraith" == "false" ]]; then
    eval $kissreadsCmd
fi
if [ $? -ne 0 ]
then
    echo "${red}there was a problem with kissreads2$reset":
    exit 1
fi
echo $reset
T="$(($(date +%s)-T))"
# echo "Kissreads (mapping reads on bubbles) time in seconds: ${T}"


#######################################################################
#################### SORT AND FORMAT  RESULTS #########################
#######################################################################

echo "${yellow}     ###############################################################"
echo "     #################### SORT AND FORMAT  RESULTS #################"
echo "     ###############################################################$reset"
if [[ "$wraith" == "false" ]]; then
    sort -rg ${kissprefix}_coherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_raw.fa
fi
if [ $? -ne 0 ]
then
    echo "${red}there was a problem with the result sorting.$reset"
    exit 1
fi
if [[ "$wraith" == "false" ]]; then
    sort -rg ${kissprefix}_uncoherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_uncoherent.fa
fi
if [ $? -ne 0 ]
then
    echo "${red}there was a problem with the result sorting$reset"
    exit 1
fi

rm -f ${read_sets}_${kissprefix}_removemeplease 
rm -f $kissprefix.fa ${kissprefix}_coherent ${kissprefix}_uncoherent
rm -f ${kissprefix}_r.fa
rm -f ${kissprefix}_uncoherent.fa

#######################################################################
#################### DISCOSNP FINISHED ###############################
#######################################################################



##################################################################################
#################### Deal with Downstream analyses ###############################
##################################################################################

echo "${yellow}     ###############################################################"
echo "     ######## CLUSTERING PER LOCUS AND/OR FORMATTING ###############"
echo "     ###############################################################$reset"

T="$(date +%s)"
if [[ "$clustering" == "true" ]]; then
    if [[ "$wraith" == "false" ]]; then
        echo "${yellow}Clustering and vcf formmatting$reset"
    fi
    final_output="${kissprefix}_clustered.vcf"
    cmd="$EDIR/clustering_scripts/discoRAD_clustering.sh -f ${kissprefix}_raw.fa -s $src_file -o ${final_output} -m ${max_missing} -r ${min_rank} -c ${max_size_cluster}"
    echo $green$cmd$cyan$reset
    if [[ "$wraith" == "false" ]]; then
        eval $cmd
    fi  
    T="$(($(date +%s)-T))"
    if [[ "$wraith" == "false" ]]; then
        echo "${yellow}RAD clustering per locus time in seconds: ${T}$reset"
    fi
    echo $reset
else
    if [[ "$wraith" == "false" ]]; then
        echo "${red}NO CLUSTERING (missing -S option)"
        echo "IF YOU WANT TO CLUSTERIZE RESULTS, RUN: "
        echo "  $EDIR/clustering_scripts/discoRAD_clustering.sh -f ${kissprefix}_raw.fa -s short_read_connector_path"
        #echo "  With short_read_connector_path indicating the directory containing short_read_connector.sh command "
        echo "Filtering and vcf formatting$reset"
    fi
    final_output="${kissprefix}.vcf"
    cmd="python3 $EDIR/../scripts/create_filtered_vcf.py -i ${kissprefix}_raw.fa -o ${final_output} -m ${max_missing} -r ${min_rank}"
    echo $green$cmd$cyan$reset
    if [[ "$wraith" == "false" ]]; then
        eval $cmd
    fi
    T="$(($(date +%s)-T))"
    if [[ "$wraith" == "false" ]]; then
        echo "${red}Filtering and vcf formatting time in seconds: ${T}$reset"
    fi
fi

if [[ "$wraith" == "false" ]]; then
    echo "${yellow}     ###############################################################"
    echo "     #################### DISCOSNPRAD FINISHED ######################"
    echo "     ###############################################################"
    Ttot="$(($(date +%s)-Ttot))"
    echo "DiscoSnpRad total time in seconds: ${Ttot}"
    echo "     ################################################################################################################"
    echo "      fasta of predicted variant is \""${kissprefix}_raw.fa"\""
    echo "      Ghost VCF file (1-based) is \""${final_output}"\""
    echo "      Thanks for using discoSnpRad - http://colibread.inria.fr/discoSnp/ - Forum: http://www.biostars.org/t/discoSnp/"
    echo "     ################################################################################################################${reset}"
fi


