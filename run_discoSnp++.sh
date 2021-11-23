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
read_sets_kissreads=""
prefix="discoRes" # all intermediate and final files will be written will start with this prefix
k=31 # size of kmers
b=0 # smart branching approach: bubbles in which both paths are equaly branching are  discarded, all others are accepted
c=3 # minimal coverage
C=$max_C # maximal coverage
d=1 # estimated number of error per read (used by kissreads only)
D=100 # maximal size of searched deletions
max_ambigous_indel=20
P=3 # number of polymorphsim per bubble
option_max_symmetrical_crossroads=""
l="-l"
extend="-t"
x="" # if set to -x : authorize non closed bubbles (used mainly in RNA seq)
output_coverage_option=""
genotyping="-genotype"
verbose=1
stop_after_kissnp=0
aav=0 # if set to not zero : considers the aav context: 1/ no uncoherent 2/ authorize -x option 
e="" # if set to -e: Map variant predictions on reference genome with their unitig or contig extensions
graph_reused="Egg62hdS7knSFvF3" # with -g option, we use a previously created graph. 

#EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
EDIR=$( python -c "import os.path; print(os.path.dirname(os.path.realpath(\"${BASH_SOURCE[0]}\")))" ) # as suggested by Philippe Bordron 


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

chmod +x $EDIR/scripts/*.sh $EDIR/run_discoSnpRad.sh 2>/dev/null # Usefull for binary distributions

useref=""
wraith="false"
genome=""
bwa_path_option=""
option_phase_variants=""
bwa_distance=4

#######################################################################
#################### END HEADER                 #######################
#######################################################################

function help {
    echo " ************"
    echo " *** HELP ***"
    echo " ************"
    echo "run_discoSnp++.sh, a pipelining kissnp2 and kissreads for calling SNPs and small indels from NGS reads without the need of a reference genome"
    echo "Version "$version
    echo "Usage: ./run_discoSnp++.sh -r read_file_of_files [OPTIONS]"
    echo -e "MANDATORY"
    echo -e "\t -r|--fof <file name of a file of file(s)>"
    echo -e "\t\t The input read files indicated in a file of file(s)"
    echo -e "\t\t Example: -r bank.fof with bank.fof containing the two lines \n\t\t\t data_sample/reads_sequence1.fasta\n\t\t\t data_sample/reads_sequence2.fasta.gz"

    echo -e "\nOPTIONS"
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
    echo -e "\t -b | --branching value. "
    echo -e "\t\t 0: forbid variants for which any of the two paths is branching (high precision, lowers the recall in complex genomes)."
    echo -e "\t\t Default value"
    echo -e "\t\t 1: (smart branching) forbid SNPs for which the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position"
    echo -e "\t\t2: No limitation on branching (lowers the precision, high recall)"
    echo -e "\t -s | --symmetrical value <int value>"
    echo -e "\t\t In -b 2 mode only: maximal number of symmetrical crossroads traversed while trying to close a bubble. Default: no limit"
    echo -e "\t -g | --graph <file name>"
    echo -e "\t\t Reuse a previously created graph (.h5 file)"
    echo -e "\t -X\t Stop discoSnp++ right after variant calling - the output is only a fasta file with no coverage information."
    echo -e "\t -D | --deletion_max_size <int>"
    echo -e "\t\t discoSnp++ will search for deletions of size from 0 to D included. -D 0 means no deletion detected. Default=100."
    echo -e "\t -a | --ambiguity_max_size <int>"
    echo -e "\t\t Maximal size of ambiguity of INDELs. INDELS whose ambiguity is higher than this value are not output  [default '20']"
    echo -e "\t -P | --max_snp_per_bubble <int>"
    echo -e "\t\t discoSnp++ will search up to P SNPs in a unique bubble. Default=3"
    echo -e "\t --fof_mapping <file name of a file of file(s)>"
    echo -e "\t\t If this option is used this fof is used when mapping back reads on the predicted variants instead of the original fof file provided by -r|--fof option"     
    echo -e "\t -p | --prefix <string>"
    echo -e "\t\t All out files will start with this prefix. Default=\"discoRes\""
    echo -e "\t -l | --no_low_complexity"
    echo -e "\t\t Remove low complexity bubbles"
    echo -e "\t -T | --contigs"
    echo -e "\t\t Extend found polymorphisms with contigs (default: extend with unitigs)"
    echo -e "\t -d | --max_substitutions <int>"
    echo -e "\t\t Set the number of authorized substitutions used while mapping reads on found SNPs (kissreads). Default=1"
    echo -e "\t -n | --no_genotype"
    echo -e "\t\t Do not compute the genotypes"
    echo -e "\t -u | --max_threads <int>"
    echo -e "\t\t Max number of used threads. 0 means all threads"
    echo -e "\t\t default 0"


    echo -e "\nREFERENCE GENOME AND/OR VCF CREATION OPTIONS"
    echo -e "\t -G | --reference_genome <file name>"
    echo -e "\t\t Reference genome file (fasta, fastq, gzipped or nor). In absence of this file the VCF created by VCF_creator won't contain mapping related results."
    echo -e "\t -R"
    echo -e "\t\t Use the reference file also in the variant calling, not only for mapping results"
    echo -e "\t -B | --bwa_path <directory name>"
    echo -e "\t\t bwa path. e.g. /home/me/my_programs/bwa-0.7.12/ (note that bwa must be pre-compiled)"
    echo -e "\t\t Optional unless option -G used and bwa is not in the binary path."
    echo -e "\t -e\t Map variant predictions on reference genome with their unitig or contig extensions."
    echo -e "\t\t Useless unless mapping on reference genome is required (option -G). "

    echo -e "\nCONTEXT OPTION"
    echo -e "\t -z\t Considers AAV context: all predictions are coherents and authorize non closed bubbles (-x option)" 
    
    echo -e "\nMISC."
    echo -e "\t -w\t Wraith mode: only show all discoSnp++ commands without running them"
    echo -e "\t -v <0 or 1>"
    echo -e "\t\t Verbose 0 (avoids progress output) or 1 (enables progress output) -- default=1."
    echo -e "\t -h | --help"
    echo -e "\t\t Prints this message and exist\n"
    
    echo "Any further question: read the readme file or contact us via the Biostar forum: https://www.biostars.org/t/discosnp/"
}




echo "${yellow}"

while :; do
    case $1 in
    -A) 
        option_phase_variants="-phasing"
        echo "Will phase variants during kissreads process - WARNING this option is too experimental and thus not described in the help message"
        echo "You can obtain clusters using script : \"script/from_phased_alleles_to_clusters.sh file_name_of_phased_alleles\" (the filename(s) is/are given during kissreads process"
        ;;
    -X)
        stop_after_kissnp=1
        ;;

    -z)
        aav=1
        ;;
    -w)
        wraith="true"
        ;;
    -R)
        useref="true"
        output_coverage_option="-dont_output_first_coverage"
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

    -s|--symmetrical)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            option_max_symmetrical_crossroads="-max_symmetrical_crossroads "$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;



    -T|--contigs)
        extend="-T"
        ;;

    -g|--graph)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            graph_reused=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;


    -n|--no_genotype)
        genotyping=""
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

    --fof_mapping)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            read_sets_kissreads=$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;
        
    -b|--branching)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            b=$2
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

    B|--bwa_path)
        if [ "$2" ] && [ ${2:0:1} != "-" ] ; then
            bwa_path_option="-B "$2
            shift
        else
            die 'ERROR: "'$1'" option requires a non-empty option argument.'
        fi
        ;;

    -x)
        x="-x" ##Authorise non closed bubbles
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
        e="-e" # map with extensions to the reference genome
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

if [ -z "$read_sets" ]; then
    echo "$red You must provide at least one read set (-r|--fof)"
    help
    echo $reset
    exit 1
fi


######### CHECK THE k PARITY ##########
rest=$(( $k % 2 ))
if [ $rest -eq 0 ]
then
    echo "$red k=$k is even number, to avoid palindromes, we set it to $(($k-1)) $reset"
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
mapped_readsFilesDump=${prefix}_mapped_read_files_correspondance.txt


#######################################
c_dbgh5=$c
rm -f ${read_sets}_${kissprefix}_removemeplease
if [[ "$useref" == "true" ]]; then

    if [ -z "$genome" ]; then
        echo "$red You can't use option -R without providing a reference genome (-G) $reset"
        help
        exit 1
    fi
       
    myrealpath $genome > ${read_sets}_${kissprefix}_removemeplease
    c_dbgh5="1,"$c

fi
cat $read_sets >> ${read_sets}_${kissprefix}_removemeplease

if [ $aav -eq 1 ]; then
    x="-x"
fi


#######################################################################
#################### OPTIONS SUMMARY            #######################
#######################################################################
if [[ "$wraith" == "false" ]]; then
    echo -e "$yellow Running discoSnp++ "$version", in directory "$EDIR" with following parameters:"
    echo -e "\t read_sets="$read_sets
    echo -e "\t prefix="$h5prefix
    echo -e "\t c="$c
    echo -e "\t C="$C
    echo -e "\t k="$k
    echo -e "\t b="$b
    echo -e "\t d="$d
    echo -e "\t D="$D
    echo -e "\t s="$option_max_symmetrical_crossroads
    echo -e "\t P="$P
    if [ ! -z "${read_sets_kissreads}" ]; then
        echo -e "\t fof_mapping read_file_of_files="${read_sets_kissreads}
    fi
    if [ -f ${graph_reused} ]; then
        echo -e "\t reuse graph="${graph_reused}
    fi
    echo -e "\t p="$prefix
    echo -e "\t G="$genome
    echo -e "\t e="$e
    echo -e "\t x="$x
    
    
    
    echo -e -n "\t starting date="
    date
    echo
fi
echo $reset
#######################################################################
#################### END OPTIONS SUMMARY        #######################
#######################################################################

#############################################################
#################### DUMP READ FILES  #######################
#############################################################
${read_file_names_bin} -in $read_sets > $readsFilesDump
if [ ! -z "${read_sets_kissreads}" ]; then
    ${read_file_names_bin} -in ${read_sets_kissreads} > $mapped_readsFilesDump
fi


############################################################
#################### GRAPH CREATION  #######################
############################################################
if [ ! -f ${graph_reused} ]; then # no graph was given or the given graph was not a file. 
    T="$(date +%s)"
    echo -e "$yellow ############################################################"
    echo -e " #################### GRAPH CREATION  #######################"
    echo -e " ############################################################$reset"

    graphCmd="${dbgh5_bin} -in ${read_sets}_${kissprefix}_removemeplease -out $h5prefix -kmer-size $k -abundance-min ${c_dbgh5} -abundance-max $C -solidity-kind one ${option_cores_gatb} -verbose $verbose  -skip-bcalm -skip-bglue -no-mphf -histo-max 1000000"
    echo $green${graphCmd}$cyan
    if [[ "$wraith" == "false" ]]; then
        ${graphCmd}
    fi
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with graph construction$ reset"
        exit 1
    fi

    T="$(($(date +%s)-T))"
    if [[ "$wraith" == "false" ]]; then
        echo "$yellow Graph creation time in seconds: ${T}$reset"
    fi
    graph_reused=$h5prefix.h5
else
    if [[ "$wraith" == "false" ]]; then
        echo -e "$yellow File ${graph_reused} exists. We use it as input graph$reset"
    fi
fi

echo $reset
cleanCmd="rm -rf trashme_*"
echo $green${cleanCmd}$cyan
if [[ "$wraith" == "false" ]]; then
    ${cleanCmd}
fi
echo $reset

######################################################
#################### KISSNP2   #######################
######################################################
T="$(date +%s)"
echo -e "$yellow ############################################################"
echo -e " #################### KISSNP2 MODULE  #######################"
echo -e " ############################################################$reset"
kissnp2Cmd="${kissnp2_bin} -in ${graph_reused} -out $kissprefix  -b $b $l $x -P $P  -D $D $extend $option_cores_gatb -coverage_file ${h5prefix}_cov.h5 -max_ambigous_indel ${max_ambigous_indel} ${option_max_symmetrical_crossroads}  -verbose $verbose"
echo $green${kissnp2Cmd}$cyan
if [[ "$wraith" == "false" ]]; then
    ${kissnp2Cmd}
fi

if [ $? -ne 0 ]
then
    echo "$red there was a problem with kissnp2$reset"
    exit 1
fi

T="$(($(date +%s)-T))"
if [[ "$wraith" == "false" ]]; then
    echo "$yellow Bubble detection time in seconds: ${T}$reset"
fi

if [ ! -f $kissprefix.fa ]
then
    if [[ "$wraith" == "false" ]]; then
        echo "$yellow No polymorphism predicted by discoSnp++"
        echo -e -n "\t ending date="
        date
        echo -e " Thanks for using discoSnp++ - http://colibread.inria.fr/discoSnp/$reset"
        exit 
    fi
fi

if [ $stop_after_kissnp -eq 1 ]; then
    if [[ "$wraith" == "false" ]]; then
        echo "$yellow -X option detected, computation stopped after variant detection."
        echo "Results (with no read coverage) are located here: "$kissprefix.fa
        echo -e -n "\t ending date="
        date
        echo -e " Thanks for using discoSnp++ - http://colibread.inria.fr/discoSnp/ $reset"
        exit 
    fi  
fi

#######################################################################
#################### KISSREADS                  #######################
#######################################################################

T="$(date +%s)"
    echo -e "$yellow #############################################################"
    echo -e " #################### KISSREADS MODULE #######################"
    echo -e " #############################################################$reset"


smallk=$k
if (( $smallk>31 ))  ; then
    smallk=31
fi
i=5 #avoid modidy this (or increase this if memory needed by kissread is too high. Min 1. Large i (7-10) decreases memory and increases time).
index_stride=$(($i+1)); size_seed=$(($smallk-$i)) # DON'T modify this.

if [ ! -z "${read_sets_kissreads}" ]; then
    read_sets=${read_sets_kissreads}
fi

if [ $aav -eq 1 ]; then
    rmcmd="rm -f ${h5prefix}_cov.h5"
    echo $green${rmcmd}$cyan
    $rmcmd # erase this file so all predictions are coherent after kissreads
fi
kissreadsCmd="${kissreads2_bin} -predictions $kissprefix.fa -reads  ${read_sets}_${kissprefix}_removemeplease -co ${kissprefix}_coherent -unco ${kissprefix}_uncoherent -k $k -size_seeds ${size_seed} -index_stride ${index_stride} -hamming $d  $genotyping -coverage_file ${h5prefix}_cov.h5 $option_cores_gatb  -verbose $verbose ${option_phase_variants}"

echo $green $kissreadsCmd$cyan
if [[ "$wraith" == "false" ]]; then
$kissreadsCmd
fi

if [ $? -ne 0 ]
then
    echo "$red there was a problem with kissreads2$reset":
    exit 1
fi
echo $reset
T="$(($(date +%s)-T))"
# echo "Kissreads (mapping reads on bubbles) time in seconds: ${T}"


#######################################################################
#################### SORT AND FORMAT  RESULTS #########################
#######################################################################

if [[ "$wraith" == "false" ]]; then
echo -e "$yellow ###############################################################"
echo -e " #################### SORT AND FORMAT  RESULTS #################"
echo -e " ###############################################################$reset"
fi

if [[ "$wraith" == "false" ]]; then
    sort -rg ${kissprefix}_coherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_coherent.fa
fi
if [ $? -ne 0 ]
then
    echo "$red there was a problem with the result sorting.$reset"
    exit 1
fi

if [[ "$wraith" == "false" ]]; then
    sort -rg ${kissprefix}_uncoherent | cut -d " " -f 2 | tr ';' '\n' > ${kissprefix}_uncoherent.fa
fi

if [ $? -ne 0 ]
then
    echo "$red there was a problem with the result sorting$reset"
    exit 1
fi

rm -f $kissprefix.fa ${kissprefix}_coherent ${kissprefix}_uncoherent
rm -rf ${read_sets}_${kissprefix}_removemeplease 

#######################################################################
#################### DISCOSNP FINISHED ###############################
#######################################################################






#######################################################################
#################### Deal with VCF ###############################
#######################################################################

T="$(date +%s)"
echo -e "$yellow ###############################################################"
echo -e " #################### CREATE VCF         #######################"
echo -e " ############################################################### $reset"

if [ -z "$genome" ]; then #  NO reference genome use, vcf creator mode 1
    vcfCreatorCmd="$EDIR/scripts/run_VCF_creator.sh -p ${kissprefix}_coherent.fa -o ${kissprefix}_coherent.vcf"
    echo $green$vcfCreatorCmd$cyan
    if [[ "$wraith" == "false" ]]; then
        $vcfCreatorCmd
    fi
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with VCF creation. See how to use the \"run_VCF_creator.sh\" alone.$reset"
        exit 1
    fi
else # A Reference genome is provided, vcf creator mode 2
    vcfCreatorCmd="$EDIR/scripts/run_VCF_creator.sh $bwa_path_option -G $genome $bwa_path_option -p ${kissprefix}_coherent.fa -o ${kissprefix}_coherent.vcf  -I $option_cores_post_analysis $e"
    echo $green$vcfCreatorCmd$cyan
    if [[ "$wraith" == "false" ]]; then
        $vcfCreatorCmd
    fi 
    
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with VCF creation. See how to use the \"run_VCF_creator.sh\" alone.$reset"
        exit 1
    fi
fi
echo $reset
T="$(($(date +%s)-T))"
if [[ "$wraith" == "false" ]]; then
    echo "$yellow Vcf creation time in seconds: ${T}"
    
    echo -e " ###############################################################"
    echo -e " #################### DISCOSNP++ FINISHED ######################"
    echo -e " ###############################################################"
    Ttot="$(($(date +%s)-Ttot))"
    echo "DiscoSnp++ total time in seconds: ${Ttot}"
    echo -e "################################################################################################################"
    echo -e " fasta of predicted variant is \""${kissprefix}_coherent.fa"\""
    
    if [ -z "$genome" ]; then
        echo -e " Ghost VCF file (1-based) is \""${kissprefix}_coherent.vcf"\""
    else
        echo -e " VCF file (1-based) is \""${kissprefix}_coherent.vcf"\""
        echo -e " An IGV ready VCF file (sorted by position, only mapped variants, 0-based) is \""${kissprefix}_coherent_for_IGV.vcf"\""
    fi
    echo -e " Thanks for using discoSnp++ - http://colibread.inria.fr/discoSnp/ - Forum: http://www.biostars.org/t/discoSnp/"
    echo -e "################################################################################################################$reset"
fi