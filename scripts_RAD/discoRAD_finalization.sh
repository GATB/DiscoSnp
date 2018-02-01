#!/bin/sh

# REQUIRES:
## Python 3
## short_read_connector: installed and compiled: https://github.com/GATB/short_read_connector
## quick_hierarchical_clustering compiled:
### c++ -std=gnu++11 quick_hierarchical_clustering.cpp -o quick_hierarchical_clustering     # WITH LINUX
### clang++ -std=gnu++11 quick_hierarchical_clustering.cpp -o quick_hierarchical_clustering # WITH MACOS

# echo "WARNING1: short_read_connector must have been compiled"
# echo "WARNING2: quick_hierarchical_clustering.cpp must have been compiled :"
# echo "  c++ -std=gnu++11 quick_hierarchical_clustering.cpp -o quick_hierarchical_clustering     # WITH LINUX"
# echo "  clang++ -std=gnu++11 quick_hierarchical_clustering.cpp -o quick_hierarchical_clustering # WITH MACOS"


function help {
echo "====================================================="
echo "Filtration of $rawdiscofile"
echo "====================================================="
echo "run discoRAD.sh, this script manages bubble clustering from a discofile.fa file, and the integration of cluster informations in a disco.vcf file"
echo " 1/ Remove variants with more than 95% genotypes"
echo " 2/ Clustering variants"
echo " 3/ Removing heterozygous clusters (paralogs filter)"
echo " 4/ Removing low ranked variants (those whose rank is < 0.2 \n"
echo "Usage: ./discoRAD.sh -f discofile -s SRC_directory/"
echo "nb: all options are MANDATORY\n"
echo "OPTIONS:"
echo "\t -f: DiscoSnp output containing coherent predictions"
echo "\t -s: Path to Short Read Connector"
}


function sumup {
echo "====================================================="
echo "Filtration of $rawdiscofile"
echo "====================================================="
echo " 1/ Remove variants with more than ${percent_missing} genotypes"
echo " 2/ Clustering variants (sharing at least a ${usedk}-mers)"
echo " 3/ Removing clusters with more than 100*${max_hetero}% heterozygous variants in more than 100*${max_indivs}% individuals"
echo " 4/ Removing low ranked variants (those whose rank is < ${min_rank}"
echo " Resulting file is ${rawdiscofile_base}_sorted_with_clusters.vcf"
}


if [ "$#" -ne 4 ]; then
help
exit
fi

while getopts "f:s:h" opt; do
    case $opt in
        f)
        rawdiscofile=$OPTARG
        ;;

        s)
        short_read_connector_directory=$OPTARG
        ;;

        h)
        help
        exit
        ;;
    esac
done

# Detect the directory path

EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

rawdiscofile_base=$( basename  "${rawdiscofile}" .fa)

#################### PARAMETERS VALUES #######################

#Get k value (for clustering purpose)
originalk=$( echo $rawdiscofile | awk -F k_ '{ print $2 }' | cut -d "_" -f 1)
usedk=$((originalk-1))

# rank filter parameter
min_rank=0.2

# paralogous filter parameters
max_hetero=0.1
max_indivs=0.5

percent_missing=0.95

echo "############################################################"
echo "############### MISSING DATA FILTERING  ####################"
echo "############################################################"

echo "Filtering variants with more than 0.95 missing data ..."

python3 ${EDIR}/filter_missgeno.py ${rawdiscofile} ${percent_missing}

mv ${percent_missing}missing_${rawdiscofile_base}.fa ERASEME_${percent_missing}missing_${rawdiscofile_base}.fa
original_disco=ERASEME_${percent_missing}missing_${rawdiscofile_base}

before=$( grep 'higher' $rawdiscofile | wc -l )
after=$( grep 'higher' ${original_disco}.fa | wc -l )

echo "Done, "$after"/"$before" variants conserved for clustering"


######################### Clustering ###########################

echo "############################################################"
echo "###################### CLUSTERING ##########################"
echo "############################################################"

# Simplify headers (for dsk purposes)
cat ${original_disco}.fa | cut -d "|" -f 1 | sed -e "s/^ *//g" > ${original_disco}_simpler.fa
discofile=${original_disco}_simpler
ls ${discofile}.fa > ${discofile}.fof

# Compute sequence similarities
cmdSRC="${short_read_connector_directory}/short_read_connector.sh -b ${discofile}.fa -q ${discofile}.fof -s 0 -k ${usedk} -a 1 -l -p ${discofile}"
$cmdSRC

if [ $? -ne 0 ]
then
    echo "there was a problem with Short Read Connector, exit"
    exit 1
fi

# Compute the clustering
${EDIR}/quick_hierarchical_clustering ${discofile}.txt > ${discofile}.cluster

# Generate a .fa file with clustering information
python3 ${EDIR}/clusters_and_fasta_to_fasta.py ${original_disco}.fa ${discofile}.cluster > ${original_disco}_with_clusters.fa


######################### VCF generation with cluster information ###########################

echo "############################################################"
echo "###################### OUTPUT VCF ##########################"
echo "############################################################"

cmdVCF="${EDIR}/../scripts/run_VCF_creator.sh -p  ${original_disco}_with_clusters.fa -o ${original_disco}_with_clusters.vcf"
$cmdVCF

if [ $? -ne 0 ]
then
    echo "there was a problem with vcf creation, exit"
    exit 1
fi

######################### Filter suspicious paralogous clusters ###########################

echo "###########################################################"
echo "################### FILTER PARALOGS ########################"
echo "############################################################"

cmdpara="python3 ${EDIR}/filter_paralogs.py ${original_disco}_with_clusters.vcf ${max_hetero} ${max_indivs}"
$cmdpara

if [ $? -ne 0 ]
then
    echo "there was a problem with paralogs fitlering"
    exit 1
fi

# Remove low ranked variants

cmdrk="python3 ${EDIR}/filter_rank_vcf.py para_${max_hetero}_${max_indivs}_${original_disco}_with_clusters.vcf ${min_rank}"
$cmdrk

if [ $? -ne 0 ]
then
    echo "there was a problem with rank filtering"
    exit 1
fi

# Sort the .vcf file
grep ^# ${min_rank}rk_para_${max_hetero}_${max_indivs}_${original_disco}_with_clusters.vcf > ${original_disco}_with_sorted_clusters.vcf; grep -v ^# ${original_disco}_with_clusters.vcf | sort >> ${original_disco}_with_sorted_clusters.vcf;

# Format chromosome Names in the VCF
python3 ${EDIR}/format_VCF_with_cluster_ids.py ${original_disco}_with_sorted_clusters.vcf > ${original_disco}_with_sorted_formatted_clusters.vcf

# Clean results
mv ${original_disco}_with_sorted_formatted_clusters.vcf ${rawdiscofile_base}_sorted_with_clusters.vcf

rm -f *ERASEME*
sumup > log_${rawdiscofile_base}_sorted_with_clusters.txt

echo "============================"
echo " DISCORAD FINALIZATION DONE"
echo "============================"
echo " Results in ${rawdiscofile_base}_sorted_with_clusters.vcf"
echo " Logs in log_${rawdiscofile_base}_sorted_with_clusters.txt"

