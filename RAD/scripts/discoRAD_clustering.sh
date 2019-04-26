#!/bin/sh

# REQUIRES:
## Python 3
## short_read_connector: installed and compiled: https://github.com/GATB/short_read_connector
# echo "WARNING: short_read_connector must have been compiled"


function help {
echo "====================================================="
echo "Filtering, clustering per locus and vcf formatting of $rawdiscofile"
echo "====================================================="
echo "this script manages bubble clustering from a discofile.fa file, and the integration of cluster informations in a vcf file"
echo " 1/ Remove variants with more than 95% missing genotypes and low rank (<0.4)"
echo " 2/ Cluster variants per locus"
echo " 3/ Format the variants in a vcf file with cluster information"
echo "Usage: ./discoRAD_clustering.sh -f discofile -s SRC_directory/ -o output_file.vcf"
echo "nb: all options are MANDATORY\n"
echo "OPTIONS:"
echo "\t -f: DiscoSnp fasta output containing coherent predictions"
echo "\t -s: Path to Short Read Connector"
echo "\t -o: output file path (vcf)"
}


function sumup {
echo "====================================================="
echo "Filtering of $rawdiscofile"
echo "====================================================="
echo " 1/ Remove variants with more than ${percent_missing} missing genotypes and low rank (<${min_rank})"
echo " 2/ Clustering variants (sharing at least a ${usedk}-mers)"
#echo " 3/ Removing variant in clusters with a size above ${max_cluster_size}"
#echo " 4/ Removing low ranked variants (those whose rank is < ${min_rank}"
echo " Resulting file is ${rawdiscofile_base}_clustered.vcf"
}


if [ "$#" -ne 6 ]; then
help
exit
fi

while getopts "f:s:o:h" opt; do
    case $opt in
        f)
        rawdiscofile=$OPTARG
        ;;

        s)
        short_read_connector_directory=$OPTARG
        ;;

        o)
        output_file=$OPTARG
        ;;

        h)
        help
        exit
        ;;
    esac
done

# Detect the directory path

EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
BINDIR=$EDIR"/../../build/bin"
rawdiscofile_base=$( basename  "${rawdiscofile}" .fa)

#################### PARAMETERS VALUES #######################

#Get k value (for clustering purpose)
originalk=$( echo $rawdiscofile | awk -F k_ '{ print $2 }' | cut -d "_" -f 1)
usedk=$((originalk-1))

# rank filter parameter
min_rank=0.4

# cluster size parameter  -- no longer used
# max_cluster_size=3000000

percent_missing=0.95

echo "############################################################"
echo "######### MISSING DATA AND LOW RANK FILTERING  #############"
echo "############################################################"

echo "Filtering variants with more than 0.95 missing data and rank<0.4 ..."

disco_filtered=${rawdiscofile_base}_filtered
python3 ${EDIR}/../../scripts/create_filtered_vcf.py -i ${rawdiscofile} -f -o ${disco_filtered}.fa -m ${percent_missing} -r ${min_rank}


######################### Clustering ###########################

echo "############################################################"
echo "###################### CLUSTERING ##########################"
echo "############################################################"

# Simplify headers (for dsk purposes)
disco_simpler=${disco_filtered}_simpler
cat ${disco_filtered}.fa | cut -d "|" -f 1 | sed -e "s/^ *//g" > ${disco_simpler}.fa

ls ${disco_simpler}.fa > ${disco_simpler}.fof

# Compute sequence similarities
cmdSRC="${short_read_connector_directory}/short_read_connector.sh -b ${disco_simpler}.fa -q ${disco_simpler}.fof -s 0 -k ${usedk} -a 1 -l -p ${disco_simpler}"
echo $cmdSRC
eval $cmdSRC

if [ $? -ne 0 ]
then
    echo "there was a problem with Short Read Connector, exit"
    exit 1
fi

# Format one line per edge
cmd="python3 ${EDIR}/from_SRC_to_edges.py ${disco_simpler}.txt"
echo $cmd "> ${disco_simpler}_edges.txt"
eval $cmd "> ${disco_simpler}_edges.txt"


# Compute the clustering
cmdqhc="${BINDIR}/quick_hierarchical_clustering ${disco_simpler}_edges.txt"
echo $cmdqhc " > ${disco_simpler}.cluster"
eval $cmdqhc "> ${disco_simpler}.cluster"

if [ $? -ne 0 ]
then
    echo "there was a problem with quick_hierarchical_clustering, exit"
    exit 1
fi
# Generate a .fa file with clustering information
disco_final="${disco_filtered}_with_clusters"
cmd="python3 ${EDIR}/clusters_and_fasta_to_fasta.py ${disco_filtered}.fa ${disco_simpler}.cluster"
echo $cmd " > ${disco_final}.fa"
eval $cmd "> ${disco_final}.fa"
if [ $? -ne 0 ]
then
    echo "there was a problem when generating a .fa file with clustering information, exit"
    exit 1
fi

######################### VCF generation with cluster information ###########################

echo "############################################################"
echo "###################### OUTPUT VCF ##########################"
echo "############################################################"

cmdVCF="python3 ${EDIR}/../../scripts/create_filtered_vcf.py -i ${disco_final}.fa -o ${output_file}"
echo $cmdVCF
eval $cmdVCF

if [ $? -ne 0 ]
then
    echo "there was a problem with vcf creation, exit"
    exit 1
fi

######################### Filter large clusters NO LONGER USED ###########################

#echo "############################################################"
#echo "################### FILTER CLUSTERS ########################"
#echo "############################################################"
#
#cmdclustsize="${EDIR}/filter_by_cluster_size.sh ${original_disco}_with_clusters.vcf ${max_cluster_size}"
#echo $cmdclustsize
#eval $cmdclustsize
#
#if [ $? -ne 0 ]
#then
#    echo "there was a problem with cluster size fitlering"
#    exit 1
#fi


echo "#######################################################################"
echo "#########################  CLEANING  ####################################"
echo "#######################################################################"

#TO UNCOMMENT
#rm -f ${disco_simpler}*
#rm -f ${disco_filtered}.fa
#rm -f ${disco_final}.fa


sumup > log_${rawdiscofile_base}_sorted_with_clusters.txt

echo "============================"
echo " DISCORAD FINALIZATION DONE "
echo "============================"
echo " Results in ${rawdiscofile_base}_sorted_with_clusters.vcf"
echo " Logs in log_${rawdiscofile_base}_sorted_with_clusters.txt"

