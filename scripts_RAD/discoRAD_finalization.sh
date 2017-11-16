# REQUIRES:
## short_read_connector: installed and compiled: https://github.com/GATB/short_read_connector

function sumup {
    echo ""
    echo "  ============================================================================================================"
    echo "  =  Filtration of $rawdiscofile"
    echo "  ============================================================================================================"
    echo "  = 1/ Remove variants with more than ${percent_missing} genotypes"
    echo "  = 2/ Clustering variants (sharing at least a ${usedk}-mers)"
    echo "  = 3/ Removing clusters with more than 100*${max_hetero}% heterozygous variants in more than 100*${max_indivs}% individuals"
    echo "  = 4/ Removing low ranked variants (those whose rank is < ${min_rank})"
    echo "  = Resulting file are "
    echo "  =   - ${rawdiscofile_base}_with_clusters.fa: fasta file after steps 1/ and 2/"
    echo "  =   - raw_${original_disco}_with_sorted_clusters.vcf: vcf file after steps 1/ and 2/"
    echo "  =   - filtered_${original_disco}_with_sorted_clusters.vcf: vcf file after steps 1/, 2/, 3/ and 4/"
    # echo "  =   - log_${rawdiscofile_base}_sorted_with_clusters.txt: contains logs of clustering and filtrations"
    echo "  ============================================================================================================"
}


if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters. Usage: "
    echo "sh discoRAD_finalization.sh discoSnpRAD_file  short_read_connector_path."
    echo "  Example: : sh discoRAD_finalization.sh /tmp/mydata_k_31_c_3_D_10_P_5_b_1_coherent.fa ~/workspace/short_read_connector/"
    echo ""
    echo "short_read_connector sources:  https://github.com/GATB/short_read_connector"
    exit
fi

# Detect the directory path
EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
rawdiscofile=$1
rawdiscofile_base=$( basename  "${rawdiscofile}" .fa)
short_read_connector_directory=$2

#Get k value (for clustering purpose)
originalk=$( echo $rawdiscofile | awk -F k_ '{ print $2 }' | cut -d "_" -f 1)
usedk=$((originalk-1))

# rank filter parameter
min_rank=0.2

# paralogous filter parameters
max_hetero=0.1
max_indivs=0.5

percent_missing=0.95
# Filter missing data
cmd="python3 ${EDIR}/filter_missgeno.py ${rawdiscofile} ${percent_missing}"

echo -e "\t\t$cmd"
$cmd

cmd="mv ${percent_missing}missing_${rawdiscofile_base}.fa ERASEME_${percent_missing}missing_${rawdiscofile_base}.fa"
echo -e "\t\t$cmd"
$cmd

original_disco=ERASEME_${percent_missing}missing_${rawdiscofile_base}

# Simplify headers (for dsk purposes)
cat ${original_disco}.fa | cut -d "|" -f 1 | sed -e "s/^ *//g" > ${original_disco}_simpler.fa
discofile=${original_disco}_simpler

ls ${discofile}.fa > ${discofile}.fof

# Compute sequence similarities
cmd="${short_read_connector_directory}/short_read_connector.sh -b ${discofile}.fa -q ${discofile}.fof -s 0 -k ${usedk} -a 1 -l -p ${discofile}" 
echo -e "\t\t$cmd"
$cmd

# Compute the clustering
cmd="${EDIR}/../build/bin/quick_hierarchical_clustering ${discofile}.txt"
echo -e "\t\t$cmd"
$cmd > ${discofile}.cluster

# Generate a .fa file with clustering information
cmd="python3 ${EDIR}/clusters_and_fasta_to_fasta.py ${original_disco}.fa ${discofile}.cluster"
echo -e "\t\t$cmd"
$cmd > ${original_disco}_with_clusters.fa

# Generate a .vcf file with clustering information
cmd="${EDIR}/../scripts/run_VCF_creator.sh -p ${original_disco}_with_clusters.fa -o ${original_disco}_with_clusters.vcf"
echo -e "\t\t$cmd"
$cmd

# Filter suspicious paralogous clusters
cmd="python3 ${EDIR}/filter_paralogs.py ${original_disco}_with_clusters.vcf ${max_hetero} ${max_indivs}"
echo -e "\t\t$cmd"
$cmd

# Remove low ranked variants
cmd="python3 ${EDIR}/filter_rank_vcf.py para_${max_hetero}_${max_indivs}_${original_disco}_with_clusters.vcf ${min_rank}"
echo -e "\t\t$cmd"
$cmd

# Sort the .vcf files
grep ^# ${min_rank}rk_para_${max_hetero}_${max_indivs}_${original_disco}_with_clusters.vcf > filtered_${original_disco}_with_sorted_clusters.vcf
grep -v ^# ${min_rank}rk_para_${max_hetero}_${max_indivs}_${original_disco}_with_clusters.vcf | sort >> filtered_${original_disco}_with_sorted_clusters.vcf

grep ^# ${original_disco}_with_clusters.vcf > raw_${original_disco}_with_sorted_clusters.vcf
grep -v ^# ${original_disco}_with_clusters.vcf | sort >> raw_${original_disco}_with_sorted_clusters.vcf

# Format chromosome Names in the VCF
for type in 'raw' 'filtered'
    do
        cmd="python3 ${EDIR}/format_VCF_with_cluster_ids.py "$type"_${original_disco}_with_sorted_clusters.vcf "
        echo -e "\t\t$cmd"
        $cmd > "$type"_${original_disco}_with_sorted_formatted_clusters.vcf
    done

# Clean results
for type in 'raw' 'filtered' 
do
    cmd="mv "$type"_${original_disco}_with_sorted_formatted_clusters.vcf "$type"_${rawdiscofile_base}_sorted_with_clusters.vcf"
    echo -e "\t\t$cmd"
    $cmd
done

mv ${original_disco}_with_clusters.fa ${rawdiscofile_base}_with_clusters.fa
rm -f *ERASEME*

sumup

# cmd="sumup > log_${rawdiscofile_base}_sorted_with_clusters.txt"
# echo -e "\t\t$cmd"
# $cmd
#
# echo
# echo "============================"
# echo " DISCORAD FINALIZATION DONE"
# echo "============================"
# echo " Filtered results in filtered_{rawdiscofile_base}_sorted_with_clusters.vcf"
# echo " Raw results in raw_${rawdiscofile_base}_sorted_with_clusters.vcf"
# echo " Logs in log_${rawdiscofile_base}_sorted_with_clusters.txt"
