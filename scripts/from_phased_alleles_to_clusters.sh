# FROM A FILE GENERATED USING OPTION -phasing FROM KISSREADS, GENERATE CLUSTERS OF VARIANTS WHO TRANSITIVELY CO-OCCUR IN READS

echo "This is \"from_phased_alleles_to_clusters.sh\". Generates clusters of variants who transitively co-occur in reads from a file generated using option -phasing from kissreads"

file=$1
if [ ! -f $file ]; then
    echo "In \"from_phased_alleles_to_clusters.sh\": file \"$file\" does not exist"
   exit 0
fi
filename=$(basename -- "$file")
path=$(dirname "${file}")

# FIND THE PATH CONTAINING THE SCRIPT: 
EDIR=$( python -c "import os.path; print(os.path.dirname(os.path.realpath(\"${BASH_SOURCE[0]}\")))" ) 
echo $EDIR


# 1000547h;-2286435h; -1330792h;1152525l; => 1
# FORMAT PHASED ALLELE IDS INTO SIMPLER FORMAT FOR CONNECTED COMPONENT DETECTION
cmd="cat ${file} | tr -d \"-\" | sed '1d' | cut -d \"=\" -f 1 | python3 ${EDIR}/from_path_to_edges.py | cut -f 1,2 | tr -d \"l\" | tr -d \"h\" | sort -u "
# 1/ suppress the '-' sign occurrences 
# 2/ suppres the first line 
# 3/ remove what exists after => (included) 
# 4/ transform n-uplets into couples: 
#   a;b;c; d;e
#  becomes (r means mapped by a unique gene and p means mapped by a pair of reads). Each line is ordered (smallest id first)
#   a b r
#   b c r
#   c d p
#   d e r
# 5/ remove r/p info
# 6/ remove h/l info
# 7/ remove duplicates
echo $cmd "> edge_${filename}"
eval $cmd "> edge_${filename}"

if [ $? -ne 0 ]
then
    echo "In \"from_phased_alleles_to_clusters.sh\": there was a problem with file $file"
    exit 1
fi

# COMPUTE CONNECTED COMPONENTS
cmd="${EDIR}/../build/bin/quick_hierarchical_clustering edge_${filename}"
echo $cmd "> connected_components_${filename}"
eval $cmd > connected_components_${filename}
if [ $? -ne 0 ]
then
    echo "In \"from_phased_alleles_to_clusters.sh\": there was a problem while computing connected components of edge_${filename}"
    exit 1
fi

# REMOVE INTERMEDIATE FILE
# cmd="rm -f edge_${filename}"
# echo $cmd
# eval $cmd

# REPLACE THE CREATED FILE IN TIS ORIGINAL DIRECTORY
cmd="mv connected_components_${filename} $path/connected_components_${filename}"
echo $cmd
eval $cmd
if [ $? -ne 0 ]
then
    echo "In \"from_phased_alleles_to_clusters.sh\": there was a problem while moving connected components file "
    exit 1
fi

echo "Connected components (clusters of variants) from file $file are in $path/connected_components_${filename}"