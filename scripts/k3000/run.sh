phased_allele_file=$1
disco_fa_file=$2
EDIR=$( python -c "import os.path; print(os.path.dirname(os.path.realpath(\"${BASH_SOURCE[0]}\")))" ) # as suggested by Philippe Bordron 
if [ -z "$phased_allele_file" ]; then
    echo "You must provide a phased file name"
    exit 1
fi

if [ -z "$disco_fa_file" ]; then
    echo "You must provide a disco file (\"...coherent.fa\")"
    exit 1
fi


# Creating a file where simple paths are compacted
echo "           ### Creating a file where simple paths are compacted"
python3 ${EDIR}/K3000.py ${phased_allele_file} >  compacted_facts_int.txt

# Creating a file with sequences of the compacted paths and removing uncoherent compactions
echo "           ### Creating a fasta file from compacted facts"
python3 ${EDIR}/K3000_compacted_paths_to_fa.py ${disco_fa_file} compacted_facts_int.txt > compacted_facts.fa

/usr/bin/grep ">" compacted_facts.fa | cut -d ">" -f 2 > compacted_facts_int.txt


### Creating a GFA graph
echo "           ### Creating a GFA graph"
python3 ${EDIR}/K3000_msr_to_gfa.py compacted_facts_int.txt > compacted_facts.gfa 

# Adding paired edges and counting of compacted facts
echo "           ### Adding paired edges and counting of compacted facts"
python3 ${EDIR}/enhance_gfa.py compacted_facts.gfa ${phased_allele_file} > graph.gfa

# Detecting snp succesion
echo "           ### Detecting snp succession"
python3 ${EDIR}/find_unitig_connected_pairs_of_facts.py graph.gfa ${disco_fa_file} > graph_plus.gfa

# cleanup
#rm -f compacted_facts_int.txt compacted_facts.gfa graph.gfa
