red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

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
echo "${red}           ### Creating a file where simple paths are compacted${reset}"
python3 ${EDIR}/K3000.py ${phased_allele_file} >  compacted_facts_int.txt

# Creating a file with sequences of the compacted paths and removing uncoherent compactions
echo "${red}           ### Creating a fasta file from compacted facts${reset}"
python3 ${EDIR}/K3000_compacted_paths_to_fa.py ${disco_fa_file} compacted_facts_int.txt > compacted_facts.fa

/usr/bin/grep ">" compacted_facts.fa | cut -d ">" -f 2 > compacted_facts_int.txt


### Creating a GFA graph
echo "${red}           ### Creating a GFA graph${reset}"
python3 ${EDIR}/K3000_msr_to_gfa.py compacted_facts_int.txt > compacted_facts.gfa 

# Adding paired edges and counting of compacted facts
echo "${red}           ### Adding paired edges and counting of compacted facts${reset}"
python3 ${EDIR}/enhance_gfa.py compacted_facts.gfa ${phased_allele_file} > graph.gfa

# Detecting snp succesion
echo "${red}           ### Detecting snp succession${reset}"
python3 ${EDIR}/find_unitig_connected_pairs_of_facts.py graph.gfa ${disco_fa_file} > graph_plus.gfa


echo "${red}           ### Create final graph with sequence content${reset}"
python3 ${EDIR}/K3000_node_ids_to_node_sequences.py graph_plus.gfa compacted_facts.fa > graph_final.gfa

# cleanup
rm -f compacted_facts_int.txt compacted_facts.gfa graph.gfa compacted_facts.fa graph_plus.gfa 

echo "${green}           ### Phasing ended, the final graph is graph_plus.gfa${reset}"
