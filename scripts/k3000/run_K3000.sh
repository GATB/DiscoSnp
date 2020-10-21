red=`tput setaf 1`
green=`tput setaf 2`
cyan=`tput setaf 6`
bold=`tput bold`
reset=`tput sgr0`
underline=`tput smul`
no_underline=`tput rmul`


wraith=false #$4 # set to true if you only want to see commands without executing them
phased_allele_file=$1
disco_cofa_file=$2
disco_uncofa_file=$3
read_set_id=$4
EDIR=$( python -c "import os.path; print(os.path.dirname(os.path.realpath(\"${BASH_SOURCE[0]}\")))" ) # as suggested by Philippe Bordron 
if [ -z "$phased_allele_file" ]; then
    echo "${red}           You must provide a phased file name${reset}"
    echo "${red}           Usage: $0 phased_allele_file disco_coherent.fa disco_uncoherent.fa read_set_id${reset}"
    exit 1
fi

if [ -z "$disco_cofa_file" ]; then
    echo "${red}           You must provide a disco (coherent) file (\"...coherent.fa\")${reset}"
    echo "${red}           Usage: $0 phased_allele_file disco_coherent.fa disco_uncoherent.fa read_set_id${reset}"
    exit 1
fi

if [ -z "$disco_uncofa_file" ]; then
    echo "${red}           You must provide a disco (uncoherent) file (\"...uncoherent.fa\")${reset}"
    echo "${red}           Usage: $0 phased_allele_file disco_coherent.fa disco_uncoherent.fa read_set_id${reset}"
    exit 1
fi

if [ -z "$read_set_id" ]; then
    echo "${red}           You must provide a read set id (integer value from 1 to the number of read set used to create file $disco_cofa_file ${reset}"
    echo "${red}           Usage: $0 phased_allele_file disco_coherent.fa disco_uncoherent.fa read_set_id${reset}"
    exit 1
fi



# Creating a file where simple paths are compacted

echo "${green}${bold}           ### EXPLOITATION OF PHASING INFORMATION OBTAINED FROM DISCOSNP${reset}"
echo "${green}           ### Input phased_allele_file: `sha1sum ${phased_allele_file}`${reset}"
echo "${green}           ### Input disco_cofa_file:      `sha1sum ${disco_cofa_file}`${reset}"
echo "${green}           ### Input disco_uncofa_file:      `sha1sum ${disco_uncofa_file}`${reset}"
echo "${green}           ### Input read_set_id:        ${read_set_id}${reset}"
echo ""

# Prefiltering: removing non perfectly overlapping facts
echo "${green}           ### Removing non perfectly overlapping facts"
cmd="python3 ${EDIR}/K3000_filter_badly_overlapping_variants.py ${disco_cofa_file} ${disco_uncofa_file} ${phased_allele_file}"
echo "           "$cmd "> filtered_${phased_allele_file}${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > filtered_${phased_allele_file}
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi



# Determining working zones: 
echo "${green}           ### Determining working zones"
cmd="python3 ${EDIR}/K3000_working_zone_no_redundant_edges.py filtered_${phased_allele_file}"
echo "           "$cmd " > wz_filtered_${phased_allele_file}${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > wz_filtered_${phased_allele_file}
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


# Creating a file where simple paths are compacted
echo "${green}           ### Creating a file where simple paths are compacted"
cmd="python3 ${EDIR}/K3000.py wz_filtered_${phased_allele_file}"
echo "           "$cmd "> compacted_facts_int_${read_set_id}.txt${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > compacted_facts_int_${read_set_id}.txt
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi




# Creating a file with sequences of the compacted paths and removing uncoherent compactions
echo "${green}           ### Creating a fasta file from compacted facts"
cmd="python3 ${EDIR}/K3000_paths_to_fa.py ${disco_cofa_file} ${disco_uncofa_file} compacted_facts_int_${read_set_id}.txt" 
echo "           "$cmd "> compacted_facts_${read_set_id}.fa${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > compacted_facts_${read_set_id}.fa
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi

echo "${green}           ### Select only valid facts and add their positions"
# /usr/bin/grep ">" compacted_facts.fa | cut -d ">" -f 2 | cut -d " " -f 1 > compacted_facts_int.txt
cmd="/usr/bin/grep \">\" compacted_facts_${read_set_id}.fa | cut -d \">\" -f 2"
echo "           "$cmd "> compacted_facts_int_${read_set_id}.txt${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > compacted_facts_int_${read_set_id}.txt
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


### Creating a GFA graph
echo "${green}           ### Creating a GFA graph"
# python3 ${EDIR}/K3000_facts_to_gfa.py compacted_facts_int.txt > compacted_facts.gfa 
cmd="python3 ${EDIR}/K3000_facts_to_gfa.py compacted_facts_int_${read_set_id}.txt"
echo "           "$cmd "> compacted_facts_${read_set_id}.gfa${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > compacted_facts_${read_set_id}.gfa
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi



# Adding paired edges and counting of compacted facts
echo "${green}           ### Adding paired edges and counting of compacted facts"
cmd="python3 ${EDIR}/K3000_enhance_gfa.py compacted_facts_${read_set_id}.gfa wz_filtered_${phased_allele_file} ${disco_cofa_file} ${disco_uncofa_file} ${read_set_id}"
echo "           "$cmd" > graph_${read_set_id}.gfa${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > graph_${read_set_id}.gfa
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


# Detecting snp succesion
echo "${green}           ### Detecting snp succession"
# python3 ${EDIR}/find_unitig_connected_pairs_of_facts.py graph.gfa ${disco_cofa_file} > graph_plus.gfa
cmd="python3 ${EDIR}/K3000_find_unitig_connected_pairs_of_facts.py graph_${read_set_id}.gfa ${disco_cofa_file} ${disco_uncofa_file}"
echo "           "$cmd" > graph_plus_${read_set_id}.gfa${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > graph_plus_${read_set_id}.gfa
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


echo "${green}           ### Create final graph with sequence content"
# python3 ${EDIR}/K3000_node_ids_to_node_sequences.py graph_plus.gfa compacted_facts.fa > graph_final.gfa
cmd="python3 ${EDIR}/K3000_node_ids_to_node_sequences.py graph_plus_${read_set_id}.gfa compacted_facts_${read_set_id}.fa"
echo "           "$cmd" > graph_final_${read_set_id}.gfa${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd > graph_final_${read_set_id}.gfa
fi
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi

echo "${green}           ### Create stats (requires mathplotlib)"
cmd="python3 ${EDIR}/stats.py graph_plus_${read_set_id}.gfa graph_final_${read_set_id}.gfa ${read_set_id}"
echo "           "$cmd"${cyan}"
if [[ "$wraith" == "false" ]]; then
    eval $cmd
fi
if [ $? -ne 0 ]
then
    echo "${red}           ### Non critical problem detected, check logs.${reset}"
    echo "${green}           ### You may remove useless files: rm -f compacted_facts_int_${read_set_id}.txt compacted_facts_${read_set_id}.gfa graph_${read_set_id}.gfa compacted_facts_${read_set_id}.fa "
    echo "${green}${bold}           ### EXPLOITATION OF PHASING INFORMATION OBTAINED FROM DISCOSNP ENDED, the final graph is $(tput blink)${underline}graph_final_${read_set_id}.gfa${no_underline}${reset} "
    exit 0
fi

# cleanup
echo echo "${green}           ### You may remove useless files: rm -f compacted_facts_int_${read_set_id}.txt compacted_facts_${read_set_id}.gfa graph_${read_set_id}.gfa compacted_facts_${read_set_id}.fa "
#rm -f compacted_facts_int.txt compacted_facts.gfa graph.gfa compacted_facts.fa graph_plus.gfa 

echo "${green}${bold}           ### EXPLOITATION OF PHASING INFORMATION OBTAINED FROM DISCOSNP ENDED, the final graph is ${underline}graph_final_${read_set_id}.gfa${no_underline}${green}${reset} ${green}${bold} stats are available in ${underline}distributions_${read_set_id}.png${no_underline}${reset}"
