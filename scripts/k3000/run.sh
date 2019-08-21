red=`tput setaf 1`
green=`tput setaf 2`
cyan=`tput setaf 6`
bold=`tput bold`
reset=`tput sgr0`
underline=`tput smul`
no_underline=`tput rmul`

phased_allele_file=$1
disco_fa_file=$2
EDIR=$( python -c "import os.path; print(os.path.dirname(os.path.realpath(\"${BASH_SOURCE[0]}\")))" ) # as suggested by Philippe Bordron 
if [ -z "$phased_allele_file" ]; then
    echo "${red}           You must provide a phased file name${reset}"
    exit 1
fi

if [ -z "$disco_fa_file" ]; then
    echo "${red}           You must provide a disco file (\"...coherent.fa\")${reset}"
    exit 1
fi


# Creating a file where simple paths are compacted
echo "${green}           ### Creating a file where simple paths are compacted"
cmd="python3 ${EDIR}/K3000.py ${phased_allele_file}"
echo "           "$cmd "> compacted_facts_int.txt${cyan}"
eval $cmd > compacted_facts_int.txt
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


# Creating a file with sequences of the compacted paths and removing uncoherent compactions
echo "${green}           ### Creating a fasta file from compacted facts"
cmd="python3 ${EDIR}/K3000_compacted_paths_to_fa.py ${disco_fa_file} compacted_facts_int.txt" 
echo "           "$cmd "> compacted_facts.fa${cyan}"
eval $cmd > compacted_facts.fa
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi

echo "${green}           ### Select only valid facts"
# /usr/bin/grep ">" compacted_facts.fa | cut -d ">" -f 2 | cut -d " " -f 1 > compacted_facts_int.txt
cmd="/usr/bin/grep \">\" compacted_facts.fa | cut -d \">\" -f 2 | cut -d \" \" -f 1 "
echo "           "$cmd "> compacted_facts_int.txt${cyan}"
eval $cmd > compacted_facts_int.txt
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


### Creating a GFA graph
echo "${green}           ### Creating a GFA graph"
# python3 ${EDIR}/K3000_msr_to_gfa.py compacted_facts_int.txt > compacted_facts.gfa 
cmd="python3 ${EDIR}/K3000_msr_to_gfa.py compacted_facts_int.txt"
echo "           "$cmd "> compacted_facts.gfa${cyan}"
eval $cmd > compacted_facts.gfa
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi



# Adding paired edges and counting of compacted facts
echo "${green}           ### Adding paired edges and counting of compacted facts"
# python3 ${EDIR}/enhance_gfa.py compacted_facts.gfa ${phased_allele_file} > graph.gfa
cmd="python3 ${EDIR}/enhance_gfa.py compacted_facts.gfa ${phased_allele_file}"
echo "           "$cmd" > graph.gfa${cyan}"
eval $cmd > graph.gfa
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


# Detecting snp succesion
echo "${green}           ### Detecting snp succession"
# python3 ${EDIR}/find_unitig_connected_pairs_of_facts.py graph.gfa ${disco_fa_file} > graph_plus.gfa
cmd="python3 ${EDIR}/find_unitig_connected_pairs_of_facts.py graph.gfa ${disco_fa_file}"
echo "           "$cmd" > graph_plus.gfa${cyan}"
eval $cmd > graph_plus.gfa
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi


echo "${green}           ### Create final graph with sequence content"
# python3 ${EDIR}/K3000_node_ids_to_node_sequences.py graph_plus.gfa compacted_facts.fa > graph_final.gfa
cmd="python3 ${EDIR}/K3000_node_ids_to_node_sequences.py graph_plus.gfa compacted_facts.fa"
echo "           "$cmd" > graph_final.gfa${cyan}"
eval $cmd > graph_final.gfa
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    exit 1
fi

echo "${green}           ### Create stats (requires mathplotlib)"
cmd="python3 ${EDIR}/stats.py graph_plus.gfa graph_final.gfa"
echo "           "$cmd"${cyan}"
eval $cmd
if [ $? -ne 0 ]
then
    echo "${red}           ###Problem detected, check logs.${reset}"
    echo "${green}           ### You may remove useless files: rm -f compacted_facts_int.txt compacted_facts.gfa graph.gfa compacted_facts.fa graph_plus.gfa "
    echo "${green}${bold}           ### Phasing ended, the final graph is $(tput blink)${underline}graph_final.gfa${no_underline}.${reset}. "
    exit 0
fi

# cleanup
echo "${green}           ### You may remove useless files: rm -f compacted_facts_int.txt compacted_facts.gfa graph.gfa compacted_facts.fa graph_plus.gfa "
#rm -f compacted_facts_int.txt compacted_facts.gfa graph.gfa compacted_facts.fa graph_plus.gfa 

echo "${green}${bold}           ### Phasing ended, the final graph is $(tput blink)${underline}graph_final.gfa${no_underline}${green}${reset} ${green}${bold} stats are available in $(tput blink)${underline}distributions.png${no_underline}${reset}"
#, stats are available in $(tput blink)${underline}distributions.png${no_underline.${reset}. "
