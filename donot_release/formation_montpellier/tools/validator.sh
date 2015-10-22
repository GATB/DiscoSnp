prefix=$(basename $1 .fa)

EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/

ln -s ${EDIR}makeRoc.R .
ln -s ${EDIR}statistics.R .

echo python ${EDIR}remove_extensions_disco_file.py $prefix.fa  $prefix\_up.fa
python ${EDIR}remove_extensions_disco_file.py $prefix.fa  $prefix\_up.fa

# Mapping sur le génome de ref.
echo ${EDIR}gassst/Gassst -i ${prefix}_up.fa -d ${EDIR}humch1_first_5M.fasta -p 95 -w 15 -m 1 -l 0 -r 1 -s 5 -o $prefix.gassst
${EDIR}gassst/Gassst -i ${prefix}_up.fa -d ${EDIR}humch1_first_5M.fasta -p 95 -w 15 -m 1 -l 0 -r 1 -s 5 -o $prefix.gassst

# Analyse des données mappées
echo python ${EDIR}parse_gassst_discoMore_results_from_generic_position_list_creat_roc.py -d ${prefix}_up.fa -l ${EDIR}ref_human -g $prefix.gassst -r > $prefix.log
python ${EDIR}parse_gassst_discoMore_results_from_generic_position_list_creat_roc.py -d ${prefix}_up.fa -l ${EDIR}ref_human -g $prefix.gassst -r > $prefix.log

# Création des courbes ROC
echo Rscript makeRoc.R $prefix.png
Rscript makeRoc.R $prefix.png



# Creation des stats
echo python ${EDIR}parse_gassst_discoMore_results_from_generic_position_list_creat_roc.py -d ${prefix}_up.fa -l ${EDIR}ref_human -g $prefix.gassst -rv | grep FP | cut -d " " -f 2- > $prefix\_fp.txt 
python ${EDIR}parse_gassst_discoMore_results_from_generic_position_list_creat_roc.py -d ${prefix}_up.fa -l ${EDIR}ref_human -g $prefix.gassst -rv | grep FP | cut -d " " -f 2- > $prefix\_fp.txt 

echo grep "higher" $prefix.fa > $prefix.higher_headers
grep "higher" $prefix.fa > $prefix.higher_headers

echo Rscript statistics.R $prefix.higher_headers $prefix\_fp.txt $prefix\_stat.png
Rscript statistics.R $prefix.higher_headers $prefix\_fp.txt $prefix\_stat.png > $prefix\_stat.txt


# menage
rm -f $prefix.gassst roc_SNP roc_INDEL $prefix\_up.fa $prefix\_fp.txt  $prefix.higher_headers


