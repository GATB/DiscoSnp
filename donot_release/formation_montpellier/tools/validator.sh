prefix=$(basename $1 .fa)"_valid"

# Mapping sur le génome de ref.
~/workspace/gassst/trunk/Gassst -i $1 -d humch1_first_10M.fasta -p 95 -w 15 -m 1 -l 0 -r 1 -s 5 -o $prefix.gassst

# Analyse des données mappées
python ./parse_gassst_discoMore_results_from_generic_position_list_creat_roc.py -d $1 -l ref_human -g $prefix.gassst -r > $prefix.log

# Création des courbes ROC
Rscript makeRoc.R $prefix.png


rm -f $prefix.gassst roc_SNP roc_INDEL


