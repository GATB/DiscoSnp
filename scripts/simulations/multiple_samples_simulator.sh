#!/bin/bash

#multiple_samples_simulator.sh
#
#Bash script to simulate sequencing data from a reference genome and following various divergence statistics
#
#
#
###############################################################################
# set default options
###############################################################################

num_pop=1
num_sample=20
div_pop=0.015
#percentage of SNP shared by all samples from the same population
perc_sha=70
#percentage of shared SNP by samples
perc_spe=60
#percentage of homozygosity (1-x heterozygosity)
perc_H=90

#read lenght
read_l=150
#number of reads simulated
read_s=876960

###############################################################################
# parse options
###############################################################################


while getopts "g:p:s:d:p:n:o:l:x:h" OPTION; do
		case $OPTION in

				g) genome=${OPTARG} ;;
				p) num_pop=${OPTARG} ;;
				s) num_sample=${OPTARG} ;;
				d) div_pop=${OPTARG} ;;
				m) perc_sha=${OPTARG} ;;
				n) gperc_spe=${OPTARG} ;;
				o) perc_H=${OPTARG} ;;
				l) read_l=${OPTARG} ;;
				x) read_s=${OPTARG} ;;
				h)
						echo "Usage:"
						echo ""
						echo "   -g	 fasta of the genome"
						echo "   -p	 number of population to simulate"
						echo "   -s	 number of sample by population to simulate"	
						echo "   -d	 population divergence from reference genome"	
						echo "   -m	 percentage of population specific polymorphism"	
						echo "   -n	 percentage of sample specific polymorphism"	
						echo "   -o	 percentage of homozygosity (1-x heterozygosity)"
						echo "   -l	 read lenght"
						echo "   -x	 number of reads simulated"	 #warning must be multiplied by 2			
						echo "   -h	 help"
						exit 0
						;;
		esac
done

if [ -z ${genome} ]
then
		echo "missing genome"
		exit 0
fi

for p in `seq 1 $num_pop`
	do
	python ./random_mut_fasta.py $genome $div_pop > ERASEME_pos_mut_pop"$p"

	#pos random ordering
	sort -R ERASEME_pos_mut_pop"$p" > ERASEME_pos_mut_random_pop"$p"
	#nb all mut
	nb_line_all=`grep "." -c ERASEME_pos_mut_random_pop"$p"`
	#nb_mut_pop_specific
	nb_shared_all=`echo $(($nb_line_all/100*$perc_sha ))`
	head -n +"$nb_shared_all" ERASEME_pos_mut_random_pop"$p" > ERASEME_pos_mut_random_shared_allpop"$p"
	nb_spe=`echo $(($nb_shared_all+1 ))`
	tail -n +"$nb_spe" ERASEME_pos_mut_random_pop"$p" > ERASEME_pos_mut_random_to_separate_by_pop"$p"


	for i in `seq 1 $num_sample`
		do
		nb_line=`grep "." -c ERASEME_pos_mut_random_to_separate_by_pop"$p"`
		nb_shared=`echo $(($nb_line/100*$perc_spe ))`
		#selection rondomly 80% of the mutations for each of the 20 samples
		sort -R ERASEME_pos_mut_random_to_separate_by_pop"$p" | head -n +"$nb_shared" > ERASEME_pos_mut_random_shared_pop"$p"_s"$i"
		#homozygotes mutations
		cat ERASEME_pos_mut_random_shared_allpop"$p" ERASEME_pos_mut_random_shared_pop"$p"_s"$i" > ERASEME_pos_mut_random_spe_pop"$p"_and_s"$i"
		#mutation inducing
		python ./targeted_mut_fasta_corrected.py "$genome" ERASEME_pos_mut_random_spe_pop"$p"_and_s"$i"
		mv "$genome"_mut ERASEME_"$genome"_pop"$p"_s"$i"_withhetero.fasta
		#homozygote and heterozygote mutations
		nb_line2=`grep "." -c ERASEME_pos_mut_random_shared_pop"$p"_s"$i"`
		nb_homo=`echo $(($nb_line2/100*$perc_H ))`
		#homo mutations
		head -n +"$nb_homo" ERASEME_pos_mut_random_shared_pop"$p"_s"$i" > ERASEME_pos_mut_random_shared_pop"$p"_s"$i"_homo
		#population mutationq
		cat ERASEME_pos_mut_random_shared_allpop"$p" ERASEME_pos_mut_random_shared_pop"$p"_s"$i"_homo > ERASEME_pos_mut_random_spe_pop"$p"_and_s"$i"_homo
		python ./targeted_mut_fasta_corrected.py "$genome" ERASEME_pos_mut_random_spe_pop"$p"_and_s"$i"_homo
		mv "$genome"_mut ERASEME_"$genome"_pop"$p"_s"$i"_homo.fasta
		#READS SIMULATION
		mutareads_forward ERASEME_"$genome"_pop"$p"_s"$i"_withhetero.fasta pop"$p"_ind"$i"_allele1_err_reads $read_s $read_l 0.01 0 0
		mutareads_forward ERASEME_"$genome"_pop"$p"_s"$i"_homo.fasta pop"$p"_ind"$i"_allele2_err_reads $read_s $read_l 0.01 0 0
		cat pop"$p"_ind"$i"_allele{1..2}_err_reads.fasta > pop"$p"_ind"$i"_err_reads.fasta
		rm pop"$p"_ind"$i"_allele1_err_reads.fasta
		rm pop"$p"_ind"$i"_allele2_err_reads.fasta
		done
	done
#ghost only position vcf creator
cat ERASEME_pos_mut_random_shared_allpop* ERASEME_pos_mut_random_spe_pop*_and_s* | awk ' !x[$0]++' > only_position.vcf
#clean
rm ERASEME_*
rm genome_mut.fasta
rm genome_ref.fasta


