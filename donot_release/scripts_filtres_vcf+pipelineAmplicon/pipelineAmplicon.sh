#!/bin/bash
#V30112015 
#bash ../Script/pipelineAmplicon.sh <fof.txt> <prefix_disco> <individual_name> <b> <c> <MAF>
#bash ../Script/pipelineAmplicon.sh fofTPCHOR1c.txt TPCHOR1cWithRefOR TPCHOR1c 2 10 0.01
#Run pipeline on each sample
#for i in $(ls /home/Workhan_Data/Expes_Run_15S42A_INCA/RUN_15S42A/|cut -f1 -d "-" |sort|uniq);do bash ../Script/pipelineAmplicon.sh fof$i\.txt $i\WithRefOR $i 2 10 0.01 > log$i\MAF_001_b2_c10.txt ; done
fof=$1
prefix=$2
nameInd=$3

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DIRDISCO=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && cd .. && pwd )

#PATHS 
disco_path=$DIRDISCO
compare_vcf_path="$DIR/compare_vcf_Vamplicon.py"
sub_vcf_to_full_vcf_path="$DIR/sub_vcf_to_full_vcf.py"
filtre_uniq_pos="$DIR/filterOnBestDP_multiple_variant_at_same_pos.py"


#If you want to change the filter for discoSnp filter : put filter_MAF_path="$DIRDISCO/filter_out_using_MAF.py" and name_filter="filter_out_using_MAF.py"
filter_MAF_path="$DIR/filter_out_using_MAF_up_MAF_low_Amplicon_VCF.py"
name_filter="filter_out_using_MAF_up_MAF_low_Amplicon_VCF.py"

######################################################################
#TODO Put your own paths
bwa_path="" #"-B /home/Workhan_Data/Documents/Programmes/bwa-0.7.10/"
fofFile_path="/Users/ppeterlo/workspace/gatb-tools/gatb-tools/tools/DiscoSNP/donot_release/scripts_filtres_vcf+pipelineAmplicon/fof_amplicons/" #/home/Workhan_Data/Expes_Run_15S42A_INCA/fofFiles/"
truth_File_path="/Users/ppeterlo/workspace/gatb-tools/gatb-tools/tools/DiscoSNP/donot_release/clefs_chloe_1_dec_2015/truth_files/"
G="-G  /Users/ppeterlo/workspace/gatb-tools/gatb-tools/tools/DiscoSNP/donot_release/scripts_filtres_vcf+pipelineAmplicon/amplicons_hg19.fa"
######################################################################

#DiscoSNP++ options
D=100
P=3
R="-R"
k=31
b=$4
c=$5
if [ -n "$6" ] ; then
	MAF=$6
else 
	MAF=""
fi

#To do it on the cluster do not forget :
#source /local/env/envgcc-4.9.1.sh
#source /local/env/envbwa-0.7.10.sh


echo "#$filter_MAF_path"
echo "#$name_filter"
echo "#sub_vcf_to_full_vcf.py"
echo "#$compare_vcf_path"
echo "MAF : $MAF"

#Disco
$disco_path\/run_discoSnp++.sh -r $fofFile_path$fof -c $c -D $D -P $P -p $prefix -b $b $R $G -k $k $bwa_path -g -d 2



#-------------------------------------------------------------------------------------------------------------------------------------
#filter on MAF : DiscoSnp++ version

if [ "$name_filter" == "filter_out_using_MAF.py" ]; then
        if [ "$MAF" != ""  ]; then #filter on the MAF in the fasta file 
	        python $filter_MAF_path $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent.fa $MAF > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_MAF.fa

        else # With no filtering on MAF 
	        cat $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent.fa > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_MAF.fa
        fi
        rm $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_*.vcf #rm previous VCF file 

        #Run VCF_creator => create new vcf file after filtering 
        $disco_path\run_VCF_creator.sh $G -p $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_MAF.fa -o $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent.vcf -I $bwa_path

        #Conversion positions into real positions on the reference genome 
        python $sub_vcf_to_full_vcf_path $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV.vcf > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations.vcf
        
        cat $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations.vcf |grep "PASS" > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locationsPASS.vcf

#---------------------------------------------------------------------------------------------------------------------------------------------
#filter on MAF amplicon version
elif [ "$name_filter" == "filter_out_using_MAF_up_MAF_low_Amplicon_VCF.py" ]; then
#Conversion positions into real positions on the reference genome 
        python $sub_vcf_to_full_vcf_path $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV.vcf > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations.vcf

        if [ "$MAF" != ""  ]; then
                python $filter_MAF_path $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations.vcf $MAF > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations_MAF.vcf

        else
                cat $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations.vcf > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations_MAF.vcf
        fi
        cat $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locations_MAF.vcf |grep "PASS" > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locationsPASS.vcf
fi


#---------------------------------------------------------------------------------------------------------------------------------------------
#Results

#filter uniq position on best DP
python $filtre_uniq_pos $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locationsPASS.vcf > $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locationsPASSuniq.vcf

python $compare_vcf_path $truth_File_path$nameInd\_controle_positifs.csv $prefix\_k_$k\_c_$c\_D_$D\_P_$P\_b_$b\_coherent_for_IGV_real_genomic_locationsPASSuniq.vcf 1









