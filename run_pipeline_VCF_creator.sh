#!/bin/bash


DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
remove=1
PATH_VCF_creator=""
samfile=""
vcffile=""
genome=""
PATH_BWA=""
discoSNPs=""
l=""
n=""
s=""


function help {
echo " ##############################"
echo "   Run VCF_creator pipeline     "
echo " ##############################"
echo "Usage : ./run_pipeline_VCF_creator.sh OPT"
echo -e "##Ghost mode : Create a vcf file without alignment : ./run_pipeline_VCF_creator.sh -p <disco_file>.fa -o <output>.vcf"
echo -e "##Pipeline mode : Alignment +VCF Creation : ./run_pipeline_VCF_creator.sh -b <path_bwa> -g <ref>.fasta -p <disco_file>.fa -o <output>.vcf -n <mismatch_number>"
echo -e "##VCF creation : samfile already there, goes to : ./run_pipeline_VCF_creator.sh -f <file>.sam -n <mismatch_number> -o <output>.vcf"
echo -e "\t-h : print this message"
echo -e "\t-b : path where bwa is"
echo -e "\t-c : path where VCF_creator is"
echo -e "\t-g : reference genome (.fasta) for alignment with snps from discosnp++ with path"
echo -e "\t-p : snps from discosnp++ with path"
echo -e "\t-s : bwa option : distance with the seed"
echo -e "\t-n : bwa option : allowed distance with the reference for a snp"
echo -e "\t-l : bwa option : lenght of the seed for alignment"
echo -e "\t-f : <file>.sam : skip the algnment phases to create the vcf file "
echo -e "\t-o : output <file>.vcf"
echo -e "\t-w : remove tmp files (index files)"

}

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

while getopts "hb:c:g:p:l:n:s:wf:o:" opt; do
case $opt in

	w)
	remove=0
	;;

	h)
	help
	exit 1
	;;

	b)
	echo -e "\t##BWA directory: $OPTARG" >&2
	PATH_BWA=$OPTARG
	;;

	c)
	echo -e "\t##VCF_creator directory: $OPTARG" >&2
	PATH_VCF_creator=$OPTARG
	;;

	g)
	echo -e "\t##use genome : $OPTARG" >&2
	genome=$OPTARG
	;;

	p)
	echo -e "\t##use disco SNPS : $OPTARG" >&2
	discoSNPs=$OPTARG
	;;

	s)
	echo -e "\t##use distance with the seed : $OPTARG" >&2
	l=$OPTARG
	;;

	n)
	echo -e "\##use number of mismatches : $OPTARG" >&2
	n=$OPTARG
	;;

	l)
	echo -e "\t##use seed length : $OPTARG" >&2
	s=$OPTARG
	;;

	f)
	echo -e "\t##use directly samfile : $OPTARG" >&2
	samfile=$OPTARG
	;;
		
	o)
	echo -e "\t##output : $OPTARG" >&2
	vcffile=$OPTARG
	;;
	
	\?)
	echo -e "##Invalid option: -$OPTARG" >&2
	exit 1
	;;

	:)
	echo "##Option -$OPTARG requires an argument." >&2
	exit 1
	;;
esac
done
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
###Tests
if [ -z "$PATH_VCF_creator" ];then
	PATH_VCF_creator=$DIR"/tools"
fi
if [ ! -e  $PATH_VCF_creator/VCF_creator.py ]; then
	echo "...Unable to find VCF_creator..."
	exit 1
fi
if [ -z "$s" ];then
	s=0
fi
if [ -z "$l" ];then
	l=10
fi

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
####Use the pipeline of aligement
if [ -z "$samfile" ];then
        #Ghost mode
        if [ -z "$genome" ]; then
                echo -e "...Ghost mode..."
                echo -e "...Creation of a vcf without alignment..."
                if [ -z "$discoSNPs" ] && [ -z "$vcffile" ];then
                       echo -e "...To create a vcf without alignment ..."
                       echo -e "...You must provide an output <file>.vcf : option -o..."
                       echo -e "...And the file disco : option -p..."
                       exit 1 
                else
                        n=3
                        vcf=$(basename $vcffile .vcf)"_"$(basename $discoSNPs .fa)".vcf"
                        python $PATH_VCF_creator/VCF_creator.py -s $discoSNPs -n $n -o $vcf
		        echo -e "... Creation of the vcf file : done ...==> $vcf" 
		        exit 
                fi    
        fi
	if [ -z "$PATH_BWA" ] ;then
		IS_BWA=$(command -v bwa)
	
	
	
		if [ -z "$IS_BWA" ];then
			echo -e "... BWA not found... add bwa to \$PATH or give directly the path (-b)"
			exit 1
		else 
			PATH_BWA=$(dirname $IS_BWA)
		fi
	fi
	
	
	if [ -z "$vcffile" ] || [[ "$vcffile" =~ *.vcf ]]; then
		echo -e "...You must provide an output <file>.vcf : option -o (for help -h)..."
		echo -e "...Usage : ./run_pipeline_VCF_creator.sh OPT..."
		exit 1
	fi
	if [ -z "$genome" ]; then
		echo -e "...You must provide a genome of reference : option -g (for help -h)..."
		exit 1
	fi
	
	if [ -z "$discoSNPs" ];then
		echo "... Error : file disco is missing : option -p (for help -h)..."
		exit 1
	else
		if [ ! -e $PATH_VCF_creator/remove_extensions_disco_file.py ];then
			echo "...Unable to find remove_extensions_disco_file.py..."
			exit 1
		else
			discoSNPsbis=$(basename $discoSNPs .fa)"bis.fasta"
			python $PATH_VCF_creator/remove_extensions_disco_file.py $discoSNPs $discoSNPsbis
			if [ -z "$discoSNPsbis" ];then
				echo "...Error with the script remove_extensions_disco_file.py..."
			fi
		fi
	fi
	if [ -z "$n" ];then
		echo -e "\t##Default value for the number of mismatches allowed in alignment : 3 (to change it -n)"
		n=3
	fi
###BWA 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
	#BWA files
	vcf=$(basename $vcffile .vcf)"_"$(basename $discoSNPs .fa)"_n"$n"_l"$l"_s"$s".vcf"
	samfile=$(basename $vcffile .vcf)"_"$(basename $discoSNPs .fa)"_n"$n"_l"$l"_s"$s".sam"
	saifile=$(basename $vcffile .vcf)"_"$(basename $discoSNPs .fa)"_n"$n"_l"$l"_s"$s".sai"
	indexamb=$genome".amb"
	indexann=$genome".ann"
	indexbwt=$genome".bwt"
	indexpac=$genome".pac"
	indexsa=$genome".sa"
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
	##Indexation of reference genome
	if [ -e $indexamb ] && [ -e $indexann ] && [ -e $indexbwt ] && [ -e $indexpac ] && [ -e $indexsa ]; then
		echo -e "...Indexation : Using the existing index..."
	else
		$PATH_BWA/bwa index $genome
	fi
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
	##Alignment disco snps on the reference genome
	if [ -e $saifile ]; then
		echo -e "...Alignment : Using the existing file : $saifile"
	else
		$PATH_BWA/bwa aln -N -k $s -n $n -l $l $genome $discoSNPsbis > $saifile
	fi
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
	##Creation of the samfile
	if [ -e $samfile ]; then
		echo -e "...Samfile already exists. Skipped to vcf creation..."
	else
		$PATH_BWA/bwa samse -n 1000 $genome $saifile $discoSNPsbis > $samfile
	fi 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
	##Creation of the vcf file
	if [ -e $vcf ]; then
		echo -e "...VCF file for this alignment already exists...==> $vcf"
	else 
		python $PATH_VCF_creator/VCF_creator.py -s $samfile -n $n --o $vcf
		echo -e "... Creation of the vcf file : done ...==> $vcf"
	fi	
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
	
else
####Skip alignment phase to create a vcf file
##Test to execute VCF_creator
	echo -e "...Skipping alignment phase..."
	if [ -z "$vcffile" ] || [[ "$vcffile" =~ *.vcf ]]; then
		echo -e "...You must provide an output <file>.vcf..."
		exit 1
	fi
	if [ -z "$n" ]; then
		echo -e "...You must provide the number of mismatch allowed during the alignment phase.."
		exit 1
	fi
	##Creation of the vcf file
	python $PATH_VCF_creator/VCF_creator.py -s $samfile -n $n -output $vcffile	
fi

if [ $remove=0 ];then
	rm -r $saifile $discoSNPsbis
else
	rm -r $indexamb $indexann $indexbwt $indexpac $indexsa $saifile $discoSNPsbis
fi	


