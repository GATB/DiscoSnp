# FIND THE PATH CONTAINING THE SCRIPT: 
EDIR=$( python -c "import os.path; print(os.path.dirname(os.path.realpath(\"${BASH_SOURCE[0]}\")))" ) 

previous=$1



# removing a ghost from a line may create an empty line or an orphan SNP in a line. Thus we need an iterative process that removes ghost SNPs until nothing is done
while : ; do
    python3 ${EDIR}/remove_ghost_phased_SNPs.py $previous > afac_next.txt
    
    diff $previous afac_next.txt > /dev/null
    if [ $? -eq 0 ] ; then
        break
    fi
    
    cp afac_next.txt previous.txt
    previous="previous.txt"
done
rm -f previous.txt


# after removing ghost SNPs, some lines become equal, thus we have to reduce this redundancy. A sort -u is not sufficient as one need to cumulate the coverages. 
s=$1
basename=${s##*/}
sort afac_next.txt > afac_sorted.txt
rm -f afac_next.txt
python3 ${EDIR}/reduce_phased_allelle_redundancy.py afac_sorted.txt 
rm -f afac_sorted.txt

