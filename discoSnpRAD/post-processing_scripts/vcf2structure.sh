#! /bin/bash


out_name=`echo $1 | sed -e 's/.vcf//g'`
#extract genotypes
grep -v "#" $1 | sed -e 's/\t\.:/\t\.\/\.:/g' | cut -f 10- | sed -e 's/[:][[:graph:]]*//g' | sed -e 's/\//\t/g' -e 's/|/\t/g' -e 's/\./-9/g' > $1_matrix

#transpose matrix
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' $1_matrix > $out_name.str

rm $1_matrix


