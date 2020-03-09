# DiscoSnp-RAD cookbook

{:toc}

Datasets `set1.fastq.gz`, `set2.fastq.gz`, ..., `set5.fastq.gz`  for this cookbook can be downloaded as following: 

```bash
for i in 1 2 3 4 5; do wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/set$i.fastq.gz; done
```

For each example, the full commands are proposed at the end of the section. 

## 1. No reference genome - Using only reads 1 from pairs - No clustering

This is the most classical usage of DiscoSnp-RAD. 

Consider one has $n$ rad datasets composed only of reads 1. Suppose those sets are called `set1.fastq.gz`, `set2.fastq.gz`, ..., `setn.fastq.gz`.

**First**, one need to create a *file of file* (fof), used as input for DiscoSnp-RAD. In this case, we want each data set to be considered individually in the displayed results. Hence, each line of the input fof contains exactly one dataset name. We may create this fof as following: 

```bash
ls set*.fastq.gz > my_fof.txt
```

> **Reminder**
>
> DiscoSnp-RAD can analyse fastq or fasta files, gzipped or not.

> **Note**
>
> If you wish to run DiscoSnp-RAD from a directory distinct from the one containing the read datasets, you have to indicate the absolute paths in the fof: 
>
> ```bash
> ls -d $PWD/set*.fastq.gz > my_fof.txt
> ```

**Second**, run discoSnp-RAD, using the input fof:

```bash
/my/discoSnp/path/discoSnpRAD/run_discoSnpRad.sh -r my_fof.txt 
```

That's it! Results are available in files 

* `discoRad_k_31_c_3_D_0_P_5_m_5_raw.fa` that contains the raw bubble sequences detected by discoSnp. For instance the first variant of this file is :

```
>SNP_higher_path_9965|P_1:30_C/G|high|nb_pol_1|left_unitig_length_83|right_unitig_length_6|left_contig_length_83|right_contig_length_6|C1_0|C2_22|C3_0|C4_0|C5_0|Q1_0|Q2_71|Q3_0|Q4_0|Q5_0|G1_1/1:364,58,5|G2_0/0:5,70,444|G3_1/1:384,61,5|
G4_1/1:424,67,5|G5_1/1:364,58,5|rank_1
atgtggcctgccgaggtggaggcggtcatcgacgagctgccggaggtgaagcgggtgtgcgtgatcggggtttacgacgagacCCAGGGAGATGTGCCTGGTGCCCTGGTTGTCCGGGAGGATAATGCCACTCTGACCGCACAGcaggtg
>SNP_lower_path_9965|P_1:30_C/G|high|nb_pol_1|left_unitig_length_83|right_unitig_length_6|left_contig_length_83|right_contig_length_6|C1_18|C2_0|C3_19|C4_21|C5_18|Q1_71|Q2_0|Q3_71|Q4_71|Q5_71|G1_1/1:364,58,5|G2_0/0:5,70,444|G3_1/1:384,
61,5|G4_1/1:424,67,5|G5_1/1:364,58,5|rank_1
atgtggcctgccgaggtggaggcggtcatcgacgagctgccggaggtgaagcgggtgtgcgtgatcggggtttacgacgagacCCAGGGAGATGTGCCTGGTGCCCTGGTTGTGCGGGAGGATAATGCCACTCTGACCGCACAGcaggtg
```

​	It is a SNP for which one can find **1/** its id: `9965`, **2/** its variable nucleotides `A/C` **3/** its nucleotidic complexity `high` (this is not a low complexity sequence), **4/** the length of its left and right unitigs `83` and `6` (corresponding resp. to the leftmost and rightmost lowercase acgt chatacters). For the higher path (first sequence corresponding to the `C` allele - the lower path corresponding to the `G`allele) **5/**  the abundance of the allele is provided in the Ci fields (`0` for first dataset, `22` for the second and so on), **6/** the average phred quality of mapped reads on the variable locus is provided (`0`for first dataset, `71` for the second and so on). For the both paths, **7/** the estimated genotype of the variant is provided for each dataset (`G1` to `G5`), and **8/** its rank (see manuscript for detailed explanation). 

​	<u>Note1</u>: information **1/**, **2/**, **3/**, **4/**, **7/** and **8/** are redundant for the two lines, while information **5/** and **6/** are specific to each allele, that is to say to each line. 

​	<u>Note2</u>: higher case characters correspond to the sequences of the bubble, and lower case sequences correspond to the left and right unitigs surrounding the bubble. 

* `discoRad_k_31_c_3_D_0_P_5_m_5.vcf`. For instance the first variant is this file is 

```
SNP_higher_path_9965	113	9965	C	G	.	.	Ty=SNP;Rk=1.0;UL=83;UR=6;CL=83;CR=6;Genome=.;Sd=.;Cluster=.;ClSize=.	GT:DP:PL:AD:HQ	1/1:18:364,58,5:0,18:0,71	0/0:22:5,70,444:22,0:71,0	1/1:19:384,
61,5:0,19:0,71	1/1:21:424,67,5:0,21:0,71	1/1:18:364,58,5:0,18:0,71
```

​	This line contains the same information as the one existing in the comment of the fasta file previously presented. Additional fields are proposed. They are empty as they are related to the usage of a reference genome that we have not done in this example or to the clustering of the variants that has not been done neither. 

* `discoRad_read_files_correspondance.txt`. This files simply recall the correspondence between the `C1`, ... `C5` fields with their datasets names. For instance, first line is `C_1 set1.fastq.gz`



**Full commands:**

```bash
disco_path=/my/discoSnp/path/	
```

```
for i in 1 2 3 4 5; do wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/set$i.fastq.gz; done
ls set*.fastq.gz > my_fof.txt
${disco_path}/discoSnpRAD/run_discoSnpRad.sh -r my_fof.txt
```



## 2. No reference genome - Using only reads 1 from pairs - With clustering

Given the fof file created as previously, run discoSnp-RAD indicating the `short_read_connector` installation path. 

 ```bash
/my/discoSnp/path/discoSnpRAD/run_discoSnpRad.sh -r my_fof.txt -S /my/SRC/path/
 ```

​	In this case we retrieve the `discoRad_k_31_c_3_D_0_P_5_m_5_raw.fa` and the `discoRad_read_files_correspondance.txt` files as previously exposed.

​	The vcf file is now called `discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf` as the `Cluster` field contains for each variant the cluster id of the variant and the `ClSize` field contains the size of the corresponding cluster.

​	In addition a second .fa file is provided: `discoRad_k_31_c_3_D_0_P_5_m_5_raw_filtered.fa`. In this file,  variants with more than 0.4 missing data and rank<0.4 are filtered out. 

**Full commands:**

```bash
disco_path=/my/discoSnp/path/	
src_path=/my/short_short_read_connector/path
```

```
for i in 1 2 3 4 5; do wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/set$i.fastq.gz; done
ls set*.fastq.gz > my_fof.txt
${disco_path}/discoSnpRAD/run_discoSnpRad.sh -r my_fof.txt -S ${src_path}
```



## 3. Using a reference genome - Using only reads 1 from pairs - With clustering

If one disposes for a reference genome, it can be used for determine the position of each predicted variant (without using the reference genome for prediction) on the genome. 

We may use the proposed reference genome as following:

```bash
wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/ref.fa
```

Variants from `discoRad_k_31_c_3_D_0_P_5_m_5_raw_filtered.fa` can be mapped to the reference genome as following:

```bash
/my/discoSnp/path/scripts/run_VCF_creator.sh -G ref.fa -p discoRad_k_31_c_3_D_0_P_5_m_5_raw_filtered.fa -e -o temp.vcf
```

<u>Note</u>: bwa must be installed and in the path. 

In order to retrieve the clustering information already computed, the `add_cluster_info_to_mapped_vcf.py `script may be used:

```
python /my/discoSnp/path/discoSnpRAD/post-processing_scripts/add_cluster_info_to_mapped_vcf.py -m temp.vcf -u discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf -o discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf
```

**Full commands:**

```bash
disco_path=/my/discoSnp/path/	
src_path=/my/short_short_read_connector/path
```

```
for i in 1 2 3 4 5; do wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/set$i.fastq.gz; done
wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/ref.fa
ls set*.fastq.gz > my_fof.txt
${disco_path}/discoSnpRAD/run_discoSnpRad.sh -r my_fof.txt -S ${src_path}
${disco_path}/scripts/run_VCF_creator.sh -G ref.fa -p discoRad_k_31_c_3_D_0_P_5_m_5_raw_filtered.fa -e -o temp.vcf
python ${disco_path}/discoSnpRAD/post-processing_scripts/add_cluster_info_to_mapped_vcf.py -m temp.vcf -u discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf -o discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf
```



## 4. Using forward and reverse reads. 

Whatever the wanted usage (with or without reference genome, with or without clustering) one may use pairend data. 

Imagine one disposes from 5 pairend read sets

* `set1_1.fastq.gz`, `set1_2.fastq.gz`
* `set2_1.fastq.gz`, `set2_2.fastq.gz` 
* ...
* `set5_1.fastq.gz`, `set5_2.fastq.gz`



1. **If our aim is to consider individually each file** (hence considering each file as a set), then we can simply create a *file of files* (fof) in which each line is a .fastq.gz file:

```bash
ls *.fastq.gz > my_fof.txt
```

2. **If our aim is to virtually concatenate forward and reverse reads for each sample**, we have to create as many fof as samples: 

```bash
for (( i=1; i<=5; i++ ));
	do 
	ls set${i}_*.fastq.gz > my_fof_set${i}.txt
	; 
done
```

​	Finally the fof provided to `run_discoSnpRad.sh` script (`-r` option) is a file in which each line is a fof file for one sample:

```bash
ls my_fof_set*.txt > my_fof.txt
```

