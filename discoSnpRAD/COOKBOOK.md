# DiscoSnp-RAD cookbook

**Table of Contents**

* [1. No reference genome - Using only reads 1 from pairs - No clustering](#1)
* [2. No reference genome - Using only reads 1 from pairs - With clustering](#2)
* [3. Using a reference genome - Using only reads 1 from pairs - With clustering](#3)
* [4. Using forward and reverse reads.](#4)
* [5. Post-processing](#5)

- - - -

**Prerequisite**

Datasets `set1.fastq.gz`, `set2.fastq.gz`, ..., `set5.fastq.gz`  for this cookbook can be downloaded as following:

```bash
for i in 1 2 3 4 5; do wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/set$i.fastq.gz; done
```

**Full commands:** For each example, the full commands are proposed at the end of the section.

**Computation time** is approximately 10 to 15 minutes per example (once dataset is dowloaded).

***

**Note about multiplexed data**

To date, DiscoSnp-RAD is not able to consider multiplexed data. Hence, input files need to be demultiplexed to samples.

- - - -

## 1. <a name="1"> No reference genome - Using only reads 1 from pairs - No clustering</a>

This is the most classical usage of DiscoSnp-RAD.

Consider one has $n$ rad datasets composed only of reads 1. Suppose those sets are called `set1.fastq.gz`, `set2.fastq.gz`, ..., `setn.fastq.gz`.

**First**, one need to create a *file of file* (fof), used as input for DiscoSnp-RAD. In this case, we want each data set to be considered individually in the displayed results. Hence, each line of the input fof contains exactly one dataset name. We may create this fof as following:

```bash
ls set*.fastq.gz > my_fof.txt
```

> **Reminder** DiscoSnp-RAD can analyse fastq or fasta files, gzipped or not.  

> **Note** If you wish to run DiscoSnp-RAD from a directory distinct from the one containing the read datasets, you have to indicate the absolute paths in the fof:
> 
> 
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

​    It is a SNP for which one can find **1/** its id: `9965`, **2/** its variable nucleotides `A/C` **3/** its nucleotidic complexity `high` (this is not a low complexity sequence), **4/** the length of its left and right unitigs `83` and `6` (corresponding resp. to the leftmost and rightmost lowercase acgt chatacters). For the higher path (first sequence corresponding to the `C` allele - the lower path corresponding to the `G`allele) **5/**  the abundance of the allele is provided in the Ci fields (`0` for first dataset, `22` for the second and so on), **6/** the average phred quality of mapped reads on the variable locus is provided (`0`for first dataset, `71` for the second and so on). For the both paths, **7/** the estimated genotype of the variant is provided for each dataset (`G1` to `G5`), and **8/** its rank (see manuscript for detailed explanation).

​    Note1: information **1/**, **2/**, **3/**, **4/**, **7/** and **8/** are redundant for the two lines, while information **5/** and **6/** are specific to each allele, that is to say to each line.

​    Note2: higher case characters correspond to the sequences of the bubble, and lower case sequences correspond to the left and right unitigs surrounding the bubble.

* `discoRad_k_31_c_3_D_0_P_5_m_5.vcf`. For instance the first variant is this file is

```
SNP_higher_path_9965    113    9965    C    G    .    .    Ty=SNP;Rk=1.0;UL=83;UR=6;CL=83;CR=6;Genome=.;Sd=.;Cluster=.;ClSize=.    GT:DP:PL:AD:HQ    1/1:18:364,58,5:0,18:0,71    0/0:22:5,70,444:22,0:71,0    1/1:19:384,
61,5:0,19:0,71    1/1:21:424,67,5:0,21:0,71    1/1:18:364,58,5:0,18:0,71
```

​    This line contains the same information as the one existing in the comment of the fasta file previously presented. Additional fields are proposed. They are empty as they are related to the usage of a reference genome that we have not done in this example or to the clustering of the variants that has not been done neither.

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

## 2. <a name="2">No reference genome - Using only reads 1 from pairs - With clustering</a>

Given the fof file created as previously, run discoSnp-RAD indicating the `short_read_connector` installation path.

```bash
/my/discoSnp/path/discoSnpRAD/run_discoSnpRad.sh -r my_fof.txt -S /my/SRC/path/
```

​    In this case we retrieve the `discoRad_k_31_c_3_D_0_P_5_m_5_raw.fa` and the `discoRad_read_files_correspondance.txt` files as previously exposed.

​    The vcf file is now called `discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf` as the `Cluster` field contains for each variant the cluster id of the variant and the `ClSize` field contains the size of the corresponding cluster.

​    In addition a second .fa file is provided: `discoRad_k_31_c_3_D_0_P_5_m_5_raw_filtered.fa`. In this file,  variants with more than 0.4 missing data and rank<0.4 are filtered out.

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

## 3. <a name="3">Using a reference genome - Using only reads 1 from pairs - With clustering</a>

If one disposes for a reference genome, it can be used for determine the position of each predicted variant (without using the reference genome for prediction) on the genome.

We may use the proposed reference genome as following:

```bash
wget http://bioinformatique.rennes.inria.fr/data_cookbook_discoSnp-RAD/ref.fa
```

Variants from `discoRad_k_31_c_3_D_0_P_5_m_5_raw_filtered.fa` can be mapped to the reference genome as following:

```bash
/my/discoSnp/path/scripts/run_VCF_creator.sh -G ref.fa -p discoRad_k_31_c_3_D_0_P_5_m_5_raw_filtered.fa -e -o temp.vcf
```

Note: bwa must be installed and in the path.

In order to retrieve the clustering information already computed, the `add_cluster_info_to_mapped_vcf.py`script may be used:

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

## 4. <a name="4">Using forward and reverse reads.</a>

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

​    Finally the fof provided to `run_discoSnpRad.sh` script (`-r` option) is a file in which each line is a fof file for one sample:

```bash
ls my_fof_set*.txt > my_fof.txt
```

## 5. <a name="5">Post-processing</a>

A bench of post-processing scripts can be found in the dedicated [directory](https://github.com/GATB/DiscoSnp/tree/master/discoSnpRAD/post-processing_scripts).

We consider here that [case 3](#3.%20Using%20a%20reference%20genome%20-%20Using%20only%20reads%201%20from%20pairs%20-%20With%20clustering) (Using a reference genome - Using only reads 1 from pairs - With clustering) was performed, hence obtaining the file: `discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf`

### 5.1 Filtering scripts

#### **script** `filter_by_cluster_size_and_rank.py`

* Filtering on cluster size. 
  Use case: Need clusters of size between 2 to 100:
  
  ```bash
  python filter_by_cluster_size_and_rank.py -i discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf -o filtered_on_cluster_size.vcf -m 2 -M 100
  ```

* Filtering on rank.
  Use case: Need variants with a rank higher than 0.8:
  
  ```bash
  python filter_by_cluster_size_and_rank.py -i discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf -o filtered_on_rank.vcf -r 0.8
  ```

#### script `filter_vcf_by_indiv_cov_max_missing_and_maf.py`

* Filtering on read coverage.
  Use case: replace by `./.`original genotypes of variants whose total read coverage is below 20:
  
  ```bash
  python filter_vcf_by_indiv_cov_max_missing_and_maf.py -i discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf -o non_genotyped_low_covered.vcf -c 20
  ```

* Filtering on missing genotypes.
  Use case: Remove variants whose fraction of missing genotypes is greater than 70%:
  
  ```bash
  python filter_vcf_by_indiv_cov_max_missing_and_maf.py -i discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf -o filtered_on_missing_geno.vcf -m 0.7
  ```

* Filtering on minor allele frequency.
  Use case: Remove variants whose minor allele frequency is smaller than 0.2:
  
  ```bash
  python filter_vcf_by_indiv_cov_max_missing_and_maf.py -i discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf -o filtered_on_maf.vcf -f 0.2
  ```

#### script `filter_paralogs.py`

* Filtering on paralogs.
  Use case: Remove variants that belong to a cluster such that more than 50% of its variants  have each more than 10% of heterozygous genotypes:
  
  ```bash
  python filter_paralogs.py -i discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf -o filtered_on_paralogs.vcf -x 0.1 -y 0.5
  ```

### 5.2 Scripts for STRUCTURE analyses

#### script `1SNP_per_cluster.py`

* Conserve one variant per cluster (the one with less missing genotypes).
  Use case: prepare the vcf file to be used by STRUCTURE:
  
  ```bash
  python 1SNP_per_cluster.py -i discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf -o one_variant_per_cluster.vcf
  ```

#### script `vcf2structure.sh`

* Changes the vcf format to a Structure format (input of the software Structure).
  Use case: prepare the vcf file to be used by STRUCTURE:
  
  ```bash
  sh vcf2structure.sh discoRad_k_31_c_3_D_0_P_5_m_5_mapped.vcf > file.str
  ```

### 5.3 Mapping to a reference, and keeping the cluster information.

The script `add_cluster_info_to_mapped_vcf.py` can be used in case one wants to map predicted variants to any available genome. This use case was described in [case 3](#3.%20Using%20a%20reference%20genome%20-%20Using%20only%20reads%201%20from%20pairs%20-%20With%20clustering) (Using a reference genome - Using only reads 1 from pairs - With clustering) 
