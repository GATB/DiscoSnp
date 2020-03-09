# Directory containing scripts called by [run_discoSnpRad.sh](https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/run_discoSnpRad.sh)



These scripts are not intended to be used by a user. They are called by [run_discoSnpRad.sh](https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/run_discoSnpRad.sh). 

1. **script** `discoRAD_clustering.sh`
   * This script 
     * manages bubble clustering from a discofile.fa file, 
     * manages the integration of cluster informations in a vcf file
     * remove variants with more than 95% missing genotypes 
     * remove low rank (<0.4) variants
   * This script is automatically used in [run_discoSnpRad.sh](https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/run_discoSnpRad.sh)
   * It presents no interest to be used alone
2. **script** ` from_SRC_to_edges.py `
   * Parses the short read connector output to a format readable by the clustering algorithm. 
   * This script is used in [discoRAD_clustering.sh](https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/clustering_scripts/discoRAD_clustering.sh)
   * It presents no interest to be used alone
3. **script** `fasta_and_cluster_to_filtered_vcf`
   * This script is used in [discoRAD_clustering.sh](https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/clustering_scripts/discoRAD_clustering.sh)
   * It presents no interest to be used alone



