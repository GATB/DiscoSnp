# DiscoSnpRAD: small variant discovery and genotyping for RAD-seq data


DiscoSnpRAD is a pipeline based on discoSnp++ to discover small variants in RAD-like sequencing data. The differences with respect to using directly discoSnp++ lies in three main features:   
* an enhanced bubble model to deal with RAD-like sequences 
* using specific discoSnp++ parameters and filters, adapted to RAD-like data
* clustering the called variants into loci

**Reference:**   
Gauthier, J., Mouden, C.,  Suchan, T., Alvarez, N., Arrigo, N., Riou, C., Lemaitre, C., Peterlongo, P. (2019). [DiscoSnp-RAD: de novo detection of small variants for population genomics](https://www.biorxiv.org/content/10.1101/216747v2). BioRxiv

​	Simulation and validation scripts: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3724518.svg)](https://doi.org/10.5281/zenodo.3724518)



## Installation

* DiscoSnpRAD is installed via the discoSnp++ basic instal (from sources or conda install. see [../README.md](../README.md))
* `short_read_connector` must have been downloaded and installed (clustering task). [https://github.com/GATB/short_read_connector](https://github.com/GATB/short_read_connector)


## Usage

```
./run_discoSnpRad.sh -r read_file_of_files -S -p myDiscoSnpRADresult [discoSnp++ OPTIONS]
```

**Clustering option** (RAD-specific option):

```
-S|--src [src_path]
    performs clustering of variants with short_read_connector
    src_path: **absolute** path to short_read_connector directory, containing the "short_read_connector.sh" file. 
    -Note1: short read connector must be compiled.
    -Note2: if no value is given, it assumes short_read_connector.sh is in the PATH env variable.
    -Note3: with this option, discoSnpRad outputs a vcf file containing the variants clustered by locus.
```

All other options are described in [discoSnp++ README](../README.md). Note that many discoSNP++ parameters have here default values, specifically adapted to RAD-seq data.

To see all options:
```
./run_discoSnpRad.sh -h
```


## Output

When run with output prefix name `myDiscoSnpRADresult`, the main output file is :

* `myDiscoSnpRADresult_[parameter_values]_clustered.vcf`: the final set of variants, with various information, including clustering per locus information (see VCF format below).
* or `myDiscoSnpRADresult_[parameter_values].vcf` if no clustering was performed.

Additionnally, several other files are output that can be usefull :

* `myDiscoSnpRADresult_[parameter_values]_raw.fa`: the raw set of variants in fasta format, prior to any filtering and clustering steps.
* `myDiscoSnpRADresult_[graph_parameter_values].h5`: the de Bruijn graph in h5 format (reusable with any GATB tool)
* `myDiscoSnpRADresult_read_files_correspondance.txt`: the correspondence between read file names and IDs given as genotypes in the vcf
* the standard output reminds all filtering steps applied and the name of the output .vcf file

#### VCF format

Each variant is described with: 

* an ID: `ID` column, 

* two alleles (`REF` and `ALT` columns), 

* a quality value: `INFO` column, `Rk`, between 0 (bad) and 1 (best),

* some clustering information: `INFO` field: with the locus id (`Cluster`) and its number of varying sites (`ClSize`),

* and for each sample in the genotype columns (`G1`, `G2`,...): the inferred genotype (`0/0`, `0/1`, `1/1`or `./.`for missing value), the read depths (`RD` total, `AD`per allele), among others.

  

## Cookbook

This [cookbook](./COOKBOOK.md) presents classical usages. 



## Content of this directory

Additionnally to the main script of discoSnpRAD, this directory contains two sub-directories :   
* [clustering_scripts](clustering_scripts/) : it contains the scripts used by the main script of discoSnpRAD for clustering and formatting the variants.
* [post-processing_scripts](post-processing_scripts/) : it contains several scripts that can be usefull to post-process the results of discoSnpRAD, ie. filtering results according to various criteria, changing format, preparing data for Structure, etc.


