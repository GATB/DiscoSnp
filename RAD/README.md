# DiscoSnpRAD: small variant discovery and genotyping for RAD-seq data


DiscoSnpRAD is a pipeline based on discoSNP++ to discover small variants in RAD-like sequencing data. The differences with respect to using directly discoSNP++ lies in three main features:   
* an enhanced bubble model to deal with RAD-like sequences 
* using specific discoSNP++ parameters and filters, adapted to RAD-like data
* clustering the called variants into loci

**Reference:**   
Gauthier, J., Mouden, C.,  Suchan, T., Alvarez, N., Arrigo, N., Riou, C., Lemaitre, C., Peterlongo, P. (2017). [DiscoSnp-RAD: de novo detection of small variants for population genomics](https://www.biorxiv.org/content/early/2017/11/09/216747). BioRxiv

## Installation

* discoSNP++
* `short_read_connector` must have been downloaded and installed (clustering task). [https://github.com/GATB/short_read_connector](https://github.com/GATB/short_read_connector)


## Usage

```
./run_discoSnpRad.sh --fof read_file_of_files --src_path <directory> [discoSnp++ OPTIONS]
```

Clustering
```
-S|--src_path <directory>
    **absolute** path to short_read_connector directory, containing the "short_read_connector.sh" file. 
    -Note1: short read connector must be compiled.
    -Note2: with this option, discoSnpRad provide a vcf file containing SNPs and INDELS, clustered by locus
```

All other options are described in [discoSnp++ README](../README.md). Note that many discoSNP++ parameters have here default values, specifically adapted to RAD-seq data.

To see all options:
```
./run_discoSnpRad.sh -h
```


## Output

* a log file reminds all filtering steps applied and the name of the output .vcf file
* a vcf file containing results of filtering and clustering


