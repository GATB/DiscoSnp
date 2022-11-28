## V2.6.2 
* Add information of phased SNPs in the VCF when using discoSnpRad

## V2.6.1

* Update python version in vcf creator, fixed CI issues

## V2.6.0

- VCF file is now one-based (previously zero-based)

## V2.4.5

* `-g` option takes a graph as parameter (instead of determining itself a graph name)
* update of gatb-core (compiles with gcc >= 7.3.0)
* Using xhash - kissreads is 2 times faster

## v2.4.4

* Fixes a bug in option parser for discoSnp-RAD
* Improves option presentations for discoSnp-RAD
* Includes a cookbook for discoSnp-RAD





## Older notes:

* Jan 2020
	+ Conda Install
	+ discoRad (independent / reorganized)
	+ phasing in test (hidden -A option)
	+ default values modified
* 9/6/2016 2.2.9
	+ Fixed a VCF creator bug
	+ Optimizing VCF creator
* 26/04/2016 2.2.8
	+ Fixed a tiny kissread bug
	+ Adding a continuous integration large test collection
* 15/04/2016 2.2.7 (now uses the githup
  repository: http://github.com/GATB/DiscoSnp)
     + Adding the possibility to limit the number of
       symmetrically branching crossroads traversed during
       the bubble finding
     + b2 mode: explore all possible symmetrical paths, even
       in case of success on one of the paths
     + b2 mode: avoid redundancies
     + Increased the maximal number of close SNP detectable
       thanks to a non recursive part on the bubble
       enumeration.
     + Fixed VCF creation bugs
     + Fixed Prediction bugs
     + Improved VCF creator error messages in case of missing
       values from BWA results
     + Increased the max breadth for the indel detection
* 24/11/2015 DiscoSNP++-2.2.4
     + Fixes a tiny bug in the run_discoSnp++.sh script.
       Thanks Hanan
       (https://www.biostars.org/p/155781/#167002)
* 18/11/2015 DiscoSNP++-2.2.3
     + Improvements:
          o Dump de read file names. C_1, C_2, … are provided
            in the .fa and the vcf file. Now a file indicates
            the correspondence between C_i and a set of read
            files. See documentation for details.
          o Removes indels if the repeat size is higher than
            a user defined threshold (max_ambigous_indel).
            Indels with Long repeat size (eg >20 in our
            tests) very often are false positives.
     + BUG correction:
          o VCF bug described here:
            https://www.biostars.org/p/166298/ is now
            fixed.
          o kissread bug when P > 1 (segfault corrected).
* 05/10/2015 DiscoSnp++-2.2.1
     + BUG correction:
          o Kissreads module time (bug fixed by Guillaume
            Rizk) and memory.
          o VCF creator bug correction
          o Redundant bubble detection suppression
     + VCF creators uses BWA MEM by default
* 17/07/2015: Important update – DiscoSNP++-2.2.0
     + Input read set format has changed. Use now file of
       files. This provides an easier way of dealing with
       read sets composed of several read files (pair end or
       pools). See the documentation in the doc directory
     + The kmer coverage threshold can be
          o set separately for each read set
          o and/or automatically detected
     + If a reference genome is provided, it can used for
       predicting variants. For instance a unique read set
       may be compared to the reference genome. (option -R)
     + With respect to previous change (-R) the read
       coherency in kissreads has changed.
          o Before: a variant was read coherent if its two
            path were read coherent
          o After: a variant is read coherent if at least one
            of its two paths is read coherent (else all
            homozygous calls obtains comparing a read set to
            a reference would be uncoherent)
     + Kissreads parallelization had been improved. OMP is
       not used anymore, and running time are decreased.
     + Two memory bugs have been fixed. They occurred mainly
       while using large number of read sets.
     + It is now possible to detect only indels (ie. -P 0
       detects no SNP)
     + Memory issue detected:
          o All tools (kissnp2, kissreads2, VCF_creator) have
            a limited memory footprint. However, due to
            kissnp parallelization, in some cases, memory may
            increase linearly with the number of used cores.
            If the memory is too high, limit the number of
            cores with the -u option.
* 13/05/2015: DiscoSNP++-2.1.7
     + Fixes a bug with very long reads
* 04/05/2015: DiscoSNP++-2.1.6 + (mac and linux)
  binaries: http://gatb.inria.fr/binaries-url/
     + Adding the vcf doc
     + Adding the -u option (limiting the maximal number of
       used threads)
     + Fixing some compilation bugs
* 03/05/2015: DiscoSNP++-2.1.5 +
     + Fixes compilation bugs with some compilers
     + Fixes some VCF generation bugs
* 02/04/2015: DiscoSNP++-2.1.4
     + Fixes a redundancy bug with the b 2 option
     + Generates an IGV compatible VCF
     + Fixes various small VCF bugs
* 23/03/2015: DiscoSNP++-2.1.3
     + Automatically creates a VCF
     + Genotyping of results (for diploïds)
     + Warning: documentation not up to date.
* 02/03/2015: DiscoSNPpp-2.0.6.
     + Genotypes are automatically computed
* 24/02/2015:DiscoSNPpp-2.0.5, fixes a compilation bug on
  some OS.
* 24/02/2015:DiscoSNPpp-2.0.4
     + Fasta headers more informatives
     + Update of the kissread tool for close SNPs precision
     + Unique id for SNPs and indels
     + Documentation update
* 26/01/2015: DiscoSNPpp-2.0.2
     + Fixed a bug progress bar display on macos.
* 22/01/2015: DiscoSNPpp-2.0.1
     + Fixed a bug about command line parser.
* 19/01/2015: DiscoSNPpp-2.0.0
     + Corrected the discoSnp++2csv.py file
     + Faster, nicer presentation, unique .h5 graph file
     + Detects also indels and close SNPs
* 29/10/2014: [BETA] discoSnp_1.2.6: Kissreads changes:
  added minimal spanning option – changed the  strategy when
  a read has multiple hits on a fragment.
* 01/10/2014: discoSnp_1.2.5: fixed bugs of version 1.2.4
  (bugs were not affecting version <1.2.4).
* 25/09/2014: discoSnp_1.2.4: increased kissread speed
  (approx x2) –– bugged version — don’t use
* 28/07/2014: discoSnp_1.2.3: improved the documentation
  – cleaned the output messages and the output useless files
* 20/05/2014: discoSnp_1.2.2: -b 0 is the default option.
* 26/04/2014: discoSnp_1.2.1: Fixes a compilation problem
  due to a makefile typo.
* 05/03/2014: discoSnp_1.2.0: New option with regards to
  branching bubbles
* 14/11/2013: discoSnp_1.0.1: fixes a bug concerning the
  contigs extensions (note that unitig extension lengths were
  not affected)
* 14/10/2013: discoSnp_1.0.0: read2SNPs name changes.
  discoSnp_1.0.0 is a new name for read2SNPs_2.1.1.4
