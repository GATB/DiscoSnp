# DiscoSnp

| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/DiscoSnp/job/tool-discosnp-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/DiscoSnp/job/tool-discosnp-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/DiscoSnp/job/tool-discosnp-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/DiscoSnp/job/tool-discosnp-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## What is DiscoSnp?

DiscoSnp is designed for discovering all kinds of SNPs (not only isolated ones),  as well as insertions and deletions, from raw set(s) of reads. The number of input read sets is not constrained, it can be one, two, or more. No reference genome is needed.

Uricaru R., Rizk G., Lacroix V., Quillery E., Plantard O., Chikhi R., Lemaitre C., Peterlongo P. (2014). [Reference-free detection of isolated SNPs](http://nar.oxfordjournals.org/content/43/2/e11). Nucleic Acids Research 43(2):e11.

## Getting the latest source code

### Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

### Instructions

    # get a local copy of DiscoSnp source code
    git clone --recursive https://github.com/GATB/DiscoSnp.git
    
    # compile the code an run a simple test on your computer
    cd DiscoSnp
    sh INSTALL

## Getting a binary stable release

Binary release for Linux and Mac OSX are provided within the "Releases" tab on Github/DiscoSnp web page.

After downloading and extracting the content of the binary archive, please run the following command from DiscoSnp home directory:

    chmod +x run_discoSnp++.sh test/*.sh scripts/*.sh

## Quick start

Run DiscoSnp WITHOUT mapping results on a reference genome:

    ./run_discoSnp++.sh -r test/fof.txt -T

Run DiscoSnp WITH mapping results on a reference genome (requires bwa):

    ./run_discoSnp++.sh -r test/fof.txt -T  -G test/reference_genome.fa

Note: if bwa is not in you PATH, then add the option "-B path_to_bwa". For instance:

    ./run_discoSnp++.sh -r test/fof.txt -T  -G test/reference_genome.fa -B /home/me/my_programs/bwa-0.7.12/

Run DiscoSnp WITH mapping results on a reference genome AND using this reference genome for calling variants:

    ./run_discoSnp++.sh -r test/fof.txt -T  -G test/reference_genome.fa -R

## User manual

See 

    doc/discoSnp_user_guide.pdf or doc/discoSnp_user_guide.txt

## Contact

Remarks and questions: https://www.biostars.org/t/discosnp/

Contact: Pierre Peterlongo: pierre.peterlongo@inria.fr
