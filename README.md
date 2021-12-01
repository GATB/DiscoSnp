# DiscoSnp++ and DiscoSnpRad


| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/DiscoSnp-gitlab/job/tool-discosnp-build-debian7-64bits-gcc-4.7-gitlab/badge/icon)](https://ci.inria.fr/gatb-core/view/DiscoSnp-gitlab/job/tool-discosnp-build-debian7-64bits-gcc-4.7-gitlab/) | [![Build Status](https://ci.inria.fr/gatb-core/view/DiscoSnp-gitlab/job/tool-discosnp-build-macos-10.9.5-gcc-4.2.1-gitlab/badge/icon)](https://ci.inria.fr/gatb-core/view/DiscoSnp-gitlab/job/tool-discosnp-build-macos-10.9.5-gcc-4.2.1-gitlab/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is DiscoSnp++?

DiscoSnp is designed for discovering all kinds of SNPs (not only isolated ones),  as well as insertions and deletions, from raw set(s) of reads. The number of input read sets is not constrained, it can be one, two, or more. No reference genome is needed.

## Publications

Uricaru R., Rizk G., Lacroix V., Quillery E., Plantard O., Chikhi R., Lemaitre C., Peterlongo P. (2014). [Reference-free detection of isolated SNPs](http://nar.oxfordjournals.org/content/43/2/e11). Nucleic Acids Research 43(2):e11.

Peterlongo, P., Riou, C., Drezen, E., Lemaitre, C. (2017). [DiscoSnp ++ : de novo detection of small variants from raw unassembled read set(s).](http://doi.org/https://doi.org/10.1101/209965) BioRxiv.

Gauthier, J., Mouden, C.,  Suchan, T., Alvarez, N., Arrigo, N., Riou, C., Lemaitre, C., Peterlongo, P. (2017). [DiscoSnp-RAD: de novo detection of small variants for population genomics](https://peerj.com/articles/9291/). *PeerJ* 8 (2020): e9291.

## DiscoSnp++ or DiscoSnpRad
We propose a DiscoSnp++ adaptation for RAD-Seq data. A script, called `run_discoSnpRad.sh`, is adapted to this kind of data. See below for more details.

# Getting the latest source code

## Requirements

CMake 2.6+; see [http://www.cmake.org/cmake/resources/software.html](http://www.cmake.org/cmake/resources/software.html)

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions

### Install from source

```bash
# get a local copy of DiscoSnp source code
git clone --recursive https://github.com/GATB/DiscoSnp.git

# compile the code an run a simple test on your computer
cd DiscoSnp
sh INSTALL
```

### Install using Conda

DiscoSnp++ and DiscoSnpRAD are also distributed as a [Bioconda package](https://anaconda.org/bioconda/discosnp):

```
conda install -c bioconda discosnp
```

The scripts `run_discoSnp++.sh` and `run_discoSnpRad.sh` are then executable.

# Getting a binary stable release

Binary release for Linux and Mac OSX are provided within the "Releases" tab on Github/DiscoSnp web page.

After downloading and extracting the content of the binary archive, please run the following command from DiscoSnp home directory:

    chmod +x run_discoSnp++.sh test/*.sh scripts/*.sh

# Quick start

Run DiscoSnp WITHOUT mapping results on a reference genome:

    ./run_discoSnp++.sh -r test/fof.txt -T

Run DiscoSnp WITH mapping results on a reference genome (requires bwa):

    ./run_discoSnp++.sh -r test/fof.txt -T  -G test/reference_genome.fa

Note: if bwa is not in you PATH, then add the option "-B path_to_bwa". For instance:

    ./run_discoSnp++.sh -r test/fof.txt -T  -G test/reference_genome.fa -B /home/me/my_programs/bwa-0.7.12/

Run DiscoSnp WITH mapping results on a reference genome AND using this reference genome for calling variants:

    ./run_discoSnp++.sh -r test/fof.txt -T  -G test/reference_genome.fa -R

# User manual

See doc/discoSnp_user_guide.pdf or doc/discoSnp_user_guide.txt

# DiscoSnpRad
When dealing with RAD-Seq data,  `run_discoSnpRad.sh` script should be used. It uses options specific to RAD-Seq data: branching strategy, kind of extensions, abundance threshold, and kind of bubbles to be found. Moreover, it clusters variants per locus, the cluster information being reported in the final provided VCF file. 

See the [discoSnpRAD README file](https://github.com/GATB/DiscoSnp/tree/master/discoSnpRAD).

# Contact

Remarks and questions: [https://www.biostars.org/t/discosnp/](https://www.biostars.org/t/discosnp/)

Contact: Pierre Peterlongo: [pierre.peterlongo@inria.fr](mailto:pierre.peterlongo@inria.fr)
