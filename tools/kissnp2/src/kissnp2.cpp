#ifdef ED
//Copyright inria / irisa (2013)
//
//
//raluca.uricaru@gmail.com
//pierre.peterlongo@inria.fr
//
//This software is a computer program whose purpose is to call SNPs from NGS reads.
//
//This software is governed by the CeCILL license under French law and
//abiding by the rules of distribution of free software.  You can  use,
//modify and/ or redistribute the software under the terms of the CeCILL
//license as circulated by CEA, CNRS and INRIA at the following URL
//"http://www.cecill.info".
//
//As a counterpart to the access to the source code and  rights to copy,
//modify and redistribute granted by the license, users are provided only
//with a limited warranty  and the software's author,  the holder of the
//economic rights,  and the successive licensors  have only  limited
//liability.
//
//In this respect, the user's attention is drawn to the risks associated
//with loading,  using,  modifying and/or developing or reproducing the
//software by the user in light of its specific status of free software,
//that may mean  that it is complicated to manipulate,  and  that  also
//therefore means  that it is reserved for developers  and  experienced
//professionals having in-depth computer knowledge. Users are therefore
//encouraged to load and test the software's suitability as regards their
//requirements in conditions enabling the security of their systems and/or
//data to be ensured and,  more generally, to use and operate it in the
//same conditions as regards security.
//
//The fact that you are presently reading this means that you have had
//knowledge of the CeCILL license and that you accept its terms.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h> // for mkdir
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>
#include <sstream> // for ostringstream

#define NNKS 4 // default minimal abundance for solidity

int max_memory; // the most memory one should alloc at any time, in MB
int order=0; // in kissnp, don't change it, it should be 0
#define MIN_CONTIG_SIZE 108

#include "../minia/Bank.h"
#include "../minia/Hash16.h"
#include "../minia/Set.h"
#include "../minia/Pool.h"
#include "../minia/Bloom.h"
#include "../minia/Debloom.h"
#include "../minia/Utils.h"
#include "../minia/rvalues.h"
#include "../minia/SortingCount.h"
#include "../minia/Terminator.h"
#include "../minia/Kmer.h"
#include "SNP.h"
#include "Kmer_for_kissnp2.h"
#include "commons.h"

int64_t genome_size = 3000000000;
Bloom * bloo1;
int threshold;
BinaryBank * SolidKmers;
BranchingTerminator * terminator;
extern int size_seeds;

char * getVersion(){
    return (char *)"DiscoSnp - kissnp submodule 1.2.0 - Copyright INRIA - CeCILL License";
}

void print_usage_and_exit(char * name){
    fprintf (stderr, "NAME\nkissnp, version %s\n", getVersion());
    fprintf (stderr, "\nSYNOPSIS\n%s <readsC1.fasta> [<readsC2.fasta> [<readsC3.fasta] ...] -o name [-t] [-e length] [-l] [-b] [-k value] [-c value] [-g value] [-h]\n", name);
    fprintf (stderr, "\t or:\n");
    fprintf (stderr, "%s input_names_file.txt -o name [-t] [-e length] [-l] [-b] [-k value] [-c value] [-g value] [-h]\n", name);
    fprintf (stderr, "\t with \"input_names_file.txt\" being a file containing on each line the name of read files\n");
    fprintf (stderr, "\nDESCRIPTION\n");
    fprintf (stderr, "\t \"kissnp2\", detects SNPs from read set(s). It should usually be followed by the \"kissreads\" for recovering coverage and quality information from reads.\n");
    fprintf (stderr, "\nMANDATORY\n");
    
    fprintf (stderr, "\t At least one read set. \n");
    fprintf (stderr, "\t -o file_name_prefix: where to write outputs and debruijn graph structure files. \n");
    
    fprintf (stderr, "\nOPTIONS\n");
    fprintf (stderr, "\t -t extend found and stop at first polymorphism (strict extension=unitigs) SNPs. Uncompatible with -T\n");
    fprintf (stderr, "\t -T extend found and stop at large polymorphism (extension=contigs) SNPs. Uncompatible with -t\n");
    fprintf (stderr, "\t -e length: extend found SNPs (option -t) and conserve only those whose min(left and right extension) is bigger or equal to \"length\"\n");
    fprintf (stderr, "\t -l conserve low complexity SNPs. Default: false (filter out low complexity results)\n");
    fprintf (stderr, "\t -b INT:\n");
    fprintf (stderr, "\t \t 0: forbid SNPs for wich any of the two paths is branching (high precision, low recall)\n");
    fprintf (stderr, "\t \t 1: forbid SNPs for wich the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position (default value)\n");
    fprintf (stderr, "\t \t 2: No limitation on branching (low precision, high recall)\n");
    fprintf (stderr, "\t -k size_seed: will use seeds of length size_seed. Default: 27.\n");
    fprintf (stderr, "\t -c min_coverage: a sequence is covered by at least min_coverage coherent reads. Default: 2\n");
    fprintf (stderr, "\t -C max_coverage: a sequence is covered by at most max_coverage coherent reads. Default: infiny (=%d on your computer :) )\n", INT_MAX);
    fprintf (stderr, "\t -g estimated_genome_size: estimation of the size of the genome whose reads come from. \n \t    It is in bp, does not need to be accurate, only controls memory usage. Default: 3 billion\n");
    fprintf (stderr, "\t -h prints this message and exit\n");
    exit(0);
}



//31 /local/ruricaru/workspace_CodeBlocks/kissnp/head_bubbles_kissnp_7aphids
//3 /local/ruricaru/workspace_CodeBlocks/kissnp/xaa /local/ruricaru/workspace_CodeBlocks/kissnp/xab /local/ruricaru/workspace_CodeBlocks/kissnp/xac  25 3 3000000 t3
int main(int argc, char *argv[])
{
    printf("This is kissnp2, version 2.1.1.5\n", getVersion());
    printf("The command line was:\n");
    for(int i=0;i<argc;i++) printf("%s ",argv[i]);
    printf("\n");
    sizeKmer=27; // let's make it even for now, because i havnt thought of how to handle palindromes (dont want to stop on them)
    int min_size_extension=-1; // TODO FINR: ajouter and constructeur SNP et utiliser lors du output.
    // solidity
    nks =NNKS;
    int low=0;
    /* authorised_branching =
    *  - 0: branching forbiden in any path
    *  - 1: same branching on both path forbiden (i.e. 2 disctinct nucelotides may be used in both paths for extension)
    *  - 2: no restriction on branching
    */
    int authorised_branching=1;
    bool extend_snps=false;
    bool print_extensions = true;

    ////////////////////////////////// WITH GET OPTS
    
    ////////////////////////////////// GETTING THE OPTIONS /////////////////////////////////////
    
    // GET ALL THE READ FILES
    // find the number of read sets
    int number_of_read_sets=0;
    while(number_of_read_sets+1<argc && argv[number_of_read_sets+1][0]!='-') number_of_read_sets++;
    if(number_of_read_sets<1) print_usage_and_exit(argv[0]);
    printf("%d read set(s):\n",number_of_read_sets); for(int i=0;i<number_of_read_sets;i++) printf("\t %s\n", argv[1+i]);
    Bank *Reads = new Bank(argv+1, number_of_read_sets);
    char * SNP_file_name = new char [4096];
    bool outputfilename=false;
    bool strict_extension=true;
    
    while (1)
    {
        int temoin = getopt (argc-number_of_read_sets, &argv[number_of_read_sets], "Tte:lb:g:c:C:o:k:hx");
        if (temoin == -1){
            break;
        }
        switch (temoin)
        {
                
            case 'x':
                extend_snps=true;
                strict_extension=true;
                print_extensions=false;
                printf(" kisnnp will extend snps but not print them\n");
                break;
            case 't':
                extend_snps=true;
                strict_extension=true;
                printf(" kisnnp will extend snps - unitig mode\n");
                break;
            case 'T':
                extend_snps=true;
                strict_extension=false;
                printf(" kisnnp will extend snps - contig mode\n");
                break;
            case 'e':
                extend_snps=true;
                min_size_extension=atoi(optarg);
                printf(" kisnnp will extend snps and conserve only those whose min extension (left or right) is bigger than %d\n", min_size_extension);
                break;
            case 'l':
                low=1;
                printf(" kisnnp will output low complexity snps\n");
                break;
            case 'b':
                authorised_branching=atoi(optarg);
                printf("authorised_branching value=%d\n", authorised_branching);
                if(authorised_branching<0 || authorised_branching>2){
                    fprintf(stderr,"Value %d forbiden for the authorised_branching value, exit\n", authorised_branching);
                    exit(1);
                }
                break;
            case 'g':
                genome_size=atoll(optarg);
                printf(" Will use genome_size=%llu\n", genome_size); //TODO format
                break;
            case 'c':
                nks = atoi(optarg);
                if (nks<1) nks=1; // min abundance can't be 0
                printf(" Will create SNPs with k-mers occurring at least %d time(s)\n", nks);
                break;
            case 'C':
                max_couv = atoi(optarg);
                if (nks<1) nks=1; // min abundance can't be 0
                printf(" Will create SNPs with k-mers occurring at most %d time(s)\n", max_couv);
                break;

            case 'o':
                strcpy (SNP_file_name, optarg);
                outputfilename=true;
                //SNP_file_name=(char *)string(optarg).c_str();
                printf(" Will output results in file: %s.fa\n", optarg);
                break;                
            case 'k':
                sizeKmer=atoi(optarg);
                printf("required size seeds = %d\n", sizeKmer);
                break;
            case 'h':
                print_usage_and_exit(argv[0]);
                break;
            case '-':
                if(strcmp("version", optarg)==0){
                    printf("kissReads version %s\n", getVersion());
                    exit(0);
                }
                printf(" what next ? %s\n",optarg);
                break;
            default:
                printf ("Unknown option %c\n", temoin);
                print_usage_and_exit(argv[0]);
        }
    }
//    printf("%d %d", argc, optind);
    if ( argc  - optind <1 || !outputfilename)
    {
        print_usage_and_exit(argv[0]);
    }
    
    
    
    ///////////////////////////////////////////////////// END GET OPTS ////// 
    
    
    
    

    if (sizeKmer%2==0)
    {
        sizeKmer-=1;
        printf("Need odd kmer size to avoid palindromes. I've set kmer size to %d.\n",sizeKmer);
    }
    if (sizeKmer<11){
        printf("I can't use k so small(=%d), I fix it to 11\n",sizeKmer);
        sizeKmer=11;
    }
    if (sizeKmer>(sizeof(kmer_type)*4))
    {
        printf("Max kmer size on this compiled version is %lu\n",sizeof(kmer_type)*2);
        exit(1);
    }
    size_seeds=sizeKmer; // TODO: conserve only one of this two redundant variables
    threshold = (sizeKmer/2-2)*(sizeKmer/2-3);
    init_static_variables();// kmer size
    kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
    
    
    sprintf(prefix,"%s_k_%d_c_%d", SNP_file_name, sizeKmer, nks);
    
    sprintf(SNP_file_name,"%s_k_%d_c_%d", SNP_file_name, sizeKmer, nks);
    if(min_size_extension>-1)
        sprintf(SNP_file_name,"%s_e_%d", SNP_file_name, min_size_extension);
    //sprintf(SNP_file_name_too_many, "%s_too_many", SNP_file_name);
    
    if(max_couv < INT_MAX) printf("will eliminate kmers with coverage > %i \n",max_couv);
       ///////////////////////////  CREATE OR LOAD THE GRAPH
    printf("Indexing reads, using minia approach, generating file prefixed with \"%s\"\n", prefix);
    
    double lg2 = log(2);
   // NBITS_PER_KMER = log(16*sizeKmer*(lg2*lg2))/(lg2*lg2); // needed to process argv[5]
    NBITS_PER_KMER = rvalues[sizeKmer][1];
    
    
    int estimated_BL1 = max( (int)ceilf(log2f(genome_size * NBITS_PER_KMER )), 1);
    
    uint64_t estimated_nb_FP =  (uint64_t)(genome_size * 4 * powf(0.6,11)); // just indicative
    
    max_memory = max( (1LL << estimated_BL1)/8LL /1024LL/1024LL, 1LL );
    printf("estimated values: nbits Bloom %i, nb FP %lld, max memory %i MB\n",estimated_BL1,estimated_nb_FP,max_memory);
    
    kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
    NBITS_PER_KMER = log(16*sizeKmer*(lg2*lg2))/(lg2*lg2); // needed to process argv[5]
    max_memory = max( (1LL << estimated_BL1)/8LL /1024LL/1024LL, 1LL );
    //  printf("estimated values: nbits Bloom %i, nb FP %lld, max memory %i MB\n",estimated_BL1,estimated_nb_FP,max_memory);
    
    
    
    
    
    // shortcuts to go directly to assembly using serialized bloom and serialized hash
    int START_FROM_SOLID_KMERS=0; // if = 0, construct the fasta file of solid kmers, if = 1, start directly from that file
    int LOAD_FALSE_POSITIVE_KMERS=0; // if = 0, construct the fasta file of false positive kmers (debloom), if = 1, load that file into the hashtable
    int NO_FALSE_POSITIVES_AT_ALL=0; // if = 0, normal behavior, if = 1, don't load false positives (will be a probabilistic de bruijn graph)
    
    printf("SNPFILE = %s\n", SNP_file_name);
    
    if( access( (string(prefix)+string(".debloom")).c_str(), F_OK ) != -1 &&
       access( (string(prefix)+string(".debloom2")).c_str(), F_OK ) != -1 &&
       access( (string(prefix)+string(".false_positive_kmers")).c_str(), F_OK ) != -1 &&
       access( (string(prefix)+string(".solid_kmers_binary")).c_str(), F_OK ) != -1){
        printf("THE INDEX %s.... were already created we use them\n", prefix);
        START_FROM_SOLID_KMERS=1;
        LOAD_FALSE_POSITIVE_KMERS=1;
    }
    else
        printf("CREATING THE index %s.... \n", prefix);
    
    
    printf("SNPFILE = %s\n", SNP_file_name);
    
    //  fprintf (stderr,"taille cell %lu \n", sizeof(cell<kmer_type>));
    

    STARTWALL(0);
    
    // counter kmers, write solid kmers to disk, insert them into bloo1
    if (!START_FROM_SOLID_KMERS)
    {
        
        int max_disk_space = 0; // let dsk decide
        int verbose = 0;
        bool write_count = false;
        
        sorting_count(Reads,prefix, max_memory, max_disk_space, write_count, verbose);
    }
    
    
    // debloom, write false positives to disk, insert them into false_positives
    if (! LOAD_FALSE_POSITIVE_KMERS)
    {
        debloom(order, max_memory);
    }
    
    bloo1 = bloom_create_bloo1((BloomCpt *)NULL);
    
    // load false positives from disk into false_positives
    false_positives = load_false_positives_cascading4();
    
    
    // load branching kmers
    BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer_type),0);
    int LOAD_BRANCHING_KMERS=0;
    if (LOAD_BRANCHING_KMERS)
    {
        BinaryBank *BranchingKmers = new BinaryBank(return_file_name(branching_kmers_file),sizeof(kmer_type),false);
        terminator = new BranchingTerminator(BranchingKmers,SolidKmers, bloo1,false_positives);
        BranchingKmers->close();
    }
    else{
        BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer_type),0);
        terminator = new BranchingTerminator(SolidKmers,genome_size, bloo1,false_positives);
        
    }

    
    ///////////////////////////  FIND BUBBLES
    

    //fprintf (stderr,"\t Finding bubbles: call \"new Bubble\" function\n");
    fprintf (stderr,"\t Finding bubbles, printing results in %s", SNP_file_name);
    Bubble *bubble = new Bubble(SolidKmers,bloo1,false_positives, extend_snps, min_size_extension,print_extensions, strict_extension);
    //fprintf (stderr,"\t Finding bubbles: call \"find_bubble\" function, printing results in %s", SNP_file_name);
    bubble->find_bubbles((string(SNP_file_name)+string(".fa")).c_str(), low, authorised_branching);
    //fprintf (stderr,"\t Finding bubbles: done... call of STOPWALL before to exit, tcho\n");

    STOPWALL(0,"Total");

    delete Reads;

    return 0;
}

#endif
