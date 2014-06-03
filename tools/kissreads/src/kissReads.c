/**
 * Copyright INRIA , contributors Peterlongo
 * pierre.peterlongo@inria.fr
 *
 *
 * This software is a computer program whose purpose is to detect the
 * presence of a sequence in a set of NGS reads, and to compute its average quality and coverage
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 */

/*
 *
 *  Created on: 28 oct. 2010
 *      Author: Pierre Peterlongo
 */


#include<fragment_index.h>
#include<list.h>
#include<commons.h>
#include<extension_algorithm.h>
#include<outputs.h>
#include<limits.h> // MAXINT
#include <unistd.h>  /* sleep(1) */
#include <string.h> // strdup
#include <time.h>
#include <libchash.h>

//#define OMP

#ifdef OMP
#include "omp.h"
#endif

char * getVersion(){
    return  "1.3.2_check - Copyright INRIA - CeCILL License";
    
}
//#define VERBOSE


void print_usage_and_exit(char * name){
	fprintf (stderr, "NAME\nkissReads, version %s\n", getVersion());
#ifdef INPUT_FROM_KISSPLICE
	fprintf (stderr, "\nSYNOPSIS\n%s <toCheck.fasta> <readsC1.fasta> [<readsC2.fasta> [<readsC3.fasta] ...] [-k value] [-c value] [-d value] [-l value] [-o name] [-u name] [-i index_stride] [-m align_file] [-s] [-f] [-h] \n", name);
#else
    fprintf (stderr, "\nSYNOPSIS\n%s <toCheck.fasta> <readsC1.fasta> [<readsC2.fasta> [<readsC3.fasta] ...] [-k value] [-c value] [-d value] [-o name] [-u name] [-n] [-I] [-i index_stride] [-m align_file] [-s] [-f] [-h] \n", name);
#endif// INPUT_FROM_KISSPLICE
	fprintf (stderr, "\nDESCRIPTION\n");
	fprintf (stderr, "Checks for each sequence contained into the toCheck.fasta if\n");
	fprintf (stderr, "it is read coherent (each position is covered by at least \"min_coverage\" read(s)) with reads from readsA.fasta or readsB.fasta\n");
	fprintf (stderr, "A sequence s from toCheck is treated as follow:\n");
	fprintf (stderr, "  if (s coherent with at least one read set): output the sequence as follows\n");
	fprintf (stderr, "  \t >original fasta comment|C1:min<avg-corr_avg<max|C2:min<avg-cor_avg<max|C3...:\n");
	fprintf (stderr, "  \t >s\n");
	fprintf (stderr, "  With A:min<avg-cor_avg<max standing for : value of the position having minimal coverage in readsA.fasta < average coverage in readsA.fasta - R-squarred corrected average in readsA.fa < value of the position having maximal coverage in readsA.fasta\n");
	fprintf (stderr, "  The coverage is the number of reads that perfectly mapped a position\n");
	fprintf (stderr, "  Any other situation (s not coherent with any): couple non read coherent, not outputed \n");
    
	fprintf (stderr, "\nOPTIONS\n");
	fprintf (stderr, "\t -k size_seed: will use seeds of length size_seed. Default: 25.\n");
	fprintf (stderr, "\t -c min_coverage: a sequence is covered by at least min_coverage coherent reads. Default: 2\n");
    fprintf (stderr, "\t -d max_substitutions: Maximal number of substitutions authorized between a read and a fragment. Note that no substitution is allowed on the central position while anaylizing the kissnp output. Default: 1.\n");
#ifdef INPUT_FROM_KISSPLICE
    fprintf (stderr, "\t -l min_overlap: Kissplice min_voerlap. Default: 3\n");
    fprintf (stderr, "\t -j counting option. 0: counts are summed and represent the whole path, 1: counts are only in the junctions, 2: all the counts are shown separately   . Default: 0\n");
#endif
	fprintf (stderr, "\t -o file_name: write read-coherent outputs. Default: standard output \n");
	fprintf (stderr, "\t -u file_name: write unread-coherent outputs. Default: standard output \n");
#ifndef INPUT_FROM_KISSPLICE
	fprintf (stderr, "\t -n the input file (toCheck.fasta) is a kissnp output (incompatible with -I option) \n");
    fprintf (stderr, "\t\t in this case: 1/ only the upper characters are considered (no mapping done on the extensions) and 2/ the central position (where the SNP occurs) is strictly mapped, no subsitution is authorized on this position.\n");
	fprintf (stderr, "\t -I the input file (toCheck.fasta) is an Intl output (incompatible with -n option) \n");
#endif // not INPUT_FROM_KISSPLICE
    fprintf (stderr, "\t -i index_stride (int value). This is a heuristic for limiting the memory footprint. Instead of indexing each kmer of the sequences contained into the toCheck.fasta, kissreads indexes kmers occurring each \"index_stride\" position. Default = 1 (no heuristic)\n");
	//fprintf (stderr, "\t -c Write the coverage of each position of the sequences. The fasta file will have 3 lines (instead of one) starting with \'>\' before the sequences themselves\n");
    //	fprintf (stderr, "\t -q kmer_size for kissplice: computes qualities for variants obtained with the given kmer_size; use this if you need the quality of a SNP position, or the average quality of an indel or a splicing event\n");
    //	fprintf (stderr, "\t -p number of sets in which to partition the set of variants \n");
    fprintf (stderr, "\t -t max number of threads (also limited by number of input files)\n");
	fprintf (stderr, "\t -m align_file, write a file of reads mapped to sequences in file align_file\n");
	fprintf (stderr, "\t -s silent mode\n");
	fprintf (stderr, "\t --version get the kissReads version and exit\n");
    fprintf (stderr, "\t -f outputs coherent events in a standard fasta file format\n");
	fprintf (stderr, "\t -h prints this message and exit\n");
	exit(0);
}


int main(int argc, char **argv) {
#ifdef LIVIUS_TESTS
    printf("Compiled with LIVIUS_TESTS: output only motifs where au-vb is specific to one datasets and av'-u'b is specific to the other\n");
#endif
#ifdef CLASSICAL_SPANNING
    printf("Compiled with CLASSICAL_SPANNING: a position is considered as covered as soon as a read maps this position (else a position is considered as covered if the kmer starting at this position is fully covered by a read)\n");
#endif
#ifdef INPUT_FROM_KISSPLICE
    printf("Compiled with INPUT_FROM_KISSPLICE: dealing with kissplice output and to count separately junctions and central portions\n");
#endif
    
    
    
    char map_snps=0; // input is a kissnp output
    char map_invs=0; // input is an intl (read2sv) output.
    
    char no_subsutitution_on_central_position=0; // By default: authorize a subsitution on the central position
    char input_only_upper=0; // By default: all characters (upper or lower) are read
    int number_paths_per_event=1; // By default (generic usage) each event is composed by a unique sequence.
#ifdef INPUT_FROM_KISSPLICE
    min_overlap=3; // default value.
    countingOption = 0 ; // default value
    no_subsutitution_on_central_position=0; // just to be sure, as normaly this is the default behaviour
    input_only_upper=1; // don't considere the lower script characters.
    number_paths_per_event=2; // each splicing event is composed by two sequences in the fasta file.
#endif
    
    int index_stride =1;
	setvbuf(stdout,(char*)NULL,_IONBF,0); // disable buffering for printf()
	time_t before_all, after_all;
	before_all = time(NULL);
	size_seeds=25;
    
	int i, k;
    //	kmer_size=0;
    //	nb_event_sets=1;
	int nb_events_per_set;
    int max_substitutions=1;
    int max_threads = 1;
#ifdef OMP
     max_threads =  omp_get_num_procs();
#endif
    //	size_before_reads_starting=16384;
	int min_coverage=2; // minimal number of reads per positions extending the initial fragment
    
	if(argc<2){
		print_usage_and_exit(argv[0]);
	}
    
	char * toCheck_file =strdup(argv[1]);
	FILE * coherent_out = stdout;
	FILE * uncoherent_out = stdout;
	FILE * sam_out = NULL;
	silent=0; verbose=0; whole_coverage=0; standard_fasta=0;
    
    nbits_nbseeds = 8*sizeof(uint64_t)- NBITS_OFFSET_SEED ;
    mask_nbseed  = ( 1ULL <<  (uint64_t) nbits_nbseeds  ) -1 ;
    mask_offset_seed = (1ULL << (NBITS_OFFSET_SEED)) -1 ;
    // printf("mask nb seed %llx  %llx \n",mask_nbseed,mask_offset_seed);
    
	// GET ALL THE READ FILES
	// find the number of read sets
	number_of_read_sets=0;
	while(number_of_read_sets+2<argc && argv[number_of_read_sets+2][0]!='-') number_of_read_sets++;
	char ** reads_file_names = (char **) malloc(sizeof(char *)*number_of_read_sets);
    
	// copy the read set file names
	number_of_read_sets=0;
	while(number_of_read_sets+2<argc && argv[number_of_read_sets+2][0]!='-'){
		reads_file_names[number_of_read_sets]=strdup(argv[number_of_read_sets+2]);
		number_of_read_sets++;
	}
    printf("This is kissreads, version %s\n", getVersion());
    printf("The command line was:\n");
    for(i=0;i<argc;i++) printf("%s ",argv[i]);
    printf("\n");
	while (1)
	{
#ifdef INPUT_FROM_KISSPLICE
        int temoin = getopt (argc-number_of_read_sets-1, &argv[number_of_read_sets+1], "c:d:k:o:u:q:m:i:vfs-:j:l:");
#else
        int temoin = getopt (argc-number_of_read_sets-1, &argv[number_of_read_sets+1], "c:d:k:o:u:q:m:i:vfsnI-:");
#endif //INPUT_FROM_KISSPLICE
		if (temoin == -1){
			break;
		}
		switch (temoin)
		{
            case 'c':
                min_coverage=atoi(optarg);
                if(min_coverage<1) min_coverage=1;
                printf(" Will consider as read coherent if coverage is at least %d\n", min_coverage);
                break;
                
#ifdef INPUT_FROM_KISSPLICE
            case 'l':
                min_overlap=atoi(optarg);
                printf(" Kissplice count uses min_overlap %d\n", min_overlap);
                break;
            case 'j':
              countingOption  = atoi(optarg);
              printf(" Kissplice countingOption %d\n", countingOption);
              break;
#endif //INPUT_FROM_KISSPLICE
            case 'o':
                coherent_out = fopen(optarg, "w");
                if(coherent_out == NULL){
                    fprintf(stderr,"cannot open %s for writing results, exit\n", optarg);
                    exit(1);
                }
                printf(" Will output read-coherent results in file: %s\n", optarg);
                break;
            case 'u':
		        uncoherent_out = fopen(optarg, "w");
                if(uncoherent_out == NULL){
                    fprintf(stderr,"cannot open %s for writing results, exit\n", optarg);
                    exit(1);
                }
                printf(" Will output unread-coherent results in file: %s\n", optarg);
                break;
            case 'd':
                max_substitutions=atoi(optarg);
                printf(" Will authorize up to %d substitutions while mapping a read on a fragment.\n", max_substitutions);
                break;
            case 't':
                max_threads=atoi(optarg);
                printf(" Will authorize up to %d threads.\n", max_threads);
                break;
            case 'i':
                index_stride=atoi(optarg);
                if(index_stride<0) index_stride=1;
                printf(" Will index a kmer each %d position\n", index_stride);
                break;
                
#ifndef INPUT_FROM_KISSPLICE
            case 'n':
                map_snps=1;
                no_subsutitution_on_central_position=1;
                input_only_upper=1;
                number_paths_per_event=2;
                
                printf(" Will treat the input file as a kissnp output\n");
                break;
            case 'I':
                map_invs=1;
                no_subsutitution_on_central_position=0;
                input_only_upper=0; // read all characters for now (still no extension pluged to intl).
                number_paths_per_event=4;
                
                printf(" Will treat the input file as an Intl output\n");
                break;
#endif
                
            case 'k':
                size_seeds=atoi(optarg);
                
                if(size_seeds<5){
                    fprintf(stderr,"%d too small to be used as seed length. kisSnpCheckReads continues with seeds of length 5\n",size_seeds);
                    size_seeds=10;
                }
                else
                    printf(" Will use seeds of length %d\n", size_seeds);
                break;
            case 'h':
                print_usage_and_exit(argv[0]);
            case 'v':
                printf(" VERBOSE mode\n");
                verbose=1;
                break;
            case 'f':
                printf("Standard fasta output mode\n");
                standard_fasta=1;
                break;
            case 's':
                printf(" SILENT mode\n");
                verbose=0;
                silent=1;
                break;
            case 'm':
                sam_out = fopen(optarg, "w");
                if(sam_out == NULL){
                    fprintf(stderr,"cannot open %s for writing sam like results, exit\n", optarg);
                    exit(1);
                }
                fprintf(sam_out,"SNP_header\tmapped_read\tid_read_file\tposition_where_read_is_mapped\n");
                printf(" Will output reads mapped results in file: %s\n", optarg);
                break;
                //            case 'q':
                //                printf(" Computes qualities for variants obtained with the given kmer_size, as used by kisSplice\n");
                //                kmer_size=atoi(optarg);
                //                quality=1;
                //                break;
                //            case 'p':
                //                nb_event_sets=atoi(optarg);
                //                printf(" Partitions the set of fasta sequences whose read coherence is to be checked in %d sets\n", nb_event_sets);
                //                break;
                
            case '-':
                if(strcmp("version", optarg)==0){
                    printf("kissReads version %s\n", getVersion());
                    exit(0);
                }
                printf(" what next ? %s\n",optarg);
                break;
            default:
                //			printf ("Unknown option %c\n", temoin);
                print_usage_and_exit(argv[0]);
		}
	}
    
	if ( argc  - optind <2)
	{
		print_usage_and_exit(argv[0]);
	}
    
#ifndef INPUT_FROM_KISSPLICE
    if(map_snps && map_invs){
        fprintf(stderr, "cannot use both options -n and -I\n");
        exit(1);
    }
    if(!map_snps && !map_invs){
        fprintf(stderr, "sorry, kissreads is not ready yet to deal with generic input (neither kissnp, kissreads, nor intl)\n");
        fprintf(stderr, " -- compile with option INPUT_FROM_KISSPLICE (make MYCFLAGS=-DINPUT_FROM_KISSPLICE)\n");
        fprintf(stderr, "   or\n");
        fprintf(stderr, " -- use either option -n or -I\n bye\n");
        exit(1);
    }
#endif 
    
    
    if(!silent){
#ifdef KMER_SPANNING
        printf("Each kmer(=%d) spanning each position should be covered by at least %d reads.\n", size_seeds, min_coverage);
#else
        printf("Each position should be covered by at least %d reads.\n", min_coverage);
#endif
    }
	init_static_variables(size_seeds);
    //	float sum_average_read_length=0;
	
	gzFile * reads_files = (gzFile*) malloc(sizeof(gzFile)*number_of_read_sets);
	p_fragment_info * results_against_set;
    
	char *found;
	// test if all formats are in fastq
	quality=1;
	for (i=0;i<number_of_read_sets;i++){
        quality &= strstr(reads_file_names[i],"fastq") || strstr(reads_file_names[i],"fq") || strstr(reads_file_names[i],"txt");
	}
	if(quality) quality=1;
    
    
	for (i=0;i<number_of_read_sets;i++){
        if (quality){
            found = strstr(reads_file_names[i],"fastq") || strstr(reads_file_names[i],"fq") || strstr(reads_file_names[i],"txt");
            if (!found ){		fprintf(stderr,"\nwith q option, fastq reads files are needed, wrong %s, exit\n",reads_file_names[i]);		exit(1);	}
        }
        else
	    {
            found = strstr(reads_file_names[i],"fastq") || strstr(reads_file_names[i],"fq") || strstr(reads_file_names[i],"txt") || strstr(reads_file_names[i],"fa") || strstr(reads_file_names[i],"fasta") || strstr(reads_file_names[i],"fna");
            if (!found ){		fprintf(stderr,"\naccepted extensions are .fa[.gz], .fasta[.gz], .fna[.gz] for fasta format and .fq[.gz], .fastq[.gz], .txt[.gz] for fastq format, wrong %s, exit\n",reads_file_names[i]);		exit(1);	}
	    }
	}
	
	file=gzopen(toCheck_file,"r");
	if(file == NULL){
        fprintf(stderr,"cannot open file %s, exit\n",toCheck_file);
        exit(1);
	}
    if(!silent)
        printf("counts the number of fragments we have\n");
	// how many fragments do we have (fragments are in fastq format)?
    char * line =  (char *)malloc(sizeof(char)*1048576); // 1048576
	number_of_starters = number_of_sequences_in_file(file,line);
    free(line);
    //
    //    if(map_snps)
    //        nb_events_per_set = (number_of_starters/2)/nb_event_sets;
    //    else if(map_invs)
    //        nb_events_per_set = (number_of_starters/4)/nb_event_sets;
    //    else
    //        nb_events_per_set = number_of_starters/nb_event_sets; // TODO FINISH GENERIC
    
    nb_events_per_set = (number_of_starters/number_paths_per_event);
    
    
    if(!silent){
        //        if(map_snps) printf("\nThe number of sequences in %s file is %d, therefore the number of mouths is %d\n", toCheck_file, number_of_starters, nb_events_per_set);
        //        else if(map_invs) printf("\nThe number of sequences in %s file is %d, therefore the number of inversions is %d\n", toCheck_file, number_of_starters, nb_events_per_set);
        //        else printf("\nThe number of sequences in %s file is %d (%d events)\n", toCheck_file, number_of_starters, nb_events_per_set);
        printf("\nThe number of sequences in %s file is %d (%d events)\n", toCheck_file, number_of_starters, nb_events_per_set);
    }
    
    
    
    //	for ( j = 1; j<=nb_event_sets; j++ )
    //    {
    //	    if ( j == nb_event_sets ){ // FIXME: ?
    //            if(map_snps) nb_events_per_set = (number_of_starters/2) - nb_events_per_set * (nb_event_sets-1);
    //            else if(map_invs) nb_events_per_set = (number_of_starters/4) - nb_events_per_set * (nb_event_sets-1);
    //            else nb_events_per_set = (number_of_starters) - nb_events_per_set * (nb_event_sets-1);
    //        }
    
    if(!silent)	printf("\nIndexing\n");
    results_against_set = index_starters_from_input_file (size_seeds, nb_events_per_set, number_paths_per_event, input_only_upper, index_stride);
    
     int nbthreads = number_of_read_sets;
     if (nbthreads > max_threads) nbthreads = max_threads;
#ifdef OMP
    silent=1; // avoids melting messages
#endif
    
#ifdef OMP
#pragma omp parallel for if(nbthreads>1 && sam_out==NULL) num_threads(nbthreads) private(i)
#endif
    for (i=0;i<number_of_read_sets;i++){
        reads_files[i]=gzopen(reads_file_names[i],"r");
        if(reads_files[i] == NULL){		fprintf(stderr,"cannot open reads file %s, exit\n",reads_file_names[i]);		exit(1);	}
        
        if(!silent) printf("\nCheck read coherence... vs reads from %s (set %d)\n", reads_file_names[i], i);
        
        
        
        read_coherence(reads_files[i], reads_file_names[i], size_seeds,  min_coverage, results_against_set,  number_of_starters, i, quality, nb_events_per_set,  number_paths_per_event, sam_out, max_substitutions, no_subsutitution_on_central_position);
        
        
        gzclose(reads_files[i]);
    }
    
    
    if(number_paths_per_event==2)
        print_results_2_paths_per_event(coherent_out, uncoherent_out, results_against_set, number_of_read_sets, nb_events_per_set, quality);
    else if(map_invs)
        print_results_invs(coherent_out, uncoherent_out, results_against_set, number_of_read_sets, nb_events_per_set, quality);
    //    else   print_generic_results(coherent_out, uncoherent_out, results_against_set, number_of_read_sets, nb_events_per_set, quality); // TODO
    
    free(seed_table);
    FreeHashTable((struct HashTable*)seeds_count);
    for (i=0;i<nb_events_per_set;i++)
    {
        free(results_against_set[i]->read_coherent);
        free(results_against_set[i]->number_mapped_reads);
        for (k=0;k<number_of_read_sets;k++)
        {
            free(results_against_set[i]->read_coherent_positions[k]);
            free(results_against_set[i]->sum_quality_per_position[k]);
        }
        free(results_against_set[i]->read_coherent_positions);
        free(results_against_set[i]->sum_quality_per_position);
        free(results_against_set[i]);
    }
    free(results_against_set);
    //} // nb event set
    
    
	//average_size_reads=sum_average_read_length/number_of_read_sets;
    
	gzclose(file);
	fclose(coherent_out);
	fclose(uncoherent_out);
	if(sam_out) fclose(sam_out);
	if(!silent) {
		after_all = time(NULL);
		printf("Total time: %.0lf secs\n",  difftime(after_all, before_all));
	}
	return 0;
}
