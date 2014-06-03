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
 * extension_algorithm.c
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */

#include <extension_algorithm.h>
#include<coherence_algorithm.h>
#include <fragment_index.h>
#include <stdio.h>
#include <commons.h>
#include<string.h>
#include<stdlib.h>
#include<hash.h>
#include<couple.h>
#include<list.h>
#include<limits.h>
#include<assert.h>


//#define DEBUG_MAPPING
//#define DEBUG_QUALITY
#define min(a, b) ((a) < (b) ? (a) : (b))
int number_of_reads;


// checks if, for a fragment, a tuple read/position mapped already (return 1) or none (return 0)
int exists_tuple_read_position(const hash_t map, const char * starter, const int pos){
	list * read_coherent_positions;                 //for one couple (read/fragment) : list of coherent mapped positions (usually, one position)
	if(hash_entry_by_key(map, starter, (void **) (&read_coherent_positions)) !=0){
		cell *c = ((list *)read_coherent_positions)->first;
		while(c != NULL){
			if((*(int *) (c->val))== pos){
				return 1;}
			c = c->prox;
		}
	}
	return 0;
}

void print_mapping(const int pos_on_fragment, const char * fragment, const char *read){
	int i;
	printf("pos: %d\n", pos_on_fragment);
	for(i=0;i<60;i++) printf(" ");
	printf("%s\n", fragment);
    
	for(i=0;i<60+pos_on_fragment;i++) printf(" ");
	printf("%s\n", read);
    
}



// TODO: stop when all positions are bigger than a
void feed_coherent_positions(p_fragment_info the_starter, const int start, const int length_read, const int qual, char *quality, int start_on_read, int read_file){
    int i, start_on_starter, stop_on_starter;
	if(start<0) start_on_starter=0; else start_on_starter=start;
	if(start+length_read<strlen(the_starter->w)) stop_on_starter=start+length_read; else stop_on_starter=strlen(the_starter->w);
	
    
    
    
	// FEED THE STARTING POSITIONS:
	//  ---------------[]------------**************************
	//    before fragment            |  starter
	//    size_before_reads_starting |    strlen(starter)
	//
	// starter[i] corresponds to the_starter->reads_starting[i+size_before_reads_starting]
	// TOO HEAVY,
	//the_starter->reads_starting[start+size_before_reads_starting]++;
	// REPLACED BY:
	the_starter->number_mapped_reads[read_file]++;
    
    
#ifdef INPUT_FROM_KISSPLICE
    if (the_starter->upperpath){
        // *************************
      // alice 2/10/13 : new kissplice output, 1 nt of context added
        // <-k--><----?----><-k->
        //    A        S        B
      // overlap  on junctions : need one more

        char overlap_AS=0;
        char overlap_SB=0;
        
        if(start_on_starter<=size_seeds-min_overlap && stop_on_starter>=size_seeds+min_overlap) overlap_AS=1;
        const int pos_SB = strlen(the_starter->w)-(size_seeds);
        if(start_on_starter<=pos_SB-min_overlap && stop_on_starter>=pos_SB+min_overlap) overlap_SB=1;
        if(overlap_AS) the_starter->nb_reads_overlapping_AS[read_file]++;
        if(overlap_SB) the_starter->nb_reads_overlapping_SB[read_file]++;
        if(overlap_AS && overlap_SB) the_starter->nb_reads_overlapping_both_AS_and_SB[read_file]++;
        if(start_on_starter>size_seeds-min_overlap && stop_on_starter<pos_SB+min_overlap) the_starter->nb_reads_fully_in_S[read_file]++;
    }
    else{// lower path
      if(stop_on_starter-start_on_starter>size_seeds+min_overlap)
            the_starter->nb_reads_overlapping_both_AS_and_SB[read_file]++;
//            else
//                printf("%s\n%d\n",  the_starter->comment, the_starter->nb_reads_overlapping_both_AS_and_SB[read_file]);
    }
    
#endif
    
#ifdef CHARQUAL
    
    if ( qual )
    {
        //GR 06/11/2012 : cumulative moving average, allow quals to be stored in uchar
        //but average is computed only over 255 first read mappeds, because read_coherent_positions also changed to uchar
        unsigned char  nbreads;
        float new_val;
        for(i=start_on_starter;i<stop_on_starter;i++)
        {
            nbreads=the_starter->read_coherent_positions[read_file][i];
			if (nbreads == 0) //first read mapped
				the_starter->sum_quality_per_position[read_file][i] = quality[start_on_read + i - start_on_starter];
			else //at least 1 read mapped already
            {
                new_val =  ((float)the_starter->sum_quality_per_position[read_file][i]) *(nbreads-1) +quality[start_on_read + i - start_on_starter] ;
                new_val = new_val / nbreads;
                the_starter->sum_quality_per_position[read_file][i] = (unsigned char) new_val ;
            }
        }
    }
#else
    if ( qual )
        for(i=start_on_starter;i<stop_on_starter;i++)
            if (the_starter->read_coherent_positions[read_file][i] == 0) //first read mapped
                the_starter->sum_quality_per_position[read_file][i] = quality[start_on_read + i - start_on_starter];
            else //at least 1 read mapped already
                the_starter->sum_quality_per_position[read_file][i] += (int) quality[start_on_read + i - start_on_starter];
#endif
    
    
#ifdef KMER_SPANNING
    // the position i is contained into a kmer fully contained into only 1 mapped read, return 1
    // for doing this we stored on each position of the fragment the number of k-mers starting at this position that fully belong to a read that was mapped
	
    //  -------------------------------------------------------------- starter
    //            ************************
    //                       <-----k----->
    //  00000000001111111111110000000000000000000000000000000000000000 the_starter->read_coherent_positions[read_file]
	if(start+length_read-size_seeds<strlen(the_starter->w)) stop_on_starter=start+length_read-size_seeds; else stop_on_starter=strlen(the_starter->w);
	//	printf("KMER SPANNING, STOP = %d\n", stop_on_starter);
	//#else
	//	printf("FULL SPANNING, STOP = %d\n", stop_on_starter);
#endif
    
	for(i=start_on_starter;i<stop_on_starter;i++) Sinc8(the_starter->read_coherent_positions[read_file][i]);
	


}



// For each position in the starter:
// find the kmer supported by the smallest number of reads, e.g.:
//                                i
//  ------------------------------X------------------------------- starter
// a         ************************
// b         **************************
// c                    *******************
// d                       *************************
//                     <-----k-----> : 2 reads (a and b)
//                      <-----k-----> : 3 reads (a,b and b)
//                       <-----k-----> : 2 reads (b and c)
//                        <-----k-----> : 1 reads (c)
//                         <-----k-----> : 2 reads (c and d)
//                            ...
// the position i is contained into a kmer fully contained into only 1 mapped read, return 1
// for doing this we stored on each position of the fragment the number of k-mers starting at this position that fully belong to a read that was mapped.
int minimal_kmer_coverage(p_fragment_info the_starter, int read_file){
    
    int i, val_min=INT_MAX;
    const  int stopi=strlen(the_starter->w);
    for(i=0;i<stopi;i++){ // for each position on the read
        val_min=min(val_min, the_starter->read_coherent_positions[read_file][i]);
    }
    return val_min;
}

void set_read_coherent(p_fragment_info the_starter, const int min_coverage, int read_file){
    int i;
    // V1: the whole fragment has to be k_read coherent or V2 where the last k positions have no influence on the coherency of the fragment.
    // V2 is appropriate for the cases where the fragment is the end of a sequence (transcript, chromosome) and thus, no read are "longuer" than the sequence:
    //    ----------------- fragment
    //    °°°°°°°°°°°°        read
    //    °°°°°°°°°°°°     read
    //         °°°°°°°°°°°° read
    //         °°°°°°°°°°°° read
#ifdef KMER_SPANNING
    const int stop=strlen(the_starter->w)-size_seeds;
#else
    const int stop=strlen(the_starter->w);
#endif
    for(i=0;i<stop;i++) if(the_starter->read_coherent_positions[read_file][i]<min_coverage) {the_starter->read_coherent[read_file]=0; return;}
    the_starter->read_coherent[read_file]=1;
}


///////////////////////////////////////////////////////////////////////////////////////////////
///// we use a boolean vector to know if a read was alread mapped at one position on a fragment
///////////////////////////////////////////////////////////////////////////////////////////////
char mask [8];

char * new_boleanVector (long size){
    char* boolean_vector = (char*)malloc(size/8+1);
    test_alloc(boolean_vector);
    memset(boolean_vector, 0, size/8+1);
    mask[0]=1;  //00000001
    mask[1]=2;  //00000010
    mask[2]=4;  //00000100
    mask[3]=8;  //00001000
    mask[4]=16; //00010000
    mask[5]=32; //00100000
    mask[6]=64; //01000000
    mask[7]=128;//10000000
    return boolean_vector;
}


char is_boolean_vector_visited (const char * boolean_vector, long i){
    return (boolean_vector[i/8]&mask[i%8])==0?0:1;
}

void set_boolean_vector_visited(char * boolean_vector, long i){
    boolean_vector[i/8]|=mask[i%8];
}

void free_boolean_vector(char * boolean_vector){
    free(boolean_vector);
}

void reinit_boolean_vector(char * boolean_vector, long size){
    memset(boolean_vector, 0, size/8+1);
}




/**
 * Performs the first extension of the algorithm:
 * For each starter:
 *  - verify which reads maps on it with at most subst_allowed substitutions
 *  - for fragments that are fully covered by reads (each position is covered at least once):
 *    1) add their id in a list (that is returned by the function)
 *    2) detect the extending reads (at least "min_extension_coverage_depth" reads) that enable to extend right the starter
 *    3) store this extending reads in the structure of the fragments (fragment_info)
 
 * returns the average read length
 */

float read_coherence (FILE * reads_file,
                      char * reads_file_name,
                      const int k,
                      const int min_coverage,
                      p_fragment_info * all_starters,
                      const int number_starters,
                      int read_file,
                      int qual,
                      int nb_events_per_set,
                      int nb_fragment_per_event,
                      FILE * sam_out,
                      int subst_allowed,
                      const char no_subsutitution_on_central_position){
    
	//////////////////////////////////////////////////////////////////////////
	/////////////// read all reads - storing those coherent with reads ///////
	//////////////////////////////////////////////////////////////////////////
	//double cumulated_size_reads=0;
	char * seed = (char *) malloc (sizeof(char)*(k+1)); // fragment used for reading seeds
	char * starter;
    
    uint64_t offset_seed;
    uint64_t nb_seeds;
    
	// working variables
	int read_len, i, ii, j, starter_id,  pwi, stop, read_coherence;
    long int  read_number=0;
    
	// map of starter -> position (for each read and direction, stores the starter and position already tested.)
    char * boolean_vector = new_boleanVector(number_starters);
    //	hash_t already_tested_direct = hash_create(1000);
    //	hash_t already_tested_rev_comp = hash_create(1000);
    
    
	// book space for the read that going to be read.
    char * line =  (char *)malloc(sizeof(char)*1048576); // 1048576
	char * read = (char *)malloc(sizeof(char)*16384); //
	char * quality = (char *)malloc(sizeof(char)*16384);//
	test_alloc(read);
	test_alloc(quality);
//	char * rev_comp_read = (char *)malloc(sizeof(char)*16384); //
//	test_alloc(rev_comp_read);
//	char * rev_comp_quality = (char *)malloc(sizeof(char)*16384); //
//	test_alloc(rev_comp_quality);
//	char starter_key[50]; // USED as key for hash table.
    
	//number_of_reads=number_of_sequences_in_file(reads_file);
    
	//if(!silent){
	//	printf("\n");                                         // print progress
	//	for (i=0;i<100;i++) printf(" "); printf("] mapping reads\r[");       // print progress
	//}
    char is_fastq=0;
    if ( qual || strstr(reads_file_name,"fastq") || strstr(reads_file_name,"fq") || strstr(reads_file_name,"txt") ) is_fastq=1;
	// read all the reads
	do{
		//		int howmuchmapped=0; //DEB
		//if(!silent) if((++extention_number)%(1+(number_of_reads/100)) == 0) printf(".");  // print progress
		if(!silent) if((++read_number)%(100000) == 0) printf("\r %ld reads treated", read_number);  // print progress
		if ( is_fastq )
			read_len = get_next_sequence_for_fastq(reads_file, read, quality,line);
		else
			read_len = get_next_fasta_sequence(reads_file, read,line);
		if(read_len<1) break; // we have read all the file.
		//cumulated_size_reads+=read_len;
        
//		// find the reverse complement of the read
//		strcpy(rev_comp_read, read);
//		revcomp(rev_comp_read,read_len);
//		// if qual, find the reverse complement of the quality
//		if ( qual )
//		{
//			strcpy(rev_comp_quality, quality);
//			rev(rev_comp_quality,read_len);
//		}
        
		stop = read_len-k+1;
        //#ifdef DEBUG_MAPPING
        //		printf("read = %s\n", read);
        //#endif
		// read all seeds present on the read:
        int direction;
        for(direction=0;direction<2;direction++){ // try the two possible directions of the read
            for (i=0;i<stop;i++){
                for(j=0;j<k;j++) seed[j]=read[i+j]; seed[j]='\0';// read the seed
                
                // if the seed is indexed in the fragments:
                if(get_seed_info(seeds_count,seed,&offset_seed,&nb_seeds)){
                    // for each occurrence of this seed on the starter:
                    for (ii=offset_seed; ii<offset_seed+nb_seeds; ii++) {
                        //while(c!=NULL){
                        couple * value = &(seed_table[ii]);                        
                        // if the tuple read/starter/position was not already mapped nor already detected as unmappable, and far enought from previous mapping of the same read try to map.
//                        if(exists_tuple_read_position(already_tested_direct, starter_key, pwi)==0){
                        if(is_boolean_vector_visited(boolean_vector, value->a)==0){
                            set_boolean_vector_visited(boolean_vector, value->a);
                            starter = all_starters[value->a]->w;
                            //                        sprintf(starter_key, "%d", value->a);
                            pwi = value->b-i; // starting position of the read on the starter.
                            
                            // overview situation:
                            
                            //        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  starter
                            //        <---------> b
                            //                   [--------]                     seed
                            //             ******************************       read
                            //             <----> i
                            //        <---> pwi

                        
#ifdef DEBUG_MAPPING
                            printf("read = %s maps pos %d with %s (starter %d)\n", read, pwi, starter, value->a);
#endif
 
                            if(no_subsutitution_on_central_position)
                                read_coherence = read_coherent_SNP(pwi, starter, read, subst_allowed);
                            else
                                read_coherence = read_coherent_generic(pwi, starter, read, subst_allowed);

                            
                            if(read_coherence == 1){ // tuple read starter position is read coherent
                                
#ifdef DEBUG_MAPPING
                                printf("SUCCESS %d %d \n", pwi, value->a);
                                print_mapping(pwi,all_starters[value->a]->w ,read); //DEB
#endif
                                feed_coherent_positions(all_starters[value->a], pwi, (int)strlen(read), qual, quality, i, read_file);
                                if( sam_out ) fprintf(sam_out,"%s\t%s\tC%d\t%d\n", all_starters[value->a]->comment, read, read_file+1, pwi);
                                
#ifdef DEBUG_QUALITY
                                
                                if ( qual )
                                {
                                    int l;
                                    printf("\nstarter %d pwi=%d\n",value->a,pwi);
                                    for(l=0;l<70;l++)printf(" ");
                                    printf("%s\n",all_starters[value->a]->w);
                                    for(l=0;l<70+pwi;l++)printf(" ");
                                    printf("%s\n", read);
                                    for(l=0;l<70+pwi;l++)printf(" ");
                                    printf("%s\n", quality);
                                }
#endif
                            } // end tuple read starter position is read coherent
//                            hash_add_int_to_list(already_tested_direct, starter_key, pwi);
                        }
                        //c=c->prox;
                    }
                }
            } // end all seeds of the read
            revcomp(read,read_len);
            if ( qual ) rev (quality,read_len);
            reinit_boolean_vector(boolean_vector, number_starters);
        } // end both directions
        
      
	}while(1); // end all reads of the file
    if(!silent) printf("\r %ld reads treated\n", read_number);  // print progress
//	hash_delete(already_tested_direct, list_of_generic_free);
//	hash_delete(already_tested_rev_comp, list_of_generic_free);
    
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////// for each starter: check those fully coherent and store left and right reads covering them ///////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
//	if(!silent){
//		printf("]\n");                                         // print progress
//		for (i=0;i<100;i++) printf(" "); printf("] reads coherence detection\r[");       // print progress
//	}
    
	for (starter_id=0;starter_id < nb_events_per_set*nb_fragment_per_event;starter_id++){
        if(!silent && (starter_id%(1+(nb_events_per_set/50)) == 0)) printf(".");  // print progress
        set_read_coherent(all_starters[starter_id], min_coverage, read_file);
	} // end all fragments
    
    //	if(!silent) printf("]\n");  // print progress
	free(read);
    free(quality);
    free_boolean_vector(boolean_vector);
//	free(rev_comp_read);
    free(line);
    
	//return (float)(cumulated_size_reads/(float)number_of_reads);
	return 0; 
}



