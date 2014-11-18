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
 * fragment_index.c
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */

#include<fragment_index.h>
#include<fragment_info.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<list.h>
#include<commons.h>
#include<couple.h>
#include<hash.h>
#include <stdint.h>

//#define DEBUG_INDEXING

int line_num(FILE * f)
{
	rewind(f);
	char c;
	int lines = 0;
	while((c = fgetc(f)) != EOF) if(c == '\n') lines++;
	if(c != '\n') lines++;
	rewind(f);
	return lines;
}


void index_one_seed(const char * seed, const int fragment_id, const int position_on_fragment){
	hash_add_something_to_list(seeds,(char *)seed,create_couple(fragment_id,position_on_fragment));
}



////void index_one_seed_using_nodes(const char * seed, const p_node fragment_node, const int fragment_id, const int position_on_fragment){
//void index_one_seed_using_nodes(const char * seed, const p_node fragment_node, const int position_on_fragment){
//	node_couple * seed_info = create_node_couple(fragment_node, position_on_fragment);
//	hash_add_something_to_list(seeds,(char *)seed,seed_info);
//}





char *  strdup_upper_case(char * in){
    // count number of upper case letters in "in"
    int count =0;
    int i;
    for(i=0;i<strlen(in);i++) if(in[i]>='A' && in[i]<='Z') count++;
    char * temp = (char *) malloc(sizeof(char)*(count+1)); test_alloc(temp);
    int j=0;
    for(i=0;i<strlen(in);i++) if(in[i]>='A' && in[i]<='Z') temp[j++]=in[i];
    temp[j]='\0';
    return temp;
}


char * strdup_first_lower(char * in){
    // count number of first lower case letters in "in"
    int count =0;
    int i;
    for(i=0;i<strlen(in);i++)
        if(in[i]>='a' && in[i]<='z') count++;
        else break;
    char * temp = (char *) malloc(sizeof(char)*(count+1)); test_alloc(temp);
    int j=0;
    for(i=0;i<strlen(in);i++)
        if(in[i]>='a' && in[i]<='z') temp[j++]=in[i];
        else break;
    temp[j]='\0';
    return temp;
}
#include<assert.h>

char * strdup_last_lower(char * in){
    // count number of first lower case letters in "in"
    int count =0;
    int i;
    for(i=strlen(in)-1;i>=0;i--)
        if(in[i]>='a' && in[i]<='z') {
            count++;
        }
        else break;
    
    char * temp = (char *) malloc(sizeof(char)*(count+1)); test_alloc(temp);
    int j=count-1;
    for(i=strlen(in)-1;i>=0;i--)
        if(in[i]>='a' && in[i]<='z'){
            temp[j--]=in[i];
        }
        else break;
    temp[count]='\0';
    return temp;
}



// read and store all fragments presents in the pointed file.
// index by seeds of length k all these fragments.
// each fragment is stored twice: one direct, one reverse complement.
p_fragment_info * index_starters_from_input_file (const int k, int nb_events_per_set, const int nb_fragment_per_event, const char input_only_upper, const int index_stride){
	char * temp_fragment = (char *) calloc (sizeof(char)*1048576,1); // fragment used for reading line by line
	test_alloc(temp_fragment);
	char * temp_fragment2 = (char *) malloc (sizeof(char)*131072); // fragment used for reading line by line
	test_alloc(temp_fragment2);
	int witness;                                              // is a fragment was read ?
	kmer_type coded_seed;
	int i,z,stop;
    char * line = malloc(sizeof(char)*1048576);
    
    uint64_t total_seeds = 0 ;
    seeds = hash_create(100000);  	test_alloc(seeds);
    seeds_count = hash_create_binarykey(100000);  // todo  change to binary key (hash_t)AllocateHashTable(kmersize,1); //
    
	p_fragment_info * all_starters;                   // all existing starters are stored in this array.
    // this way any starter is accessed by its offset in this array = its id
	// allocate space for the main array of fragments
	all_starters = (p_fragment_info *) malloc(sizeof(p_fragment_info)*(nb_events_per_set*nb_fragment_per_event)); // fragment_index.h (
	test_alloc(all_starters);
#ifdef DEBUG_INDEXING
    printf("indexing %d*%d sequences, allocating memory for storing info about %d read sets\n", nb_fragment_per_event, nb_events_per_set, number_of_read_sets);
#endif
	// each fragment has an id. each seed pobints to couples (id, position). and the fragment is then found thanks to all_fragment[id];
	int fragment_id=0;
	do{
//        if(fragment_id%1000==0) printf("\r%d fragments stored", fragment_id);
		if ( fragment_id+1 > nb_events_per_set*nb_fragment_per_event ) break; // we read the starters we needed to read
		witness=get_next_sequence_and_comments_for_starters(temp_fragment, temp_fragment2, input_only_upper,line);
        
		if(witness<1) break;  // we have read all the file.
        
        
		// create the corresponding direct fragment
		all_starters[fragment_id] = (p_fragment_info) malloc(sizeof(fragment_info)); test_alloc(all_starters[fragment_id]);
        //#ifdef GET_ONLY_UPPER_CHARS
        if(input_only_upper){
            all_starters[fragment_id]->left_extension = strdup_first_lower(temp_fragment);
            //        printf("left ext = %s\n", all_starters[fragment_id]->left_extension); // DEB
            all_starters[fragment_id]->w = strdup_upper_case(temp_fragment);
            //        printf("center = %s\n", all_starters[fragment_id]->w); // DEB
            all_starters[fragment_id]->right_extension = strdup_last_lower(temp_fragment);
            //        printf("right ext = %s\n", all_starters[fragment_id]->right_extension); // DEB
        }
        else {
            //#else
            all_starters[fragment_id]->w = strdup(temp_fragment);  // the fragment is stored
            //        printf("center = %s\n", all_starters[fragment_id]->w); // DEB
        }
        
        
        //#endif
        all_starters[fragment_id]->mapped_with_current_read = (char *)malloc(sizeof(char)*number_of_read_sets);test_alloc(all_starters[fragment_id]->mapped_with_current_read);
        all_starters[fragment_id]->tested_pwis_with_current_read = (listint **)malloc(sizeof(listint *)*number_of_read_sets); test_alloc(all_starters[fragment_id]->tested_pwis_with_current_read);
		all_starters[fragment_id]->read_coherent = (char*) malloc(sizeof(char)*number_of_read_sets);test_alloc(all_starters[fragment_id]->read_coherent);
		all_starters[fragment_id]->number_mapped_reads = (int*) malloc(sizeof(int)*number_of_read_sets);test_alloc(all_starters[fragment_id]->number_mapped_reads);
		all_starters[fragment_id]->read_coherent_positions = (unsigned char**) malloc(sizeof(unsigned char*)*number_of_read_sets);test_alloc(all_starters[fragment_id]->read_coherent_positions);
#ifdef CHARQUAL
		all_starters[fragment_id]->sum_quality_per_position = (unsigned char**) malloc(sizeof(unsigned char*)*number_of_read_sets);test_alloc(all_starters[fragment_id]->sum_quality_per_position);
#else
        all_starters[fragment_id]->sum_quality_per_position = (int **) malloc(sizeof(int*)*number_of_read_sets); test_alloc(all_starters[fragment_id]->sum_quality_per_position);
        
#endif
		all_starters[fragment_id]->comment = format_comment(temp_fragment2);
        
		for (i=0; i<number_of_read_sets; i++)
		{
            all_starters[fragment_id]->mapped_with_current_read[i] = 0;
            all_starters[fragment_id]->tested_pwis_with_current_read[i] = listint_create();
			all_starters[fragment_id]->read_coherent_positions[i] = (unsigned char *) malloc (strlen(all_starters[fragment_id]->w)*sizeof(unsigned char)); test_alloc(all_starters[fragment_id]->read_coherent_positions[i]);
#ifdef CHARQUAL
			all_starters[fragment_id]->sum_quality_per_position[i] = (unsigned char *) malloc (strlen(all_starters[fragment_id]->w)*sizeof(unsigned char)); test_alloc(all_starters[fragment_id]->sum_quality_per_position[i]);
#else
            all_starters[fragment_id]->sum_quality_per_position[i] = (int *) malloc (strlen(all_starters[fragment_id]->w)*sizeof(int)); test_alloc(all_starters[fragment_id]->sum_quality_per_position[i]);
            
#endif
			for(z=0;z<strlen(all_starters[fragment_id]->w); z++) all_starters[fragment_id]->read_coherent_positions[i][z]=0;
			for(z=0;z<strlen(all_starters[fragment_id]->w); z++) all_starters[fragment_id]->read_coherent_positions[i][z]=0;
			for(z=0;z<strlen(all_starters[fragment_id]->w); z++) all_starters[fragment_id]->sum_quality_per_position[i][z]=0;
            
            
			all_starters[fragment_id]->number_mapped_reads[i]=0;
		}
        
#ifdef INPUT_FROM_KISSPLICE
        if(strstr(all_starters[fragment_id]->comment,"upper")) all_starters[fragment_id]->upperpath=1; // even numbers are upper.
        else all_starters[fragment_id]->upperpath=0;
        
        if(all_starters[fragment_id]->upperpath){
            all_starters[fragment_id]->nb_reads_fully_in_S = (int *) malloc(sizeof(int)*number_of_read_sets); test_alloc(all_starters[fragment_id]->nb_reads_fully_in_S);
            all_starters[fragment_id]->nb_reads_overlapping_AS = (int *) malloc(sizeof(int)*number_of_read_sets); test_alloc(all_starters[fragment_id]->nb_reads_overlapping_AS);
            all_starters[fragment_id]->nb_reads_overlapping_SB = (int *) malloc(sizeof(int)*number_of_read_sets); test_alloc(all_starters[fragment_id]->nb_reads_overlapping_SB);
        }
        all_starters[fragment_id]->nb_reads_overlapping_both_AS_and_SB = (int *) malloc(sizeof(int)*number_of_read_sets); test_alloc(all_starters[fragment_id]->nb_reads_overlapping_both_AS_and_SB);
        
        
        for (i=0; i<number_of_read_sets; i++)
		{
            if(all_starters[fragment_id]->upperpath){
                all_starters[fragment_id]->nb_reads_fully_in_S[i] = 0;
                all_starters[fragment_id]->nb_reads_overlapping_AS[i] = 0;
                all_starters[fragment_id]->nb_reads_overlapping_SB[i] = 0;
            }
            all_starters[fragment_id]->nb_reads_overlapping_both_AS_and_SB[i] = 0;
        }
#endif
		
#ifdef DEBUG_INDEXING
		printf("counting in %s\n", all_starters[fragment_id]->w);
#endif
		// read all the seeds present on the fragment
		stop=strlen(all_starters[fragment_id]->w)-k+1;
//        char validSeed;
		for (i=0;i<stop;i+= index_stride){
//            validSeed=1;
//            for(j=0;j<k && validSeed==1;j++)// read all characters to check if its a valid seed
//                if(!valid_character(all_starters[fragment_id]->w[j+i])){ // if i+j character is not valid
//                    i+=j; // don't test the next j+1 next positions (+1 will come from the 'for' loop)
//                    validSeed=0;
//                }
//            if(validSeed){
                coded_seed=codeSeed(all_starters[fragment_id]->w+i); // init the seed (as seeds are not consecutives
                hash_incr_kmer_count(seeds_count,&coded_seed);
                total_seeds++;
//            }
			
		}
		fragment_id++;
        
	}while (1);
    
#ifdef DEBUG_INDEXING
    printf("%d seeds\n", total_seeds);
#endif
    
//    printf("\r%d fragments stored\n", fragment_id);
    
    
    seed_table  = calloc(total_seeds,sizeof(couple));
    
    
    iterate_and_fill_offsets(seeds_count);
    
    
    
    
    ///second loop over fragments  : create the index
    fragment_id=0;
	do{
        
		if ( fragment_id+1 > nb_events_per_set*nb_fragment_per_event ) break; // we read the starters we needed to read
        
#ifdef DEBUG_INDEXING
		printf("indexing in %s\n", all_starters[fragment_id]->w);
#endif
		// read all the seeds present on the fragment
		stop=strlen(all_starters[fragment_id]->w)-k+1;
//        char validSeed=1;
		for (i=0;i<stop;i+= index_stride){
//            validSeed=1;
//            for(j=0;j<k && validSeed==1;j++)// read all characters to check if its a valid seed
//                if(!valid_character(all_starters[fragment_id]->w[j+i])){ // if i+j character is not valid
//                    validSeed=0;
//                    i+=j; // don't test the next j+1 next positions (+1 will come from the 'for' loop)
//                    continue;
//                }
//            if(validSeed){
                coded_seed=codeSeed(all_starters[fragment_id]->w+i); // init the seed
                hash_fill_kmer_index(seeds_count,&coded_seed,seed_table, fragment_id, i);
//            }
			
		}
		fragment_id++;
        
	}while (1);
    
	free(temp_fragment);
    
    free(line);
    return all_starters;
}


void free_seeds_index (){
	hash_delete(seeds, list_of_generic_free);
}
