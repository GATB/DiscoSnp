/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: P.Peterlongo, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

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

char prefix(const char *pre, const char *str)
{
    return strncmp(pre, str, strlen(pre)) == 0;
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
        test_alloc(all_starters[fragment_id]->mapped_with_current_read);
        all_starters[fragment_id]->tested_pwis_with_current_read = (listint **)malloc(sizeof(listint *)*number_of_read_sets); test_alloc(all_starters[fragment_id]->tested_pwis_with_current_read);
        test_alloc(all_starters[fragment_id]->tested_pwis_with_current_read)
		all_starters[fragment_id]->read_coherent = (char*) malloc(sizeof(char)*number_of_read_sets);test_alloc(all_starters[fragment_id]->read_coherent);
        test_alloc(all_starters[fragment_id]->read_coherent);
		all_starters[fragment_id]->number_mapped_reads = (int*) malloc(sizeof(int)*number_of_read_sets);test_alloc(all_starters[fragment_id]->number_mapped_reads);
        test_alloc(all_starters[fragment_id]->number_mapped_reads);
		all_starters[fragment_id]->read_coherent_positions = (unsigned char**) malloc(sizeof(unsigned char*)*number_of_read_sets);test_alloc(all_starters[fragment_id]->read_coherent_positions);
        test_alloc(all_starters[fragment_id]->read_coherent_positions);
		all_starters[fragment_id]->sum_qualities = (unsigned int*) malloc(sizeof(unsigned int)*number_of_read_sets);
        test_alloc(all_starters[fragment_id]->sum_qualities);
        all_starters[fragment_id]->nb_mapped_qualities = (unsigned int*) malloc(sizeof(unsigned int)*number_of_read_sets);
        test_alloc(all_starters[fragment_id]->nb_mapped_qualities);
		all_starters[fragment_id]->comment = format_comment(temp_fragment2);
        all_starters[fragment_id]->nbOfSnps = 0;
        if (prefix("SNP",all_starters[fragment_id]->comment)) {
            all_starters[fragment_id]->nbOfSnps=1; // We don't know yep how many, at least one.
        }
        
		for (i=0; i<number_of_read_sets; i++)
		{
            all_starters[fragment_id]->mapped_with_current_read[i] = 0;
            all_starters[fragment_id]->tested_pwis_with_current_read[i] = listint_create();
			all_starters[fragment_id]->read_coherent_positions[i] = (unsigned char *) malloc (strlen(all_starters[fragment_id]->w)*sizeof(unsigned char)); test_alloc(all_starters[fragment_id]->read_coherent_positions[i]);
            test_alloc(all_starters[fragment_id]->read_coherent_positions[i])

			for(z=0;z<strlen(all_starters[fragment_id]->w); z++) all_starters[fragment_id]->read_coherent_positions[i][z]=0;
			for(z=0;z<strlen(all_starters[fragment_id]->w); z++) all_starters[fragment_id]->read_coherent_positions[i][z]=0;
            
            
            all_starters[fragment_id]->nb_mapped_qualities[i]=0;
            all_starters[fragment_id]->sum_qualities[i]=0;
			all_starters[fragment_id]->number_mapped_reads[i]=0;
            
		}
    		
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
    
    
    
    
    ///third loop over fragments : for SNPs, store the SNP positions and the SNP qualities
    fragment_id=0;
	do{
        
		if ( fragment_id+1 > nb_events_per_set*nb_fragment_per_event ) break; // we read the starters we needed to read
        if ( all_starters[fragment_id]->nbOfSnps==0 ) {fragment_id+=2; continue;} // applies only for SNPs
        
        char * seq1 = all_starters[fragment_id]->w;
        char * seq2 = all_starters[fragment_id+1]->w;
        int size_seq = strlen(seq1);
        assert(size_seq == strlen(seq2));
        
        // compute the number of SNPs:
        int local_number_of_SNPs=0;
        for (i=0; i<size_seq; i++) {
            if (seq1[i]!=seq2[i]) {
                local_number_of_SNPs++;
            }
        }
        all_starters[fragment_id]->nbOfSnps=local_number_of_SNPs;
        all_starters[fragment_id]->SNP_positions = (char *) malloc (sizeof(char)*(local_number_of_SNPs+1)); // add a dummy SNP
        test_alloc(all_starters[fragment_id]->SNP_positions);
        
        // we do not fill the all_starters[fragment_id+1] with the same information
        local_number_of_SNPs=0;
        for (i=0; i<size_seq; i++) {
            if (seq1[i]!=seq2[i]) {
                all_starters[fragment_id]->SNP_positions[local_number_of_SNPs++]=i;
            }
        }
        all_starters[fragment_id]->SNP_positions[local_number_of_SNPs] = size_seq+1; // DUMMY SNP
		fragment_id+=2;
	}while (1);
    
    
	free(temp_fragment);
    
    free(line);
    return all_starters;
}


void free_seeds_index (){
	hash_delete(seeds, list_of_generic_free);
}
