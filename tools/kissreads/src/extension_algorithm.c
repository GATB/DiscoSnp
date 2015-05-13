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

//
//// checks if, for a fragment, a tuple read/position mapped already (return 1) or none (return 0)
//int exists_tuple_read_position(const hash_t map, const char * starter, const int pos){
//	list * read_coherent_positions;                 //for one couple (read/fragment) : list of coherent mapped positions (usually, one position)
//	if(hash_entry_by_key(map, starter, (void **) (&read_coherent_positions)) !=0){
//		cell *c = ((list *)read_coherent_positions)->first;
//		while(c != NULL){
//			if((*(int *) (c->val))== pos){
//				return 1;}
//			c = c->prox;
//		}
//	}
//	return 0;
//}

void print_mapping(const int pos_on_fragment, const char * fragment, const char *read){
	int i;
	printf("pos: %d\n", pos_on_fragment);
	for(i=0;i<60;i++) printf(" ");
	printf("%s\n", fragment);
    
	for(i=0;i<60+pos_on_fragment;i++) printf(" ");
	printf("%s\n", read);
    
}



void feed_coherent_positions(p_fragment_info * starters, const int starter_id, const int start, const int length_read, const int qual, char *quality, int start_on_read, int read_file_id){
    
    
    int i, start_on_starter, stop_on_starter;
	if(start<0) start_on_starter=0;
    else start_on_starter=start;
    
    p_fragment_info the_starter=starters[starter_id];
    p_fragment_info the_reference_starter = starters[2*(starter_id/2)]; // In case of snps, only the upper path starter contains informations such as the positions of the SNPs. This is the reference?
	
    if(start+length_read<strlen(the_starter->w)) stop_on_starter=start+length_read;
    else stop_on_starter=strlen(the_starter->w);
    
    
	the_starter->number_mapped_reads[read_file_id]++;
    
    if ( qual ){
        if (the_reference_starter->nbOfSnps>0) {
            //            printf("there are %d SNPs\n", the_reference_starter->nbOfSnps);
            int snp_id;
            for(snp_id=0;snp_id<the_reference_starter->nbOfSnps;snp_id++){ // we only add the qualities of the mapped SNPs
                //                printf("the_reference_starter->SNP_positions[%d]=%d \n", snp_id, the_reference_starter->SNP_positions[snp_id]); //DEB
                i=the_reference_starter->SNP_positions[snp_id];
                //                printf(" %d \n", start_on_read + i - start_on_starter); //DEB
                the_starter->sum_qualities[read_file_id] += (unsigned int) quality[start_on_read + i - start_on_starter];
                the_starter->nb_mapped_qualities[read_file_id] += 1;
            }
        }
        else{ // we sum all qualities and divide by the number of positions
            int sum_temp=0;
            int denom=0;
            for(i=start_on_starter;i<stop_on_starter;i++) { // to avoid to increase too much the the_starter->sum_qualities array, we add the average quality of the whole read.
                denom+=1;
                sum_temp+=(unsigned int) quality[start_on_read + i - start_on_starter];
            }
            if(denom>0){
                the_starter->sum_qualities[read_file_id] += (unsigned int)sum_temp/denom;
                the_starter->nb_mapped_qualities[read_file_id] += 1;
            }
        }
    }
    
    
    
#ifdef KMER_SPANNING
    // the position i is contained into a kmer fully contained into only 1 mapped read, return 1
    // for doing this we stored on each position of the fragment the number of k-mers starting at this position that fully belong to a read that was mapped
	
    //  -------------------------------------------------------------- starter
    //            ************************
    //                       <-----k----->
    //  00000000001111111111110000000000000000000000000000000000000000 the_starter->read_coherent_positions[read_file_id]
	if(start+length_read-minimal_read_overlap<strlen(the_starter->w)) stop_on_starter=start+length_read-minimal_read_overlap;
    else stop_on_starter=strlen(the_starter->w);
#endif
    
	for(i=start_on_starter;i<stop_on_starter;i++) Sinc8(the_starter->read_coherent_positions[read_file_id][i]);
	
    
    
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
int minimal_kmer_coverage(p_fragment_info the_starter, int read_file_id){
    
    int i, val_min=INT_MAX;
    const  int stopi=strlen(the_starter->w);
    for(i=0;i<stopi;i++){ // for each position on the read
        val_min=min(val_min, the_starter->read_coherent_positions[read_file_id][i]);
    }
    return val_min;
}

void set_read_coherent(p_fragment_info the_starter, const int min_coverage, int read_file_id){
    int i;
    // V1: the whole fragment has to be k_read coherent or V2 where the last k positions have no influence on the coherency of the fragment.
    // V2 is appropriate for the cases where the fragment is the end of a sequence (transcript, chromosome) and thus, no read are "longer" than the sequence:
    //    ----------------- fragment
    //    °°°°°°°°°°°°        read
    //    °°°°°°°°°°°°     read
    //         °°°°°°°°°°°° read
    //         °°°°°°°°°°°° read
#ifdef KMER_SPANNING
    const int stop=strlen(the_starter->w)-minimal_read_overlap;
    if(stop<=0){
        if(the_starter->read_coherent_positions[read_file_id][0]<min_coverage) {the_starter->read_coherent[read_file_id]=0; return;}
        the_starter->read_coherent[read_file_id]=1; return;
    }
#else
    const int stop=strlen(the_starter->w);
#endif
    //    printf("%d %d %d \n", strlen(the_starter->w), minimal_read_overlap, stop); //DEB
    for(i=0;i<stop;i++) if(the_starter->read_coherent_positions[read_file_id][i]<min_coverage) {the_starter->read_coherent[read_file_id]=0; return;}
    the_starter->read_coherent[read_file_id]=1;
}

void set_read_coherent_paired(p_fragment_info the_starter, const int min_coverage, int read_file_id){
    int i;
    // V1: the whole fragment has to be k_read coherent or V2 where the last k positions have no influence on the coherency of the fragment.
    // V2 is appropriate for the cases where the fragment is the end of a sequence (transcript, chromosome) and thus, no read are "longer" than the sequence:
    //    ----------------- fragment
    //    °°°°°°°°°°°°        read
    //    °°°°°°°°°°°°     read
    //         °°°°°°°°°°°° read
    //         °°°°°°°°°°°° read
#ifdef KMER_SPANNING
    const int stop=strlen(the_starter->w)-minimal_read_overlap;
    if(stop<=0){
        if(the_starter->read_coherent_positions[read_file_id][0]+the_starter->read_coherent_positions[read_file_id+1][0]<min_coverage) {
            the_starter->read_coherent[read_file_id]=0;
            the_starter->read_coherent[read_file_id+1]=0;
            return;}
        the_starter->read_coherent[read_file_id]=1;
        the_starter->read_coherent[read_file_id+1]=1;
        return;
    }
#else
    const int stop=strlen(the_starter->w);
#endif
    //    printf("%d %d %d \n", strlen(the_starter->w), minimal_read_overlap, stop); //DEB
    for(i=0;i<stop;i++)
        if(the_starter->read_coherent_positions[read_file_id][i]+the_starter->read_coherent_positions[read_file_id+1][i]<min_coverage) {
            the_starter->read_coherent[read_file_id]=0;
            the_starter->read_coherent[read_file_id+1]=0;
            return;}
    the_starter->read_coherent[read_file_id]=1;
}


///////////////////////////////////////////////////////////////////////////////////////////////
///// we use a boolean vector to know if a read was alread mapped on a fragment
// OCTOBER 2014: THIS SOLUTION IS NOT GOOD: some read maps several position of a fragment
// (more and more true with long reads)
// The new way of checking if a triplet (fragment/read/pwi) is already tested uses the fragment_info structure.
// Each fragment info houses a listint where tested pwi values are stored for a read.
// Each of these lists is emptied after each read mapping.
///////////////////////////////////////////////////////////////////////////////////////////////
//char mask [8];
//
//char * new_boleanVector (long size){
//    char* boolean_vector = (char*)malloc(size/8+1);
//    test_alloc(boolean_vector);
//    memset(boolean_vector, 0, size/8+1);
//    mask[0]=1;  //00000001
//    mask[1]=2;  //00000010
//    mask[2]=4;  //00000100
//    mask[3]=8;  //00001000
//    mask[4]=16; //00010000
//    mask[5]=32; //00100000
//    mask[6]=64; //01000000
//    mask[7]=128;//10000000
//    return boolean_vector;
//}
//
//
//char is_boolean_vector_visited (const char * boolean_vector, long i){
//    return (boolean_vector[i/8]&mask[i%8])==0?0:1;
//}
//
//void set_boolean_vector_visited(char * boolean_vector, long i){
//    boolean_vector[i/8]|=mask[i%8];
//}
//
//void free_boolean_vector(char * boolean_vector){
//    free(boolean_vector);
//}
//
//void reinit_boolean_vector(char * boolean_vector, long size){
//    memset(boolean_vector, 0, size/8+1);
//}


void map_all_reads_from_a_file(gzFile read_file,
                               char * read_file_name,
                               const int k,
                               p_fragment_info * all_starters,
                               int read_file_id,
                               int qual,
                               FILE * sam_out,
                               int subst_allowed,
                               const int minimal_read_overlap){
    
    //////////////////////////////////////////////////////////////////////////
	/////////////// read all reads - storing those coherent with reads ///////
	//////////////////////////////////////////////////////////////////////////
	char * starter;
    uint64_t offset_seed;
    uint64_t nb_seeds;
    
	// working variables
	int read_len, i, ii,  pwi, stop, read_coherence;
    long int  read_number=0;
    
	// map of starter -> position (for each read and direction, stores the starter and position already tested.)
	// book space for the read that going to be read.
    
	
    
    char * line =  (char *)malloc(sizeof(char)*1048576); // 1048576
	char * read = (char *)malloc(sizeof(char)*16384); //
	char * quality = (char *)malloc(sizeof(char)*16384);//
    
    test_alloc(line);
	test_alloc(read);
	test_alloc(quality);
    
    listint * tested_starters = listint_create(); // given a read, stores all starter ids on which a seed was seen. Enables to free quickly list of tested pwis of each tested starters.
	
    char is_fastq=0;
    
    if ( qual || strstr(read_file_name,"fastq") || strstr(read_file_name,"fq") || strstr(read_file_name,"txt") ) is_fastq=1;
    
    
	// read all the reads
	do{
		if(!silent) if((++read_number)%(100000) == 0) printf("\r %ld reads treated", read_number);  // print progress
		if ( is_fastq )
			read_len = get_next_sequence_for_fastq(read_file, read, quality,line);
		else
			read_len = get_next_fasta_sequence(read_file, read,line);
		if(read_len<0) break; // we have read all the file.
        
        const int minimal_pwi = minimal_read_overlap - read_len;
        // The read must overlap the fragment with at least minimal_read_overlap positions.
        // here is the first position on which the read may map :
        //        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  starter
        //     **********       read (10)
        //        <-----> minimal_read_overlap (7)
        //     <-> -pwi (-3)
        // -pwi+minimal_read_overlap <= |read|
        // -pwi <= |read|-minimal_read_overlap
        // pwi >= minimal_read_overlap-|read|
        // pwi >= 7-10 = -3
        // minimal_pwi = minimal_read_overlap-|read|
        
        
		stop = read_len-k+1;
#ifdef DEBUG_MAPPING
        printf("new read = X%sX Y%sY\n", line, read);
#endif
		// read all seeds present on the read:
        int direction;
        kmer_type coded_seed;
        char toinit=1;
        //        char validSeed;
        for(direction=0;direction<2;direction++){ // try the two possible directions of the read
            toinit=1; // we have to init a new seed
            for (i=0;i<stop;i++){ // for all possible seed on the read
                if(toinit) {
                    coded_seed=codeSeed(read+i); // init the seed
                    toinit=0;
                }
                else { // previous seed was correct, we extend it.
                    coded_seed=updateCodeSeed(read+i,&coded_seed); // utpdate the previous seed with a
                }
                if(get_seed_info(seeds_count,&coded_seed,&offset_seed,&nb_seeds)){
                    
                    
                    // for each occurrence of this seed on the starter:
                    for (ii=offset_seed; ii<offset_seed+nb_seeds; ii++) {
                        couple * value = &(seed_table[ii]);
                        //                        printf("value->a = %d\n", value->a); //DEB
                        starter = all_starters[value->a]->w;
                        if(all_starters[value->a]->mapped_with_current_read[read_file_id]==1) {
                            //                            printf("already mapped %d with this read\n", value->a); //DEB
                            continue; // a mach was already found with this starter, we consider only one match per couple (starter/read)
                        }
                        //                        sprintf(starter_key, "%d", value->a);
                        pwi = value->b-i; // starting position of the read on the starter.
                        //                        printf("a=%d pwi = %d\n", value->a,pwi); // DEB
                        if (listint_contains(all_starters[value->a]->tested_pwis_with_current_read[read_file_id],pwi)) {
                            //                            printf("already tested %d position %d with this read\n", value->a, pwi); //DEB
                            continue; // this reads was already (unsuccessfuly) tested with this starter at this position. No need to try it again.
                        }
                        if (numberInListint (all_starters[value->a]->tested_pwis_with_current_read[read_file_id])==0){ // this is the first time we meet this starter with this read, we add it in the list of met starters (in order to free the pwis lists of each encounted starters for this read)
                            listint_add(tested_starters,value->a);
                            //                            printf("adding %d in the list of tested starters with this read \n", *p_valuea); //DEB
                        }
                        listint_add(all_starters[value->a]->tested_pwis_with_current_read[read_file_id],pwi); // We store the fact that this read was already tested at this position on this starter.
                        
                        
                        
                        // overview general situation:
                        
                        //        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  starter
                        //        <---------> b
                        //                   [--------]                     seed
                        //             ******************************       read
                        //             <----> i
                        //        <---> pwi
                        
                        if (pwi<minimal_pwi) {
                            continue; // this read to not overlap enough with the starter.
                        }
                        const int maximal_pwi = strlen(starter)-minimal_read_overlap;
                        if (pwi > maximal_pwi) {
                            continue; // this read to not overlap enough with the starter.
                        }
                        //        ;;;;;;;;;;;  starter (11)
                        //             ******************************       read
                        //             <----> minimal_read_overlap (6)
                        //        <---> pwi (5)
                        // |starter| <= pwi+minimal_read_overlap
                        
                        
                        read_coherence=0;
                        if(all_starters[value->a]->nbOfSnps>0)
                            read_coherence = read_coherent_SNP(pwi, starter, read, subst_allowed, all_starters[value->a-value->a%2]->SNP_positions);
                        //value->a-value->a%2: the smallest id of the the fragments of this bubble (3-->2, 4-->4)
                        else
                            read_coherence = read_coherent_generic(pwi, starter, read, subst_allowed);
                        
                        
                        if(read_coherence == 1){ // tuple read starter position is read coherent
                            
                            all_starters[value->a]->mapped_with_current_read[read_file_id]=1;
                            
                            
#ifdef DEBUG_MAPPING
                            printf("SUCCESS %d %d \n", pwi, value->a);
                            print_mapping(pwi,all_starters[value->a]->w ,read); //DEB
#endif
                            feed_coherent_positions(all_starters, value->a , pwi, (int)strlen(read), qual, quality, i, read_file_id);
                            if( sam_out ) fprintf(sam_out,"%d %d %s\t%s\tC%d\t%d\n", value->a, value->b, all_starters[value->a]->comment, read, read_file_id+1, pwi);
                            
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
                        //c=c->prox;
                    }
                } // end all infos for the current seed
            } // end all seeds of the read
            revcomp(read,read_len);
            if ( qual ) rev (quality,read_len);
            
            // free the lists of tested positions:
            cellint * tested_starter=tested_starters->first;
            while (tested_starter != NULL)
            {
                listint_empty(all_starters[tested_starter->val]->tested_pwis_with_current_read[read_file_id]);
                all_starters[tested_starter->val]->mapped_with_current_read[read_file_id] = 0;
                //                printf("freeing starter %d\n", *(int *)tested_starter->val); //DEB
                tested_starter=tested_starter->prox;
            }
            listint_empty(tested_starters);
            
            
        } // end both directions
        
        
	}while(1); // end all reads of the file
    
    if(!silent) printf("\r %ld reads treated\n", read_number);  // print progress
    
    free(read);
    free(quality);
    free(line);
    listint_free(tested_starters);
    
    
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
float read_mapping (char ** read_file_names,
                    int read_set_id,
                    const int k,
                    const int min_coverage,
                    p_fragment_info * all_starters,
                    int qual,
                    FILE * sam_out,
                    int subst_allowed,
                    const int minimal_read_overlap
                    ){
    
    
	int p=0;
    gzFile read_file=gzopen(read_file_names[read_set_id+p],"r");
    if(read_file == NULL){		fprintf(stderr,"cannot open reads file %s, exit\n",read_file_names[read_set_id+p]);		exit(1);	}
    map_all_reads_from_a_file(read_file, read_file_names[read_set_id+p],
                              k,
                              all_starters,
                              read_set_id,
                              qual,
                              sam_out,
                              subst_allowed,
                              minimal_read_overlap);
    gzclose(read_file);
    
	
    
    return 0;
}


void set_read_coherency(p_fragment_info * all_starters, const int nb_events_per_set, const int nb_fragment_per_event, const char paired, const int min_coverage, const int read_set_id){
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////// for each starter: check those fully coherent and store left and right reads covering them ///////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int starter_id;
	for (starter_id=0;starter_id < nb_events_per_set*nb_fragment_per_event;starter_id++){
        if (!paired)
            set_read_coherent(all_starters[starter_id], min_coverage, read_set_id);
        else
            set_read_coherent_paired(all_starters[starter_id], min_coverage, read_set_id);
	} // end all fragments
    
	
	
    
    
}




