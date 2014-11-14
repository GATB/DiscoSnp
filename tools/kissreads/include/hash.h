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
 * hashing used by mapsembler
 * this is a generic file, which defines an interface between mapsembler and a hash table implementation
 */

#include <stddef.h>
#include <unistd.h>
#include <stdint.h>
#include "commons.h"
#ifndef _HASH_H
#define _HASH_H

/*
 * generic types: we typecast in hash_*.c functions for specific implementations
 */

typedef unsigned int *hash_t;
typedef long* hash_iter;


//typedef struct {
//	unsigned int offset_seed;
//	unsigned int nb_seeds;
//} seedinfo, * p_seedinfo;

typedef uint64_t hash_val;
#define NBITS_OFFSET_SEED 40




#include "couple.h"


/*
 * basic functions:
 *  - create hash table
 *  - insert element into table
 *  - return element from table
 *  - delete table
 */

extern hash_t hash_create(unsigned int nbuckets);
extern hash_t hash_create_binarykey(unsigned int nbuckets);

extern int hash_insert(hash_t map, const char *key, const void * data, size_t len);

extern ssize_t hash_entry_by_key(hash_t map, const char *key,
		void **data);

extern int hash_delete(hash_t map, void (*specific_free)(const void *));

void hash_clear(hash_t map, void (*specific_free)(const void *));

/*
 * iterator functions:
 *  - start iterator
 *  - return current element from iterator
 *  - test for end of iterator
 */

extern hash_iter hash_iterator_init(hash_t map);

extern ssize_t hash_iterator_return_entry(hash_t map, hash_iter iter,
		char **key, void **data);

extern hash_iter hash_iterator_next(hash_t map, hash_iter iter);

extern hash_iter hash_iterator_is_end(hash_t map, hash_iter iter);

/*
 * more advanced functions:
 *  - search for an element having specific value in the hash table
 *  - maintain a list inside each element of the hash table
 *  - return a list of sorted reads by position
 */
void hash_increase_int_value_if_exists(hash_t map, char * key);

int* hash_get_int_value(hash_t map, char * key);

void hash_set_int_to_key(hash_t map, char * key, int value);

extern ssize_t hash_search_key_int_value(hash_t map, const char *key, const int value_to_search);

void hash_add_int_to_list(hash_t map, char * key, int value);

void hash_add_something_to_list(hash_t map, const char * key, void *something);

void hash_incr_kmer_count(hash_t map, const kmer_type * key);

void iterate_and_fill_offsets(hash_t map);

//int get_seed_info(hash_t map, const char * key, uint64_t * offset_seed, uint64_t * nb_seeds);
int get_seed_info(hash_t map, const kmer_type * key, uint64_t * offset_seed, uint64_t * nb_seeds);
void hash_fill_kmer_index(hash_t map, const kmer_type * key, couple * seed_table, const int fragment_id, const int position_on_fragment);

void * hash_list_sorted_by_positions(hash_t map);
 void get_offset_and_nb_from_sinfo(hash_val  sinfo, uint64_t * offset_seed,uint64_t * nb_seeds );
 void set_offset_and_nb_into_sinfo(hash_val * sinfo, uint64_t  offset_seed,uint64_t  nb_seeds );


#endif                          /* _HASH_H */
