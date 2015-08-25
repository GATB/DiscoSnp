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

//void hash_add_int_to_list(hash_t map, char * key, int value);

void hash_add_something_to_list(hash_t map, const char * key, void *something);

void hash_incr_kmer_count(hash_t map, const kmer_type * key, GlobalValues& gv);

void iterate_and_fill_offsets(hash_t map, GlobalValues &gv );

//int get_seed_info(hash_t map, const char * key, uint64_t * offset_seed, uint64_t * nb_seeds);
int get_seed_info(hash_t map, const kmer_type * key, uint64_t * offset_seed, uint64_t * nb_seeds, GlobalValues &gv );
void hash_fill_kmer_index(hash_t map, const kmer_type * key, couple * seed_table, const int fragment_id, const int position_on_fragment, GlobalValues &gv );

void * hash_list_sorted_by_positions(hash_t map);
 void get_offset_and_nb_from_sinfo(hash_val  sinfo, uint64_t & offset_seed,uint64_t & nb_seeds, GlobalValues &gv  );
 void set_offset_and_nb_into_sinfo(hash_val * sinfo, uint64_t  offset_seed,uint64_t  nb_seeds, GlobalValues &gv  );


#endif                          /* _HASH_H */
