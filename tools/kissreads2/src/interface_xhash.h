/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2020  INRIA
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
 * hashing used by kissreads2
 * this is a generic file, which defines an interface between kissreads2 and a hash table implementation
 */

#include <stddef.h>
#include <unistd.h>
#include <stdint.h>
#include "commons.h"
#include <xhash.h>
#ifndef _HASH_H
#define _HASH_H

/*
 * generic types: we typecast in hash_*.c functions for specific implementations
 */



typedef uint64_t hash_val;


/*
 * basic functions:
 *  - create hash table
 *  - insert element into table
 *  - return element from table
 *  - delete table
 */

extern xhash xhash_create_seed_index();

void hash_incr_kmer_count(xhash * map, const kmer_type * key, GlobalValues& gv);

void iterate_and_fill_offsets( xhash * map, GlobalValues &gv );

int get_seed_info(xhash * map, const kmer_type * key, uint64_t * offset_seed, uint64_t * nb_seeds, GlobalValues &gv );

void hash_fill_kmer_index(xhash * map, const kmer_type * key, std::pair <uint64_t, int > * seed_table, const uint64_t fragment_id, const int position_on_fragment, GlobalValues &gv);

 void get_offset_and_nb_from_sinfo(hash_val  sinfo, uint64_t & offset_seed,uint64_t & nb_seeds, GlobalValues &gv  );

void set_offset_and_nb_into_sinfo(hash_val * sinfo, uint64_t  offset_seed,uint64_t  nb_seeds, GlobalValues &gv  );


#endif                          /* _HASH_H */
