/*****************************************************************************
 *   DiscoMore: discovering polymorphism from raw unassembled NGS reads
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
 * fragment_index.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef FRAGMENT_INDEX_H_
#define FRAGMENT_INDEX_H_


#include <fragment_info.h>
#include <hash.h>
#include <couple.h>


int number_of_starters;
hash_t seeds; // hash table seed -> (fragment id, position)
hash_t seeds_count;
couple * seed_table;

p_fragment_info * index_starters_from_input_file (const int k, int nb_events_per_set, const int nb_fragment_per_event, const char input_only_upper, const int index_stride);           // read and store all starters presents in the pointed file. Index by seeds of length k all these starters.

void  free_seeds_index ();

#endif /* FRAGMENT_INDEX_H_ */
