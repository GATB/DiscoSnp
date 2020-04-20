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
 * fragment_index.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef FRAGMENT_INDEX_H_
#define FRAGMENT_INDEX_H_

#include<fragment.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<commons.h>
#include<interface_xhash.h>
#include<stdint.h>
#include<assert.h>

class FragmentIndex{
public:
    
    // seed table: a seed is assigned to a fragment id (uint64_t) and the position of the seed on this fragment (int)
    // the seed table is a set of consecutive assigned seeds
    //    std::vector<std::pair <uint64_t, int >> seed_table;
    std::pair <uint64_t, int > * seed_table; // TODO change for a vector
    
    // the seeds_count enables to know the number of occurrences of each seed. Hence to know 1/ the totale number of occurrences of seekds (size of the seed table) and 2/ to know for each see where it starts in the seed table. 
    xhash seeds_count;

    u_int64_t nb_coherent;
    u_int64_t nb_uncoherent;
    
    
    vector<Fragment*> all_predictions;
    

    void index_predictions (BankFasta inputBank, GlobalValues& gv);       // read and store all starters presents in the pointed file. Index by seeds of length k all these starters.
    void empty_coverage();
  
    
    
    FragmentIndex(const int numberOfIndexedSequences){
//        seeds_count = hash_create_binarykey(); test_alloc(seeds_count);
        seeds_count = xhash_create_seed_index();
        all_predictions.reserve(numberOfIndexedSequences);
    };
};



#endif /* FRAGMENT_INDEX_H_ */
