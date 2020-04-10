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
#include<hash.h>
#include<stdint.h>
#include<assert.h>

class FragmentIndex{
public:
    hash_t seeds_count;
//    std::vector<std::pair <uint64_t, int >> seed_table;
    std::pair <uint64_t, int > * seed_table; // TODO change for a vector
    u_int64_t nb_coherent;
    u_int64_t nb_uncoherent;
    
    
    vector<Fragment*> all_predictions;
    

    void index_predictions (BankFasta inputBank, GlobalValues& gv);       // read and store all starters presents in the pointed file. Index by seeds of length k all these starters.
    void empty_coverage();
  
    
    
    FragmentIndex(const int numberOfIndexedSequences){
        seeds_count = hash_create_binarykey(); test_alloc(seeds_count);
        all_predictions.reserve(numberOfIndexedSequences);
    };
};



#endif /* FRAGMENT_INDEX_H_ */
