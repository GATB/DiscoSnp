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
 * fragment_index.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef FRAGMENT_INDEX_H_
#define FRAGMENT_INDEX_H_

#include<fragment_info.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<list.h>
#include<commons.h>
#include<couple.h>
#include<hash.h>
#include <stdint.h>
#include<assert.h>

class FragmentIndex{
public:
    hash_t seeds_count;
    couple * seed_table;
    u_int64_t nb_coherent;
    u_int64_t nb_uncoherent;
    
    
    vector<FragmentInfo*> all_predictions;
    
    void index_predictions (BankFasta inputBank, GlobalValues& gv);       // read and store all starters presents in the pointed file. Index by seeds of length k all these starters.
    void empty_coverage();
  
    
    
    FragmentIndex(const int numberOfIndexedSequences){
        seeds_count = hash_create_binarykey(100000); test_alloc(seeds_count); // todo  change to binary key (hash_t)AllocateHashTable(kmersize,1); //
        all_predictions.reserve(numberOfIndexedSequences);
        
    };
};



#endif /* FRAGMENT_INDEX_H_ */
