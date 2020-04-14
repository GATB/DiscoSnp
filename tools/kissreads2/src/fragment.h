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
 * fragment.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef FRAGMENT_H_
#define FRAGMENT_H_

#include <interface_xhash.h>
#include <gatb/gatb_core.hpp>



class Fragment{
public:
    Sequence sequence;
    string upperCaseSequence;
	// fixed once at the beggining:
    unsigned int * SNP_positions; // If the fragment is a SNP, stores the positions of the SNPs in order to avoid to authorize errors at these positions. Coded on char, the SNP positions should not be longer than 255
    unsigned int nbOfSnps;                  // if zero: the sequence is generic or an indel. Else, number of predicted SNPs
    unsigned char * local_coverage;           //  number of reads covering this position can be a char, min coverage required is low
    bool * read_coherent; // for each read set: is the fragment read coherent?
    unsigned int * sum_qualities; // sum of the mapped qualities for each read set
    unsigned int * nb_mapped_qualities; // number of quality mapped for each read set. If there is a unique read, this is the number of mapped reads. In case of close SNPs, as a unique read may cover several SNPs, this can be bigger.
    
    int * number_mapped_reads;      //for every set of reads, number of reads starting (REPLACE reads_starting)
    
    
    
    Fragment(Sequence& seq, const int number_of_read_sets);
    
    
    ~Fragment(){
    }
    
    void set_read_coherent(int read_file_id, GlobalValues gv);

    
};


#endif /* FRAGMENT_H_ */
