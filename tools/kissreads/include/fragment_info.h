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
 * fragment_info.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef FRAGMENT_INFO_H_
#define FRAGMENT_INFO_H_

#include <hash.h>

#include "list.h"

//with CHARQUAL, use an unsigned char to store qual values. To avoid overflow, compute average value 'on the fly'
//#ifndef INTQUAL  // if compilation was NOT done with CFLAGS=-DINTQUAL, use the CHARQUAL option 
#define CHARQUAL 
//#endif

#ifndef CLASSICAL_SPANNING
#define KMER_SPANNING // ask more than only all position covered by at least min_coverage reads. Each kmer spanning each position should be covered by at least min_coverage reads.
#endif


typedef struct {
	// fixed once at the beggining:
	char * w;                     // the fragment
    char * SNP_positions; // If the fragment is a SNP, stores the positions of the SNPs in order to avoid to authorize errors at these positions. Coded on char, the SNP positions should not be longer than 255
    char * left_extension;
    char * right_extension;
    char * mapped_with_current_read; // For every read set (because of parrallelization): in run time, set to 1 if a match was seen with a current read. Avoid multiple matches of the same read on this starter.
    listint ** tested_pwis_with_current_read; // For every read set (because of parrallelization): in run time, all tested positions of the current read on this starter are stored. Avoids the computation redundances.
	char * comment;                // first line of the comment (fasta format)
    int nbOfSnps;                  // if zero: the sequence is generic or an indel. Else, number of predicted SNPs
	char * read_coherent;          // =for every set of reads, 1 if the fragment is detected as read coherent, else =0
	unsigned char ** read_coherent_positions;           //  number of reads covering this position can be a char, min coverage required is low


    unsigned int * sum_qualities; // sum of the mapped qualities for each read set
    unsigned int * nb_mapped_qualities; // number of quality mapped for each read set. If there is a unique read, this is the number of mapped reads. In case of close SNPs, as a unique read may cover several SNPs, this can be bigger.
    
	int * number_mapped_reads;      //for every set of reads, number of reads starting (REPLACE reads_starting)
} fragment_info, *p_fragment_info;


#endif /* FRAGMENT_INFO_H_ */
