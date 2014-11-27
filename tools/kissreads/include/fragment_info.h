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
#ifndef INTQUAL  // if compilation was NOT done with CFLAGS=-DINTQUAL, use the CHARQUAL option 
#define CHARQUAL 
#endif

#ifndef CLASSICAL_SPANNING
#define KMER_SPANNING // ask more than only all position covered by at least min_coverage reads. Each kmer spanning each position should be covered by at least min_coverage reads.
#endif


typedef struct {
	// fixed once at the beggining:
	char * w;                     // the fragment
//#ifdef GET_ONLY_UPPER_CHARS
    char * left_extension;
    char * right_extension;
    char * mapped_with_current_read; // For every read set (because of parrallelization): in run time, set to 1 if a match was seen with a current read. Avoid multiple matches of the same read on this starter.
    listint ** tested_pwis_with_current_read; // For every read set (because of parrallelization): in run time, all tested positions of the current read on this starter are stored. Avoids the computation redundances.
//#endif
	char * comment;                // first line of the comment (fasta format)
    char isASNP;                      // if true (1), the sequence is a unique SNP. In this case we do not authorize any subsitution at this position.
	char * read_coherent;          // =for every set of reads, 1 if the fragment is detected as read coherent, else =0
	unsigned char ** read_coherent_positions;           //  number of reads covering this position can be a char, min coverage required is low
#ifdef INPUT_FROM_KISSPLICE
    char upperpath; // are we dealing with ASB path (==upper path) (else it's an AB path)
    
    int * nb_reads_fully_in_S;
    int * nb_reads_overlapping_AS;
    int * nb_reads_overlapping_SB;
    int * nb_reads_overlapping_both_AS_and_SB; // also used as junction AB for lower paths.
#endif
#ifdef CHARQUAL
	unsigned char ** sum_quality_per_position;	//for every set of reads, " " sum of the qualitites of reads mapping at that position if read coherent
#else
    int ** sum_quality_per_position;	
#endif
	int * number_mapped_reads;      //for every set of reads, number of reads starting (REPLACE reads_starting)
} fragment_info, *p_fragment_info;


#endif /* FRAGMENT_INFO_H_ */
