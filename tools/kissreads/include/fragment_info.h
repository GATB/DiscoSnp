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
 * fragment_info.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef FRAGMENT_INFO_H_
#define FRAGMENT_INFO_H_

#include <hash.h>

//with CHARQUAL, use an unsigned char to store qual values. To avoid overflow, compute average value 'on the fly'
#ifndef INTQUAL  // if compilation was NOT done with CFLAGS=-DINTQUAL, use the CHARQUAL option 
#define CHARQUAL 
#endif

#ifndef CLASSICAL_SPANNING
#define KMER_SPANNING // ask more than only all position covered by at least min_coverage reads. Each kmer spanning each position should be covered by at least min_coverage reads.
#endif

//#define INPUT_FROM_KISSPLICE

typedef struct {
	// fixed once at the beggining:
	char * w;                     // the starter
//#ifdef GET_ONLY_UPPER_CHARS
    char * left_extension;
    char * right_extension;
//#endif
	char * comment;                // first line of the comment (fasta format)
	char * read_coherent;           // =for every set of reads, 1 if the fragment is detected as read coherent, else =0
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
