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
 * coherent_algorithm.c
 *
 *  Created on: 16 sept. 2010
 *      Author: pierre peterlongo
 */

#include<coherence_algorithm.h>
#include<commons.h>
#include <string.h>
#include <stdlib.h>


/**
 *  | pwi (may be negative)
 *     --------------  fragment
 *  *************      read
 *
 *  Tests if the overlapping part between read and fragment do not have more than <code>subst_allowed</code> substitions
 *  returns 1 if true between read and fragment, 0 else
 */
char read_coherent_generic(const int pwi, const char * fragment, const char * read, const int subst_allowed){
    
//    printf("read coherent, subst_allowed = %d\n", subst_allowed);
	int substitution_seen=0; // number of seen substitutions for now

	int pos_on_read, pos_on_fragment; // where to start



	/*
	 *  | pwi (negative)
	 *     --------------  fragment
	 *  *************      read
	 *     | we start here
	 */
	if(pwi<0) {
		pos_on_fragment=0;
		pos_on_read=-pwi;
	}

	/*
	 *        | pwi (positive)
	 *     --------------  fragment
	 *        *************      read
	 *        | we start here
	 */
	else{
		pos_on_fragment=pwi;
		pos_on_read=0;

	}


	// walk the read and the fragment together, detecting substitutions.
	// stop if the number of substitution is too high
	while(fragment[pos_on_fragment]!='\0' && read[pos_on_read]!='\0'){
		//if(fragment[pos_on_fragment]!=read[pos_on_read]) && fragment[pos_on_fragment]!='*') {// one subsitution
		if(fragment[pos_on_fragment]!=read[pos_on_read] && fragment[pos_on_fragment]!='*' && fragment[pos_on_fragment]!='?' && fragment[pos_on_fragment]!='N'){ // one subsitution
			substitution_seen++;
			if(substitution_seen>subst_allowed) break; // too much subsitutions
		}
		pos_on_fragment++;
		pos_on_read++;
	}
//    printf("%d subst seen\n", substitution_seen);
	if(substitution_seen<=subst_allowed) return 1;
	return 0;
}

/**
 *  | pwi (may be negative)
 *     --------------  fragment
 *  *************      read
 *
 *  Tests if the overlapping part between read and fragment do not have more than <code>subst_allowed</code> substitions
 * In case of SNPs, we need to avoid any substitution on the central fragment position (the one containing the SNP)
 * Thus in this function, we return 0 if any substitution occurs on this central position, whatever the number of substitution_seen
 *  returns 1 if true between read and fragment, 0 else
 */
char read_coherent_SNP(const int pwi, const char * fragment, const char * read, const int subst_allowed){
	int substitution_seen=0; // number of seen substitutions for now
    
	int pos_on_read, pos_on_fragment; // where to start
    
    
    
	/*
	 *  | pwi (negative)
	 *     --------------  fragment
	 *  *************      read
	 *     | we start here
	 */
	if(pwi<0) {
		pos_on_fragment=0;
		pos_on_read=-pwi;
	}
    
	/*
	 *        | pwi (positive)
	 *     --------------  fragment
	 *        *************      read
	 *        | we start here
	 */
	else{
		pos_on_fragment=pwi;
		pos_on_read=0;
	}
    
    const int snp_pos = strlen(fragment)/2;
	// walk the read and the fragment together, detecting substitutions.
	// stop if the number of substitution is too high
	while(fragment[pos_on_fragment]!='\0' && read[pos_on_read]!='\0'){
		if(fragment[pos_on_fragment]!=read[pos_on_read] &&
           fragment[pos_on_fragment]!='*' &&
           fragment[pos_on_fragment]!='?' &&
           fragment[pos_on_fragment]!='N'){ // one subsitution
			substitution_seen++;
			if(substitution_seen>subst_allowed) return 0; // too much subsitutions
            if(pos_on_fragment==snp_pos) return 0; // substition should not be on the snp
		}
		pos_on_fragment++;
		pos_on_read++;
	}
	return 1;
}
